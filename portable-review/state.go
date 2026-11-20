package main

import (
	"bytes"
	_ "embed"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"
)

//go:embed state_bootstrap.js
var stateBootstrap string

type browserState struct {
	Schema     string            `json:"schema"`
	DatasetSHA string            `json:"dataset_sha256"`
	Values     map[string]string `json:"values"`
	Saved      string            `json:"saved_at"`
}

func (a *App) stateKeyAllowed(k string) bool {
	for _, prefix := range []string{"gcf-realign-selection-" + a.data.ReviewID, "gcf-review-workspace-2-" + a.data.ReviewID, "gcf-one-step-undo-" + a.data.ReviewID, "gcf-outgroup-" + a.data.ReviewID, "gcf-alignment-first-notes-" + a.data.ReviewID} {
		if strings.HasPrefix(k, prefix) {
			return true
		}
	}
	if a.data.SourceReviewID != "" {
		for _, prefix := range []string{"gcf-stage2-" + a.data.SourceReviewID, "gcf-stage2-views-" + a.data.SourceReviewID, "gcf-synteny-final-" + a.data.SourceReviewID, "gcf-synteny-reference-v1-" + a.data.SourceReviewID, "gcf-stop4-final-" + a.data.SourceReviewID} {
			if strings.HasPrefix(k, prefix) {
				return true
			}
		}
	}
	return false
}
func (a *App) readState() browserState {
	v := browserState{Schema: "gcf_portable_session_v2", DatasetSHA: digest(a.raw), Values: map[string]string{}}
	raw, e := os.ReadFile(filepath.Join(a.root, ".review", "session-state.json"))
	if e != nil {
		return v
	}
	var old browserState
	if json.Unmarshal(raw, &old) == nil && old.Schema == v.Schema && old.DatasetSHA == v.DatasetSHA {
		v.Saved = old.Saved
		for k, s := range old.Values {
			if a.stateKeyAllowed(k) {
				v.Values[k] = s
			}
		}
	}
	return v
}
func (a *App) handleState(w http.ResponseWriter, r *http.Request) {
	a.mu.Lock()
	defer a.mu.Unlock()
	if r.Method == "GET" {
		writeJSON(w, a.readState())
		return
	}
	if r.Method != "POST" {
		http.Error(w, "POST required", 405)
		return
	}
	raw, e := io.ReadAll(io.LimitReader(r.Body, 16<<20+1))
	if e != nil || len(raw) > 16<<20 {
		fail(w, errors.New("Review state exceeds 16 MiB"))
		return
	}
	var state browserState
	if json.Unmarshal(raw, &state) != nil || state.Schema != "gcf_portable_session_v2" || state.DatasetSHA != digest(a.raw) {
		fail(w, errors.New("Review state belongs to another dataset"))
		return
	}
	for key := range state.Values {
		if !a.stateKeyAllowed(key) {
			fail(w, errors.New("Unexpected review-state key"))
			return
		}
	}
	state.Saved = time.Now().UTC().Format(time.RFC3339Nano)
	data, _ := json.Marshal(state)
	root := filepath.Join(a.root, ".review")
	if e = os.MkdirAll(root, 0700); e != nil {
		fail(w, e)
		return
	}
	f, e := os.CreateTemp(root, ".state-")
	if e != nil {
		fail(w, e)
		return
	}
	name := f.Name()
	defer os.Remove(name)
	_, e = f.Write(data)
	if close := f.Close(); e == nil {
		e = close
	}
	if e == nil {
		e = os.Rename(name, filepath.Join(root, "session-state.json"))
	}
	if e != nil {
		fail(w, e)
		return
	}
	writeJSON(w, map[string]any{"saved": true, "saved_at": state.Saved})
}
func (a *App) serveReview(w http.ResponseWriter, r *http.Request) {
	raw, e := os.ReadFile(filepath.Join(a.root, "report.html"))
	if e != nil {
		http.Error(w, "Extract the full report package", 404)
		return
	}
	state := a.readState()
	config := map[string]any{"state": state, "review_id": a.data.ReviewID, "source_review_id": a.data.SourceReviewID, "helper_source_sha256": buildHash}
	encoded, _ := json.Marshal(config)
	script := []byte("<script id=\"gcfPortableState\" type=\"application/json\">" + string(encoded) + "</script><script>" + stateBootstrap + "</script>")
	raw = bytes.Replace(raw, []byte("</head>"), append(script, []byte("</head>")...), 1)
	w.Header().Set("Content-Type", "text/html; charset=utf-8")
	w.Header().Set("Content-Length", strconv.Itoa(len(raw)))
	if r.Method != "HEAD" {
		w.Write(raw)
	}
}
func cleanupIQTREEStaleCaches() {
	packages := map[string]Package{}
	json.Unmarshal(iqtreeManifest, &packages)
	for _, p := range packages {
		root := filepath.Join(os.TempDir(), "gcf-review-iqtree-"+p.SHA256[:16])
		raw, e := os.ReadFile(filepath.Join(root, ".gcf-cache-owner"))
		if e != nil || string(raw) != "gcf-review-cache-v1" {
			continue
		}
		lock := root + ".lock"
		if os.Mkdir(lock, 0700) != nil {
			continue
		}
		os.WriteFile(filepath.Join(lock, "pid"), []byte(fmt.Sprint(os.Getpid())), 0600)
		files, _ := os.ReadDir(filepath.Join(root, "leases"))
		live := false
		for _, f := range files {
			pid, _ := strconv.Atoi(f.Name())
			if pid > 0 && alive(pid) {
				live = true
			}
		}
		if !live {
			os.RemoveAll(root)
		}
		os.RemoveAll(lock)
	}
}
