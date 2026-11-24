package main

import (
	"bytes"
	"encoding/json"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

func TestViewDoesNotDownloadOrResolveMAFFT(t *testing.T) {
	a := &App{token: "test", origin: "http://127.0.0.1:12345", last: time.Now(), quit: make(chan struct{})}
	req := httptest.NewRequest("GET", "http://127.0.0.1:12345/api/status", nil)
	req.Header.Set("X-Review-Token", "test")
	w := httptest.NewRecorder()
	a.handle(w, req)
	if w.Code != 200 || a.mafft != "" || a.temporary != "" {
		t.Fatal("Viewing tried to prepare a runtime")
	}
	var v map[string]any
	json.Unmarshal(w.Body.Bytes(), &v)
	if v["mafft_loaded"] != false {
		t.Fatal(v)
	}
}
func TestShutdownRefusesAnActiveRecalculation(t *testing.T) {
	a := &App{token: "test", origin: "http://127.0.0.1:12345", busy: true, last: time.Now(), quit: make(chan struct{})}
	req := httptest.NewRequest("POST", "http://127.0.0.1:12345/api/shutdown", bytes.NewBufferString("{}"))
	req.Header.Set("X-Review-Token", "test")
	w := httptest.NewRecorder()
	a.handle(w, req)
	if w.Code != 409 || a.closing {
		t.Fatal("Active operation could be cleaned")
	}
	select {
	case <-a.quit:
		t.Fatal("Closed while busy")
	default:
	}
}
func TestCodonsPreserveNativeSequenceAndDoNotInventQueryCDS(t *testing.T) {
	selected := map[string]Sequence{"model": {CDS: "ATGGCTTAA"}, "query": {}}
	out, e := codingAlignment(map[string]string{"model": "M-A", "query": "MVA"}, selected)
	if e != nil || out["model"] != "ATG---GCT" {
		t.Fatal(out, e)
	}
	if _, ok := out["query"]; ok {
		t.Fatal("Invented query CDS")
	}
	_, e = codingAlignment(map[string]string{"model": "M-V"}, selected)
	if e == nil {
		t.Fatal("Accepted a protein/CDS mismatch")
	}
}
func TestOnlyOwnedStaleDirectoriesAreRemoved(t *testing.T) {
	p, e := markedTemp("gcf-review-align-")
	if e != nil {
		t.Fatal(e)
	}
	defer os.RemoveAll(p)
	cleanupStaleTemps()
	if _, e = os.Stat(p); e != nil {
		t.Fatal("Removed a live process directory")
	}
	// The archive helper refuses traversal before creating any file.
	root := t.TempDir()
	if copyFile(root, "../outside", 0600, strings.NewReader("x")) == nil {
		t.Fatal("Archive traversal allowed")
	}
	if _, e = os.Stat(filepath.Join(filepath.Dir(root), "outside")); e == nil {
		t.Fatal("Escaped root")
	}
}
