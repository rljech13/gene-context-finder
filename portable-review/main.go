// Gene-context local review runner. No network access before an explicit align.
package main

import (
	"archive/tar"
	"archive/zip"
	"bytes"
	"compress/gzip"
	"context"
	"crypto/rand"
	"crypto/sha256"
	_ "embed"
	"encoding/hex"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"net"
	"net/http"
	"net/url"
	"os"
	"os/exec"
	"os/signal"
	"path/filepath"
	"regexp"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"syscall"
	"time"
)

//go:embed mafft_manifest.json
var manifestBytes []byte

type Package struct {
	URL     string
	SHA256  string
	Archive string
	Entry   string
}
type Sequence struct {
	ID        string
	Protein   string
	CDS       string
	Group     string
	Label     string
	Reference bool
	Metadata  json.RawMessage
}
type Dataset struct {
	SourceReviewID string `json:"source_review_id"`
	Schema         string
	ReviewID       string `json:"review_id"`
	Stage          int
	Sequences      []Sequence
	History        json.RawMessage
}
type Request struct {
	Notes             string          `json:"notes"`
	Phase             string          `json:"phase"`
	AlignmentResultID string          `json:"alignment_result_id"`
	OutgroupIDs       []string        `json:"outgroup_ids"`
	TrimMode          string          `json:"trim_mode"`
	RootBasis         string          `json:"root_basis"`
	Decisions         json.RawMessage `json:"decisions"`
	IDs               []string        `json:"ids"`
	Group             string          `json:"group"`
}
type Result struct {
	Phase       string            `json:"phase,omitempty"`
	AlignmentID string            `json:"alignment_id,omitempty"`
	Phylogeny   map[string]any    `json:"phylogeny,omitempty"`
	Provenance  map[string]any    `json:"provenance,omitempty"`
	ID          string            `json:"id"`
	ReviewID    string            `json:"review_id"`
	InputSHA    string            `json:"input_sha256"`
	IDs         []string          `json:"ids"`
	Group       string            `json:"group"`
	Protein     map[string]string `json:"protein"`
	CDS         map[string]string `json:"cds"`
	MAFFT       string            `json:"mafft_version"`
	Created     string            `json:"created"`
}
type App struct {
	root, token, origin, temporary, mafft, version string
	iqtree, iqtreeVersion, iqtreeTemporary         string
	data                                           Dataset
	raw                                            []byte
	byID                                           map[string]Sequence
	mu                                             sync.Mutex
	busy                                           bool
	closing                                        bool
	last                                           time.Time
	quit                                           chan struct{}
	stop                                           sync.Once
}

func digest(b []byte) string { s := sha256.Sum256(b); return hex.EncodeToString(s[:]) }
func writeJSON(w http.ResponseWriter, v any) {
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(v)
}
func fail(w http.ResponseWriter, e error) { http.Error(w, e.Error(), http.StatusBadRequest) }
func (a *App) touch()                     { a.mu.Lock(); a.last = time.Now(); a.mu.Unlock() }
func (a *App) authorised(r *http.Request) bool {
	return r.Host == strings.TrimPrefix(a.origin, "http://") && r.Header.Get("X-Review-Token") == a.token && (r.Header.Get("Origin") == "" || r.Header.Get("Origin") == a.origin)
}
func (a *App) shutdown() { a.stop.Do(func() { close(a.quit) }) }
func fasta(records map[string]string, ids []string) []byte {
	var b strings.Builder
	for _, id := range ids {
		s, ok := records[id]
		if !ok {
			continue
		}
		fmt.Fprintf(&b, ">%s\n", id)
		for len(s) > 80 {
			b.WriteString(s[:80] + "\n")
			s = s[80:]
		}
		b.WriteString(s + "\n")
	}
	return []byte(b.String())
}
func parseFASTA(raw []byte) (map[string]string, error) {
	out := map[string]string{}
	id := ""
	for _, l := range strings.Split(strings.ReplaceAll(string(raw), "\r", ""), "\n") {
		l = strings.TrimSpace(l)
		if l == "" {
			continue
		}
		if strings.HasPrefix(l, ">") {
			words := strings.Fields(l[1:])
			if len(words) == 0 {
				return nil, errors.New("Missing FASTA ID")
			}
			id = words[0]
			if _, ok := out[id]; ok {
				return nil, errors.New("Duplicate FASTA ID")
			}
			out[id] = ""
		} else {
			if id == "" {
				return nil, errors.New("Invalid FASTA")
			}
			out[id] += strings.ToUpper(l)
		}
	}
	return out, nil
}
func safePath(root, name string) (string, error) {
	name = strings.ReplaceAll(name, "\\", "/")
	p := filepath.Join(root, filepath.FromSlash(name))
	rel, e := filepath.Rel(root, p)
	if e != nil || filepath.IsAbs(name) || rel == ".." || strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
		return "", errors.New("Unsafe archive path")
	}
	return p, nil
}
func copyFile(root, name string, mode os.FileMode, src io.Reader) error {
	p, e := safePath(root, name)
	if e != nil {
		return e
	}
	if e = os.MkdirAll(filepath.Dir(p), 0700); e != nil {
		return e
	}
	f, e := os.OpenFile(p, os.O_CREATE|os.O_TRUNC|os.O_WRONLY, mode.Perm()|0600)
	if e != nil {
		return e
	}
	_, e = io.Copy(f, io.LimitReader(src, 512<<20))
	ce := f.Close()
	if e != nil {
		return e
	}
	return ce
}
func unpack(archive, root string) error {
	if strings.HasSuffix(archive, ".zip") {
		z, e := zip.OpenReader(archive)
		if e != nil {
			return e
		}
		defer z.Close()
		for _, f := range z.File {
			if f.Mode()&os.ModeSymlink != 0 {
				return errors.New("Archive symlink rejected")
			}
			if f.FileInfo().IsDir() {
				continue
			}
			r, e := f.Open()
			if e != nil {
				return e
			}
			e = copyFile(root, f.Name, f.Mode(), r)
			r.Close()
			if e != nil {
				return e
			}
		}
		return nil
	}
	f, e := os.Open(archive)
	if e != nil {
		return e
	}
	defer f.Close()
	gz, e := gzip.NewReader(f)
	if e != nil {
		return e
	}
	defer gz.Close()
	t := tar.NewReader(gz)
	for {
		h, e := t.Next()
		if e == io.EOF {
			return nil
		}
		if e != nil {
			return e
		}
		if h.Typeflag == tar.TypeDir {
			continue
		}
		if h.Typeflag != tar.TypeReg {
			return errors.New("Non-regular archive member rejected")
		}
		if e = copyFile(root, h.Name, os.FileMode(h.Mode), t); e != nil {
			return e
		}
	}
}
func command(path string, args ...string) *exec.Cmd {
	if runtime.GOOS == "windows" && strings.HasSuffix(strings.ToLower(path), ".bat") {
		return exec.Command("cmd.exe", append([]string{"/D", "/C", path}, args...)...)
	}
	return exec.Command(path, args...)
}
func probe(path string) (string, error) {
	c := command(path, "--version")
	var out bytes.Buffer
	c.Stdout = &out
	c.Stderr = &out
	e := c.Run()
	if e != nil {
		return "", e
	}
	v := strings.TrimSpace(out.String())
	matches := regexp.MustCompile(`(?m)\bv(\d+)\.(\d+)`).FindStringSubmatch(v)
	if len(matches) != 3 {
		return "", errors.New("Unrecognised MAFFT version")
	}
	major, _ := strconv.Atoi(matches[1])
	minor, _ := strconv.Atoi(matches[2])
	if major < 7 || (major == 7 && minor < 487) {
		return "", errors.New("MAFFT 7.487 or newer is required")
	}
	return v, nil
}
func alive(pid int) bool {
	if pid == os.Getpid() {
		return true
	}
	if runtime.GOOS == "windows" {
		b, e := exec.Command("tasklist", "/FI", "PID eq "+strconv.Itoa(pid), "/NH").Output()
		return e == nil && strings.Contains(string(b), strconv.Itoa(pid))
	}
	p, e := os.FindProcess(pid)
	return e == nil && p.Signal(syscall.Signal(0)) == nil
}

// A shared cache is leased by live launcher processes. Never clean another
// report's package, or replace a live package. A crashed lease is checked by PID.
func cacheLock(root string) (func(), error) {
	lock := root + ".lock"
	for i := 0; i < 600; i++ {
		if os.Mkdir(lock, 0700) == nil {
			os.WriteFile(filepath.Join(lock, "pid"), []byte(strconv.Itoa(os.Getpid())), 0600)
			return func() { os.RemoveAll(lock) }, nil
		}
		b, _ := os.ReadFile(filepath.Join(lock, "pid"))
		pid, _ := strconv.Atoi(string(b))
		if pid > 0 && !alive(pid) {
			os.RemoveAll(lock)
		}
		time.Sleep(100 * time.Millisecond)
	}
	return nil, errors.New("Another launcher is preparing MAFFT; try again")
}
func markedTemp(prefix string) (string, error) {
	p, e := os.MkdirTemp(os.TempDir(), prefix)
	if e == nil {
		b, _ := json.Marshal(map[string]any{"schema": "gcf-review-temp-v1", "pid": os.Getpid()})
		e = os.WriteFile(filepath.Join(p, ".gcf-owner.json"), b, 0600)
	}
	return p, e
}
func cleanupStaleTemps() {
	entries, _ := os.ReadDir(os.TempDir())
	for _, f := range entries {
		if !f.IsDir() || !(strings.HasPrefix(f.Name(), "gcf-review-align-") || strings.HasPrefix(f.Name(), "gcf-review-download-")) {
			continue
		}
		p := filepath.Join(os.TempDir(), f.Name())
		var owner struct {
			Schema string
			PID    int
		}
		b, e := os.ReadFile(filepath.Join(p, ".gcf-owner.json"))
		if e == nil && json.Unmarshal(b, &owner) == nil && owner.Schema == "gcf-review-temp-v1" && owner.PID > 0 && !alive(owner.PID) {
			os.RemoveAll(p)
		}
	}
	pkgs := map[string]Package{}
	json.Unmarshal(manifestBytes, &pkgs)
	for _, p := range pkgs {
		root := filepath.Join(os.TempDir(), "gcf-review-mafft-"+p.SHA256[:16])
		if _, e := os.Stat(filepath.Join(root, ".gcf-cache-owner")); e != nil {
			continue
		}
		lock := root + ".lock"
		if os.Mkdir(lock, 0700) != nil {
			continue
		}
		os.WriteFile(filepath.Join(lock, "pid"), []byte(strconv.Itoa(os.Getpid())), 0600)
		live := false
		files, _ := os.ReadDir(filepath.Join(root, "leases"))
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
func (a *App) setMAFFT(path, version string) {
	a.mu.Lock()
	a.mafft = path
	a.version = version
	a.mu.Unlock()
}
func (a *App) ensureMAFFT() error {
	if a.mafft != "" {
		return nil
	}
	if p, e := exec.LookPath("mafft"); e == nil && os.Getenv("GCF_REVIEW_FORCE_PORTABLE") != "1" {
		if v, e := probe(p); e == nil {
			a.setMAFFT(p, v)
			return nil
		}
	}
	pkgs := map[string]Package{}
	json.Unmarshal(manifestBytes, &pkgs)
	p, ok := pkgs[runtime.GOOS]
	if !ok || runtime.GOOS == "linux" && runtime.GOARCH != "amd64" {
		return errors.New("Install a native MAFFT for this architecture; the report and local runner remain usable")
	}
	root := filepath.Join(os.TempDir(), "gcf-review-mafft-"+p.SHA256[:16])
	unlock, e := cacheLock(root)
	if e != nil {
		return e
	}
	defer unlock()
	os.MkdirAll(root, 0700)
	os.WriteFile(filepath.Join(root, ".gcf-cache-owner"), []byte("gcf-review-cache-v1"), 0600)
	a.temporary = root
	entry := filepath.Join(root, p.Entry)
	if _, e = os.Stat(entry); e != nil {
		work, e := markedTemp("gcf-review-download-")
		if e != nil {
			return e
		}
		defer os.RemoveAll(work)
		resp, e := (&http.Client{Timeout: 5 * time.Minute}).Get(p.URL)
		if e != nil {
			return fmt.Errorf("MAFFT download: %w", e)
		}
		defer resp.Body.Close()
		if resp.StatusCode != 200 {
			return fmt.Errorf("MAFFT download HTTP %d", resp.StatusCode)
		}
		file := filepath.Join(work, p.Archive)
		f, e := os.Create(file)
		if e != nil {
			return e
		}
		_, e = io.Copy(f, io.LimitReader(resp.Body, 256<<20))
		f.Close()
		if e != nil {
			return e
		}
		b, e := os.ReadFile(file)
		if e != nil {
			return e
		}
		if digest(b) != p.SHA256 {
			return errors.New("MAFFT archive checksum mismatch")
		}
		extracted := filepath.Join(work, "unpacked")
		os.Mkdir(extracted, 0700)
		if e = unpack(file, extracted); e != nil {
			return e
		}
		for _, d := range []string{"mafft-mac", "mafft-linux64", "mafft-win"} {
			src := filepath.Join(extracted, d)
			if _, e = os.Stat(src); e == nil {
				if e = os.Rename(src, filepath.Join(root, d)); e != nil {
					return e
				}
			}
		}
	}
	// Some portable archives record Unix executable modes incompletely.
	if runtime.GOOS != "windows" {
		filepath.Walk(root, func(p string, i os.FileInfo, e error) error {
			if e == nil && !i.IsDir() {
				os.Chmod(p, i.Mode()|0700)
			}
			return nil
		})
	}
	v, e := probe(entry)
	if e != nil {
		return fmt.Errorf("Portable MAFFT cannot start: %w", e)
	}
	os.MkdirAll(filepath.Join(root, "leases"), 0700)
	if e = os.WriteFile(filepath.Join(root, "leases", strconv.Itoa(os.Getpid())), []byte(a.token), 0600); e != nil {
		return e
	}
	a.temporary = root
	a.setMAFFT(entry, v)
	return nil
}
func (a *App) cleanup() {
	a.cleanupIQTREE()
	if a.temporary == "" {
		return
	}
	root := a.temporary
	unlock, e := cacheLock(root)
	if e != nil {
		return
	}
	defer unlock()
	os.Remove(filepath.Join(root, "leases", strconv.Itoa(os.Getpid())))
	files, _ := os.ReadDir(filepath.Join(root, "leases"))
	for _, f := range files {
		pid, _ := strconv.Atoi(f.Name())
		if alive(pid) {
			return
		}
		os.Remove(filepath.Join(root, "leases", f.Name()))
	}
	os.RemoveAll(root)
}

const bases = "TCAG"
const aaTable = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

func amino(c string) byte {
	n := 0
	for _, b := range c {
		i := strings.IndexRune(bases, b)
		if i < 0 {
			return 'X'
		}
		n = n*4 + i
	}
	return aaTable[n]
}
func codingAlignment(protein map[string]string, selected map[string]Sequence) (map[string]string, error) {
	out := map[string]string{}
	for id, aln := range protein {
		s := selected[id]
		if s.CDS == "" {
			continue
		}
		dna := strings.ToUpper(s.CDS)
		if len(dna)%3 != 0 {
			return nil, fmt.Errorf("%s: CDS length is not divisible by 3", id)
		}
		if len(dna) >= 3 && amino(dna[len(dna)-3:]) == '*' {
			dna = dna[:len(dna)-3]
		}
		var b strings.Builder
		pos := 0
		for _, aa := range aln {
			if aa == '-' {
				b.WriteString("---")
				continue
			}
			if pos+3 > len(dna) {
				return nil, fmt.Errorf("%s: CDS is shorter than the protein", id)
			}
			codon := dna[pos : pos+3]
			translated := amino(codon)
			if translated != byte(aa) && aa != 'X' {
				return nil, fmt.Errorf("%s: CDS/protein mismatch at codon %d", id, pos/3+1)
			}
			b.WriteString(codon)
			pos += 3
		}
		if pos != len(dna) {
			return nil, fmt.Errorf("%s: unused CDS bases", id)
		}
		out[id] = b.String()
	}
	return out, nil
}
func (a *App) align(q Request) (Result, error) {
	if q.Phase == "tree" {
		return a.tree(q)
	}
	if q.Phase != "" && q.Phase != "alignment" {
		return Result{}, errors.New("Unknown review operation")
	}
	r := Result{Phase: "alignment", ReviewID: a.data.ReviewID, Group: q.Group, Created: time.Now().UTC().Format(time.RFC3339)}
	selected := map[string]Sequence{}
	proteins := map[string]string{}
	for _, id := range q.IDs {
		s, ok := a.byID[id]
		if !ok || s.Group != q.Group {
			return r, errors.New("Unknown sequence or different gene group")
		}
		if _, ok := selected[id]; ok {
			return r, errors.New("Duplicate model")
		}
		selected[id] = s
		proteins[id] = s.Protein
	}
	if len(selected) < 2 || len(selected) > 2000 {
		return r, errors.New("Select between 2 and 2000 sequences")
	}
	if e := a.ensureMAFFT(); e != nil {
		return r, e
	}
	r.MAFFT = a.version
	r.IDs = append([]string{}, q.IDs...)
	sort.Strings(r.IDs)
	r.Provenance = map[string]any{"method": "portable_review23_v2", "helper_source_sha256": buildHash,
		"dataset_sha256": digest(a.raw), "model_choice_notes": q.Notes}
	signature, _ := json.Marshal([]any{r.Provenance, digest(q.Decisions), r.IDs, a.version, []string{"--auto", "--amino", "--thread", "2"}})
	r.InputSHA = digest(signature)
	r.ID = r.InputSHA[:20]
	dest := filepath.Join(a.root, "review_results", r.ID)
	if _, e := os.Stat(filepath.Join(dest, "result.json")); e == nil {
		return a.verifiedResult(r.ID)
	}
	work, e := markedTemp("gcf-review-align-")
	if e != nil {
		return r, e
	}
	defer os.RemoveAll(work)
	input := filepath.Join(work, "proteins.faa")
	os.WriteFile(input, fasta(proteins, r.IDs), 0600)
	c := command(a.mafft, "--auto", "--amino", "--thread", "2", input)
	var output, stderr bytes.Buffer
	c.Stdout = &output
	c.Stderr = &stderr
	// Never clean the runtime or its input while MAFFT is still running.
	e = c.Start()
	if e == nil {
		owner, _ := json.Marshal(map[string]any{"schema": "gcf-review-temp-v1", "pid": c.Process.Pid})
		os.WriteFile(filepath.Join(work, ".gcf-owner.json"), owner, 0600)
		lease := ""
		if a.temporary != "" {
			lease = filepath.Join(a.temporary, "leases", strconv.Itoa(c.Process.Pid))
			os.WriteFile(lease, []byte(a.token), 0600)
		}
		e = c.Wait()
		if lease != "" {
			os.Remove(lease)
		}
	}
	if e != nil {
		return r, fmt.Errorf("MAFFT failed: %v\n%s", e, stderr.String())
	}
	aln, e := parseFASTA(output.Bytes())
	if e != nil {
		return r, e
	}
	if len(aln) != len(selected) {
		return r, errors.New("MAFFT returned a different sequence set")
	}
	n := -1
	for id, s := range selected {
		v, ok := aln[id]
		if !ok || strings.ReplaceAll(v, "-", "") != s.Protein {
			return r, errors.New("MAFFT changed a source sequence")
		}
		if n >= 0 && len(v) != n {
			return r, errors.New("Unequal MSA row lengths")
		}
		n = len(v)
	}
	r.Protein = aln
	r.CDS, e = codingAlignment(aln, selected)
	if e != nil {
		return r, e
	}
	parent := filepath.Join(a.root, "review_results")
	if e = os.MkdirAll(parent, 0700); e != nil {
		return r, e
	}
	pending, e := os.MkdirTemp(parent, ".pending-")
	if e != nil {
		return r, e
	}
	defer os.RemoveAll(pending)
	fullCDS := map[string]string{}
	for id, s := range selected {
		if s.CDS != "" {
			fullCDS[id] = s.CDS
		}
	}
	files := map[string][]byte{"proteins.aligned.faa": fasta(r.Protein, r.IDs), "cds.codon.aligned.fna": fasta(r.CDS, r.IDs), "proteins.native.faa": fasta(proteins, r.IDs), "cds.native.fna": fasta(fullCDS, r.IDs), "mafft.log": []byte(stderr.String()), "source_history.json": a.data.History}
	b, _ := json.MarshalIndent(r, "", "  ")
	files["result.json"] = b
	if len(q.Decisions) > 0 {
		var d map[string]any
		if json.Unmarshal(q.Decisions, &d) != nil || d["review_id"] != a.data.ReviewID {
			return r, errors.New("Imported decisions belong to another review")
		}
		files["imported_decisions.json"] = q.Decisions
	}
	selection, _ := json.MarshalIndent(selected, "", "  ")
	files["selected_models.json"] = selection
	seals := map[string]string{}
	for name, data := range files {
		if e = os.WriteFile(filepath.Join(pending, name), data, 0600); e != nil {
			return r, e
		}
		seals[name] = digest(data)
	}
	seal, _ := json.MarshalIndent(seals, "", "  ")
	os.WriteFile(filepath.Join(pending, "checksums.json"), seal, 0600)
	if e = os.Rename(pending, dest); e != nil {
		return r, e
	}
	return r, nil
}
func (a *App) handle(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Cache-Control", "no-store")
	w.Header().Set("X-Content-Type-Options", "nosniff")
	if strings.HasPrefix(r.URL.Path, "/api/") {
		if !a.authorised(r) {
			http.Error(w, "Local session token required", 403)
			return
		}
		a.touch()
		switch r.URL.Path {
		case "/api/state":
			a.handleState(w, r)
		case "/api/status":
			a.mu.Lock()
			loaded := a.mafft != ""
			a.mu.Unlock()
			writeJSON(w, map[string]any{"ready": true, "review_id": a.data.ReviewID, "mafft_loaded": loaded, "trimal_supported": true, "scaled_tree_views": true})
		case "/api/data":
			writeJSON(w, a.data)
		case "/api/align":
			if r.Method != "POST" {
				http.Error(w, "POST required", 405)
				return
			}
			a.mu.Lock()
			if a.busy || a.closing {
				a.mu.Unlock()
				http.Error(w, "A recalculation is already running", 409)
				return
			}
			a.busy = true
			a.mu.Unlock()
			defer func() { a.mu.Lock(); a.busy = false; a.last = time.Now(); a.mu.Unlock() }()
			var q Request
			if e := json.NewDecoder(io.LimitReader(r.Body, 1<<20)).Decode(&q); e != nil {
				fail(w, e)
				return
			}
			v, e := a.align(q)
			if e != nil {
				fail(w, e)
				return
			}
			writeJSON(w, v)
		case "/api/checkpoint":
			if r.Method != http.MethodGet {
				http.Error(w, "GET required", http.StatusMethodNotAllowed)
				return
			}
			checkpoint, e := a.verifiedResult(r.URL.Query().Get("id"))
			if e != nil {
				fail(w, e)
				return
			}
			if e := scaledTreeViews(&checkpoint); e != nil {
				fail(w, e)
				return
			}
			writeJSON(w, checkpoint)
		case "/api/export":
			if r.Method != "POST" {
				http.Error(w, "POST required", 405)
				return
			}
			var v struct{ ID string }
			json.NewDecoder(io.LimitReader(r.Body, 1024)).Decode(&v)
			if !validResultID(v.ID) {
				fail(w, errors.New("Invalid result ID"))
				return
			}
			dir := filepath.Join(a.root, "review_results", v.ID)
			b, e := os.ReadFile(filepath.Join(dir, "checksums.json"))
			if e != nil {
				fail(w, e)
				return
			}
			var seal map[string]string
			json.Unmarshal(b, &seal)
			var buf bytes.Buffer
			z := zip.NewWriter(&buf)
			for name, h := range seal {
				safe, e := safePath(dir, name)
				if e != nil {
					fail(w, e)
					return
				}
				data, e := os.ReadFile(safe)
				if e != nil || digest(data) != h {
					z.Close()
					fail(w, errors.New("Saved result changed"))
					return
				}
				f, _ := z.Create(name)
				f.Write(data)
			}
			f, _ := z.Create("checksums.json")
			f.Write(b)
			z.Close()
			w.Header().Set("Content-Type", "application/zip")
			w.Write(buf.Bytes())
		case "/api/shutdown":
			if r.Method != "POST" {
				http.Error(w, "POST required", 405)
				return
			}
			a.mu.Lock()
			busy := a.busy
			if !busy {
				a.closing = true
			}
			a.mu.Unlock()
			if busy {
				http.Error(w, "Wait for the active alignment to finish", 409)
				return
			}
			writeJSON(w, map[string]bool{"saved": true})
			go func() { time.Sleep(200 * time.Millisecond); a.shutdown() }()
		default:
			http.NotFound(w, r)
		}
		return
	}
	if r.Host != strings.TrimPrefix(a.origin, "http://") {
		http.Error(w, "Invalid host", 403)
		return
	}
	if r.URL.Path == "/" || r.URL.Path == "/report.html" {
		a.serveReview(w, r)
		return
	}
	// No directory listings and no executable/script downloads from the helper.
	name := strings.TrimPrefix(r.URL.Path, "/")
	if strings.HasPrefix(name, ".review") || strings.HasSuffix(name, ".exe") || strings.HasSuffix(name, ".command") || strings.HasSuffix(name, ".sh") {
		http.NotFound(w, r)
		return
	}
	p, e := safePath(a.root, name)
	if e != nil {
		http.NotFound(w, r)
		return
	}
	resolved, re := filepath.EvalSymlinks(p)
	rel, pe := filepath.Rel(a.root, resolved)
	if re != nil || pe != nil || rel == ".." || strings.HasPrefix(rel, ".."+string(filepath.Separator)) {
		http.NotFound(w, r)
		return
	}
	info, e := os.Stat(p)
	if e != nil || info.IsDir() {
		http.NotFound(w, r)
		return
	}
	http.ServeFile(w, r, p)
}
func main() {
	root := flag.String("root", "", "Report folder")
	noBrowser := flag.Bool("no-browser", false, "Do not open a browser")
	flag.Parse()
	if *root == "" {
		exe, _ := os.Executable()
		*root = filepath.Dir(exe)
		if filepath.Base(*root) == "bin" {
			*root = filepath.Dir(filepath.Dir(*root))
		}
	}
	cleanupStaleTemps()
	cleanupIQTREEStaleCaches()
	absolute, _ := filepath.Abs(*root)
	raw, e := os.ReadFile(filepath.Join(absolute, "review_data.json"))
	if e != nil {
		fmt.Fprintln(os.Stderr, "Extract the complete review package, then start its launcher:", e)
		os.Exit(1)
	}
	a := &App{root: absolute, raw: raw, byID: map[string]Sequence{}, last: time.Now(), quit: make(chan struct{})}
	if e = json.Unmarshal(raw, &a.data); e != nil {
		panic(e)
	}
	for _, s := range a.data.Sequences {
		if s.ID == "" || s.Protein == "" {
			panic("Missing model sequence")
		}
		a.byID[s.ID] = s
	}
	token := make([]byte, 24)
	rand.Read(token)
	a.token = hex.EncodeToString(token)
	ln, e := net.Listen("tcp", "127.0.0.1:0")
	if e != nil {
		panic(e)
	}
	a.origin = "http://" + ln.Addr().String()
	address := a.origin + "/?review_token=" + url.QueryEscape(a.token)
	fmt.Println("Open:", address)
	fmt.Println("MAFFT is downloaded only after Recalculate. Use Finish and clean up in the report.")
	server := &http.Server{Handler: http.HandlerFunc(a.handle), ReadHeaderTimeout: 10 * time.Second}
	go server.Serve(ln)
	if !*noBrowser {
		var c *exec.Cmd
		switch runtime.GOOS {
		case "darwin":
			c = exec.Command("open", address)
		case "windows":
			c = exec.Command("rundll32", "url.dll,FileProtocolHandler", address)
		default:
			c = exec.Command("xdg-open", address)
		}
		c.Start()
	}
	signals := make(chan os.Signal, 1)
	signal.Notify(signals, os.Interrupt, syscall.SIGTERM, syscall.Signal(1))
	go func() {
		<-signals
		a.mu.Lock()
		a.closing = true
		a.mu.Unlock()
		for {
			a.mu.Lock()
			busy := a.busy
			a.mu.Unlock()
			if !busy {
				a.shutdown()
				return
			}
			time.Sleep(time.Second)
		}
	}()
	go func() {
		for {
			time.Sleep(time.Minute)
			a.mu.Lock()
			idle := !a.busy && time.Since(a.last) > 30*time.Minute
			if idle {
				a.closing = true
			}
			a.mu.Unlock()
			if idle {
				a.shutdown()
				return
			}
		}
	}()
	<-a.quit
	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()
	server.Shutdown(ctx)
	a.cleanup()
	fmt.Println("Closed. Results and decisions are retained; temporary MAFFT package released.")
}
