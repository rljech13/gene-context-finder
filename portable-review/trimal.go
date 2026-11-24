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
	"os/exec"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"time"
)

//go:embed trimal_manifest.json
var trimalManifest []byte

func canonicalResidue(c byte) bool {
	return strings.ContainsRune("ACDEFGHIKLMNPQRSTVWY", rune(c))
}

func trimALExecutable(work string) (string, error) {
	if path, err := exec.LookPath("trimal"); err == nil {
		return path, nil
	}
	key := runtime.GOOS + "-" + runtime.GOARCH
	if key == "windows-arm64" {
		key = "windows-amd64"
	}
	var packages map[string]Package
	if err := json.Unmarshal(trimalManifest, &packages); err != nil {
		return "", err
	}
	p, ok := packages[key]
	if !ok {
		return "", errors.New("Install trimAl for this architecture")
	}
	response, err := (&http.Client{Timeout: 45 * time.Second}).Get(p.URL)
	if err != nil {
		return "", err
	}
	defer response.Body.Close()
	if response.StatusCode != http.StatusOK {
		return "", fmt.Errorf("trimAl download returned HTTP %d", response.StatusCode)
	}
	raw, err := io.ReadAll(io.LimitReader(response.Body, 4<<20))
	if err != nil {
		return "", err
	}
	if digest(raw) != p.SHA256 {
		return "", errors.New("trimAl archive checksum differs")
	}
	archive := filepath.Join(work, p.Archive)
	if err := os.WriteFile(archive, raw, 0600); err != nil {
		return "", err
	}
	if err := unpack(archive, work); err != nil {
		return "", err
	}
	path := filepath.Join(work, filepath.FromSlash(p.Entry))
	if runtime.GOOS != "windows" {
		if err := os.Chmod(path, 0700); err != nil {
			return "", err
		}
	}
	return path, nil
}

func trimALColumnMap(log string, width int) ([]int, error) {
	for _, line := range strings.Split(log, "\n") {
		mapText, ok := strings.CutPrefix(line, "#ColumnsMap")
		if !ok {
			continue
		}
		mapText = strings.TrimSpace(mapText)
		columns := []int{}
		if mapText == "" {
			return columns, nil
		}
		last := -1
		for _, value := range strings.Split(mapText, ",") {
			index, err := strconv.Atoi(strings.TrimSpace(value))
			if err != nil || index <= last || index >= width {
				return nil, errors.New("Invalid trimAl column map")
			}
			columns = append(columns, index+1)
			last = index
		}
		return columns, nil
	}
	return nil, errors.New("Missing trimAl column map")
}

func trimALColumns(rows map[string]string, width int) ([]int, map[string]any, error) {
	work, err := markedTemp("gcf-review-align-")
	if err != nil {
		return nil, nil, err
	}
	defer os.RemoveAll(work)
	exe, err := trimALExecutable(work)
	if err != nil {
		return nil, nil, err
	}
	version, err := command(exe, "--version").CombinedOutput()
	if err != nil || !strings.Contains(strings.ToLower(string(version)), "trimal") {
		return nil, nil, errors.New("Compatible trimAl executable is required")
	}
	binary, err := os.ReadFile(exe)
	if err != nil {
		return nil, nil, err
	}
	ids := make([]string, 0, len(rows))
	for id := range rows {
		ids = append(ids, id)
	}
	sort.Strings(ids)
	input, out := filepath.Join(work, "source.faa"), filepath.Join(work, "trimmed.faa")
	if err := os.WriteFile(input, fasta(rows, ids), 0600); err != nil {
		return nil, nil, err
	}
	args := []string{"-in", input, "-out", out, "-automated1", "-colnumbering"}
	cmd := command(exe, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout, cmd.Stderr = &stdout, &stderr
	if err := cmd.Start(); err != nil {
		return nil, nil, err
	}
	owner, _ := json.Marshal(map[string]any{"schema": "gcf-review-temp-v1", "pid": cmd.Process.Pid})
	_ = os.WriteFile(filepath.Join(work, ".gcf-owner.json"), owner, 0600)
	if err := cmd.Wait(); err != nil {
		return nil, nil, fmt.Errorf("trimAl failed: %v: %s", err, stderr.String())
	}
	columns, err := trimALColumnMap(stdout.String(), width)
	if err != nil {
		return nil, nil, err
	}
	if len(columns) > 0 {
		raw, err := os.ReadFile(out)
		if err != nil {
			return nil, nil, err
		}
		trimmed, err := parseFASTA(raw)
		if err != nil {
			return nil, nil, err
		}
		if len(trimmed) != len(rows) {
			return nil, nil, errors.New("trimAl changed the input sequence set")
		}
		for id, source := range rows {
			projected := make([]byte, len(columns))
			for i, column := range columns {
				projected[i] = source[column-1]
			}
			if trimmed[id] != string(projected) {
				return nil, nil, fmt.Errorf("trimAl output disagrees with its column map for %s", id)
			}
		}
	}
	return columns, map[string]any{"method": "trimal_automated1", "tool_version": strings.TrimSpace(string(version)), "tool_sha256": digest(binary), "options": []string{"-automated1", "-colnumbering"}}, nil
}
