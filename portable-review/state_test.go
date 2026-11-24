package main

import (
	"bytes"
	"encoding/json"
	"net/http/httptest"
	"strings"
	"testing"
)

func TestPortableStateIncludesLargeCheckpointsAndIsolatesReviews(t *testing.T) {
	a := &App{root: t.TempDir(), raw: []byte("exact synthetic dataset"), data: Dataset{ReviewID: "review", SourceReviewID: "source"}}
	for _, key := range []string{"gcf-realign-selection-review-result-g", "gcf-review-workspace-3-review", "gcf-stage2-source"} {
		if !a.stateKeyAllowed(key) {
			t.Fatal("Rejected current review key:", key)
		}
	}
	for _, key := range []string{"gcf-realign-selection-reviewOTHER-result-g", "gcf-stage2-sourceOTHER", "gcf-stage2-unrelated"} {
		if a.stateKeyAllowed(key) {
			t.Fatal("Accepted another review key:", key)
		}
	}
	key := "gcf-realign-selection-review-saved-alignment-capacity-control"
	v := browserState{Schema: "gcf_portable_session_v2", DatasetSHA: digest(a.raw), Values: map[string]string{key: strings.Repeat("x", 18<<20)}}
	raw, _ := json.Marshal(v)
	w := httptest.NewRecorder()
	a.handleState(w, httptest.NewRequest("POST", "http://localhost/api/state", bytes.NewReader(raw)))
	if w.Code != 200 {
		t.Fatal(w.Code, w.Body.String())
	}
	restored := a.readState()
	if restored.Values[key] != v.Values[key] {
		t.Fatal("Large checkpoint did not survive restart")
	}
	a.raw = []byte("changed synthetic dataset")
	if len(a.readState().Values) != 0 {
		t.Fatal("State transferred to a changed dataset")
	}
}
