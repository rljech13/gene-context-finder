package main

import (
	"math"
	"regexp"
	"strconv"
	"strings"
	"testing"
)

func TestPhylogramScaleUsesThePlottedBranchScale(t *testing.T) {
	root, err := parseTree("(A:0.1,B:0.2);")
	if err != nil {
		t.Fatal(err)
	}
	svg := treeSVG(root, nil)
	// With the longer branch occupying 510 SVG units, the scale is 2550 units/site.
	bar := regexp.MustCompile(`data-substitutions-per-site="([^"]+)"><path d="M[^ ]+ h([^ ]+)`).FindStringSubmatch(svg)
	if len(bar) != 3 {
		t.Fatal("Scale bar is missing or unlabelled")
	}
	sites, _ := strconv.ParseFloat(bar[1], 64)
	pixels, _ := strconv.ParseFloat(bar[2], 64)
	if math.Abs(pixels-sites*2550) > .001 || pixels > 100 {
		t.Fatalf("Scale bar = %g units for %g substitutions/site", pixels, sites)
	}
	radial, err := unrootedSVG("(A:0.1,B:0.2,C:0.3,D:0.4);", nil)
	if err != nil || !strings.Contains(radial, "tree-distance-scale") {
		t.Fatalf("Unrooted scale bar missing: %v", err)
	}
}

func TestTrimOffKeepsEveryColumn(t *testing.T) {
	source := map[string]string{"a": "A-C", "b": "C--"}
	_, _, mask, err := maskedTree(source, []string{"a", "b"}, nil, "none")
	if err != nil {
		t.Fatal(err)
	}
	columns := mask["retained_columns_1based"].([]int)
	if len(columns) != 3 || columns[0] != 1 || columns[1] != 2 || columns[2] != 3 {
		t.Fatalf("Trim-off column map: %v", columns)
	}
}
