package main

import (
	"bytes"
	_ "embed"
	"encoding/json"
	"errors"
	"fmt"
	"html"
	"io"
	"math"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"time"
)

//go:embed iqtree_manifest.json
var iqtreeManifest []byte
var buildHash = "development"
var resultPattern = regexp.MustCompile(`^[a-f0-9]{20,64}$`)

func validResultID(id string) bool { return resultPattern.MatchString(id) }
func (a *App) verifiedResult(id string) (Result, error) {
	var r Result
	if !validResultID(id) {
		return r, errors.New("Invalid saved result ID")
	}
	root := filepath.Join(a.root, "review_results", id)
	b, e := os.ReadFile(filepath.Join(root, "checksums.json"))
	if e != nil {
		return r, e
	}
	var seals map[string]string
	if e = json.Unmarshal(b, &seals); e != nil {
		return r, e
	}
	if _, ok := seals["result.json"]; !ok {
		return r, errors.New("Result seal is missing")
	}
	for name, want := range seals {
		p, e := safePath(root, name)
		if e != nil {
			return r, e
		}
		raw, e := os.ReadFile(p)
		if e != nil || digest(raw) != want {
			return r, fmt.Errorf("Saved result checksum failed: %s", name)
		}
	}
	raw, e := os.ReadFile(filepath.Join(root, "result.json"))
	if e != nil {
		return r, e
	}
	if e = json.Unmarshal(raw, &r); e != nil {
		return r, e
	}
	if len(r.IDs) < 2 || r.ID != id || r.ReviewID != a.data.ReviewID {
		return r, errors.New("Saved result belongs to another review")
	}
	if e = a.validateResultSequences(r); e != nil {
		return r, e
	}
	return r, nil
}
func (a *App) validateResultSequences(r Result) error {
	n := -1
	seen := map[string]bool{}
	for _, id := range r.IDs {
		s, ok := a.byID[id]
		if !ok || s.Group != r.Group || seen[id] {
			return errors.New("Unknown, duplicate or foreign model in result")
		}
		seen[id] = true
		v, ok := r.Protein[id]
		if !ok || strings.ReplaceAll(v, "-", "") != s.Protein {
			return fmt.Errorf("Saved protein changed: %s", id)
		}
		if n >= 0 && len(v) != n {
			return errors.New("Unequal saved alignment rows")
		}
		n = len(v)
	}
	if len(seen) != len(r.Protein) {
		return errors.New("Unlisted saved alignment row")
	}
	selected := map[string]Sequence{}
	for id := range seen {
		selected[id] = a.byID[id]
	}
	cds, e := codingAlignment(r.Protein, selected)
	if e != nil {
		return e
	}
	b, _ := json.Marshal(cds)
	c, _ := json.Marshal(r.CDS)
	if !bytes.Equal(b, c) {
		return errors.New("Saved CDS differs from exact native CDS")
	}
	return nil
}
func probeIQTREE(path string) (string, error) {
	out, e := command(path, "--version").CombinedOutput()
	if e != nil || !strings.Contains(strings.ToUpper(string(out)), "IQ-TREE") {
		return "", errors.New("IQ-TREE executable cannot run")
	}
	lines := strings.Split(strings.TrimSpace(string(out)), "\n")
	return lines[0], nil
}
func (a *App) ensureIQTREE() error {
	if a.iqtree != "" {
		return nil
	}
	for _, name := range []string{"iqtree3", "iqtree2", "iqtree"} {
		if p, e := exec.LookPath(name); e == nil {
			if v, e := probeIQTREE(p); e == nil {
				a.iqtree = p
				a.iqtreeVersion = v
				return nil
			}
		}
	}
	packages := map[string]Package{}
	json.Unmarshal(iqtreeManifest, &packages)
	key := runtime.GOOS + "-" + runtime.GOARCH
	if key == "windows-arm64" {
		key = "windows-amd64"
	}
	p, ok := packages[key]
	if !ok {
		return errors.New("Install a compatible IQ-TREE executable for this architecture")
	}
	root := filepath.Join(os.TempDir(), "gcf-review-iqtree-"+p.SHA256[:16])
	unlock, e := cacheLock(root)
	if e != nil {
		return e
	}
	defer unlock()
	os.MkdirAll(root, 0700)
	os.WriteFile(filepath.Join(root, ".gcf-cache-owner"), []byte("gcf-review-cache-v1"), 0600)
	entry := filepath.Join(root, filepath.FromSlash(p.Entry))
	if _, e = os.Stat(entry); e != nil {
		work, e := markedTemp("gcf-review-download-")
		if e != nil {
			return e
		}
		defer os.RemoveAll(work)
		response, e := (&http.Client{Timeout: 5 * time.Minute}).Get(p.URL)
		if e != nil {
			return e
		}
		defer response.Body.Close()
		if response.StatusCode != 200 {
			return fmt.Errorf("IQ-TREE download HTTP %d", response.StatusCode)
		}
		raw, e := io.ReadAll(io.LimitReader(response.Body, 128<<20))
		if e != nil {
			return e
		}
		if digest(raw) != p.SHA256 {
			return errors.New("IQ-TREE download checksum mismatch")
		}
		file := filepath.Join(work, p.Archive)
		if e = os.WriteFile(file, raw, 0600); e != nil {
			return e
		}
		if e = unpack(file, root); e != nil {
			return e
		}
	}
	if runtime.GOOS != "windows" {
		if e = os.Chmod(entry, 0700); e != nil {
			return e
		}
	}
	v, e := probeIQTREE(entry)
	if e != nil {
		return e
	}
	os.MkdirAll(filepath.Join(root, "leases"), 0700)
	if e = os.WriteFile(filepath.Join(root, "leases", strconv.Itoa(os.Getpid())), []byte(a.token), 0600); e != nil {
		return e
	}
	a.iqtreeTemporary = root
	a.iqtree = entry
	a.iqtreeVersion = v
	return nil
}
func (a *App) cleanupIQTREE() {
	root := a.iqtreeTemporary
	if root == "" {
		return
	}
	unlock, e := cacheLock(root)
	if e != nil {
		return
	}
	defer unlock()
	os.Remove(filepath.Join(root, "leases", strconv.Itoa(os.Getpid())))
	entries, _ := os.ReadDir(filepath.Join(root, "leases"))
	for _, f := range entries {
		pid, _ := strconv.Atoi(f.Name())
		if alive(pid) {
			return
		}
	}
	os.RemoveAll(root)
}

type treeNode struct {
	Name     string
	Length   float64
	Children []*treeNode
}

var newickTokens = regexp.MustCompile(`\(|\)|,|:|;|[^\s(),:;]+`)

func parseTree(text string) (*treeNode, error) {
	tokens := newickTokens.FindAllString(text, -1)
	pos := 0
	var parse func() (*treeNode, error)
	parse = func() (*treeNode, error) {
		n := &treeNode{}
		if pos >= len(tokens) {
			return nil, errors.New("Truncated Newick")
		}
		if tokens[pos] == "(" {
			pos++
			for {
				child, e := parse()
				if e != nil {
					return nil, e
				}
				n.Children = append(n.Children, child)
				if pos >= len(tokens) {
					return nil, errors.New("Truncated Newick children")
				}
				if tokens[pos] == ")" {
					pos++
					break
				}
				if tokens[pos] != "," {
					return nil, errors.New("Malformed Newick")
				}
				pos++
			}
		}
		if pos < len(tokens) && !strings.Contains("(:,);", tokens[pos]) {
			n.Name = tokens[pos]
			pos++
		}
		if n.Name == "" && len(n.Children) == 0 {
			return nil, errors.New("Missing Newick leaf name")
		}
		if pos < len(tokens) && tokens[pos] == ":" {
			pos++
			if pos >= len(tokens) {
				return nil, errors.New("Missing branch length")
			}
			v, e := strconv.ParseFloat(tokens[pos], 64)
			if e != nil || v < 0 || math.IsNaN(v) || math.IsInf(v, 0) {
				return nil, errors.New("Invalid branch length")
			}
			n.Length = v
			pos++
		}
		return n, nil
	}
	root, e := parse()
	if e != nil {
		return nil, e
	}
	if pos != len(tokens)-1 || tokens[pos] != ";" {
		return nil, errors.New("Malformed Newick terminator")
	}
	return root, nil
}
func emitTree(root *treeNode) string {
	var emit func(*treeNode, bool) string
	emit = func(n *treeNode, first bool) string {
		v := n.Name
		if len(n.Children) > 0 {
			parts := []string{}
			for _, c := range n.Children {
				parts = append(parts, emit(c, false))
			}
			v = "(" + strings.Join(parts, ",") + ")" + n.Name
		}
		if !first {
			v += ":" + strconv.FormatFloat(n.Length, 'g', 12, 64)
		}
		return v
	}
	return emit(root, true) + ";"
}

type edge struct {
	Length  float64
	Support string
}
type treeGraph struct {
	Adj   map[int]map[int]edge
	Names map[int]string
}

func graphTree(root *treeNode) treeGraph {
	g := treeGraph{map[int]map[int]edge{}, map[int]string{}}
	var walk func(*treeNode, int) int
	walk = func(n *treeNode, parent int) int {
		id := len(g.Adj)
		g.Adj[id] = map[int]edge{}
		g.Names[id] = n.Name
		if parent >= 0 {
			support := ""
			if len(n.Children) > 0 {
				support = n.Name
			}
			e := edge{n.Length, support}
			g.Adj[id][parent] = e
			g.Adj[parent][id] = e
		}
		for _, c := range n.Children {
			walk(c, id)
		}
		return id
	}
	r := walk(root, -1)
	if len(g.Adj[r]) == 2 {
		ns := neighbors(g.Adj[r])
		a, b := ns[0], ns[1]
		ea, eb := g.Adj[r][a], g.Adj[r][b]
		support := ea.Support
		if support == "" {
			support = eb.Support
		}
		e := edge{ea.Length + eb.Length, support}
		delete(g.Adj[a], r)
		delete(g.Adj[b], r)
		g.Adj[a][b] = e
		g.Adj[b][a] = e
		delete(g.Adj, r)
		delete(g.Names, r)
	}
	return g
}
func neighbors(adj map[int]edge) []int {
	v := []int{}
	for id := range adj {
		v = append(v, id)
	}
	sort.Ints(v)
	return v
}
func component(g treeGraph, start, blocked int) map[string]bool {
	out := map[string]bool{}
	var walk func(int, int)
	walk = func(n, p int) {
		if len(g.Adj[n]) == 1 {
			out[g.Names[n]] = true
		}
		for _, c := range neighbors(g.Adj[n]) {
			if c != p {
				walk(c, n)
			}
		}
	}
	walk(start, blocked)
	return out
}
func sameSet(a, b map[string]bool) bool {
	if len(a) != len(b) {
		return false
	}
	for x := range a {
		if !b[x] {
			return false
		}
	}
	return true
}
func rootTree(newick string, outgroups []string) (*treeNode, error) {
	original, e := parseTree(newick)
	if e != nil {
		return nil, e
	}
	g := graphTree(original)
	wanted := map[string]bool{}
	for _, id := range outgroups {
		wanted[id] = true
	}
	a, b := -1, -1
	nodes := []int{}
	for n := range g.Adj {
		nodes = append(nodes, n)
	}
	sort.Ints(nodes)
	for _, n := range nodes {
		for _, c := range neighbors(g.Adj[n]) {
			if sameSet(component(g, n, c), wanted) {
				a, b = n, c
				break
			}
		}
		if a >= 0 {
			break
		}
	}
	if a < 0 {
		return nil, errors.New("Selected outgroup is not a separate clade; no forced root was made")
	}
	var descend func(int, int, float64, string) *treeNode
	descend = func(n, p int, length float64, support string) *treeNode {
		v := &treeNode{Length: length}
		for _, c := range neighbors(g.Adj[n]) {
			if c != p {
				x := g.Adj[n][c]
				v.Children = append(v.Children, descend(c, n, x.Length, x.Support))
			}
		}
		if len(v.Children) == 0 {
			v.Name = g.Names[n]
		} else {
			v.Name = support
		}
		return v
	}
	length := g.Adj[a][b].Length
	return &treeNode{Children: []*treeNode{descend(a, b, length/2, ""), descend(b, a, length/2, "")}}, nil
}
func phylogramScale(scale, x, y float64) string {
	target := 100 / scale
	order := math.Pow(10, math.Floor(math.Log10(target)))
	value := order
	for _, multiplier := range []float64{2, 5} {
		if multiplier*order <= target {
			value = multiplier * order
		}
	}
	end := x + value*scale
	return fmt.Sprintf(`<g class="tree-distance-scale" data-substitutions-per-site="%.8g"><path d="M%.6g,%.6g h%.6g M%.6g,%.6g v8 M%.6g,%.6g v8" stroke="#27464a" fill="none"/><text x="%.6g" y="%.6g">%.8g substitutions / site</text></g>`,
		value, x, y, value*scale, x, y-4, end, y-4, x, y+17, value)
}

func treeSVG(root *treeNode, labels map[string]string) string {
	type point struct{ x, y float64 }
	points := map[*treeNode]point{}
	tips := []*treeNode{}
	maximum := 0.
	var place func(*treeNode, float64)
	place = func(n *treeNode, x float64) {
		if len(n.Children) == 0 {
			points[n] = point{x, 35 + float64(len(tips))*28}
			tips = append(tips, n)
			maximum = math.Max(maximum, x)
			return
		}
		y := 0.
		for _, c := range n.Children {
			place(c, x+c.Length)
			y += points[c].y
		}
		points[n] = point{x, y / float64(len(n.Children))}
	}
	place(root, 0)
	scale := 510 / math.Max(maximum, .001)
	height := 70 + len(tips)*28
	var b strings.Builder
	fmt.Fprintf(&b, `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1300 %d"><rect width="100%%" height="100%%" fill="white"/><g font-family="system-ui" font-size="12" fill="#27464a">`, height)
	var draw func(*treeNode)
	draw = func(n *treeNode) {
		p := points[n]
		x := 20 + p.x*scale
		if len(n.Children) == 0 {
			label := labels[n.Name]
			if label == "" {
				label = n.Name
			}
			fmt.Fprintf(&b, `<text x="%.5f" y="%.5f"><title>%s</title>%s</text>`, x+6, p.y+4, html.EscapeString(label), html.EscapeString(label))
			return
		}
		ys := []float64{}
		for _, c := range n.Children {
			ys = append(ys, points[c].y)
		}
		sort.Float64s(ys)
		fmt.Fprintf(&b, `<path d="M%.5f,%.5f V%.5f" stroke="#527b89" fill="none"/>`, x, ys[0], ys[len(ys)-1])
		for _, c := range n.Children {
			v := points[c]
			fmt.Fprintf(&b, `<path d="M%.5f,%.5f H%.5f" stroke="#527b89" fill="none"/>`, x, v.y, 20+v.x*scale)
			draw(c)
		}
		if n.Name != "" {
			fmt.Fprintf(&b, `<text x="%.5f" y="%.5f" font-size="9">%s</text>`, x+3, p.y-4, html.EscapeString(n.Name))
		}
	}
	draw(root)
	b.WriteString(phylogramScale(scale, 20, float64(height-28)))
	b.WriteString("</g></svg>")
	return b.String()
}
func unrootedSVG(newick string, labels map[string]string) (string, error) {
	r, e := parseTree(newick)
	if e != nil {
		return "", e
	}
	g := graphTree(r)
	root := -1
	nodes := []int{}
	for n := range g.Adj {
		nodes = append(nodes, n)
	}
	sort.Ints(nodes)
	for _, n := range nodes {
		if len(g.Adj[n]) >= 3 {
			root = n
			break
		}
	}
	if root < 0 {
		return "", errors.New("Too few leaves for unrooted display")
	}
	type pt struct{ x, y float64 }
	pos := map[int]pt{root: {0, 0}}
	var count func(int, int) int
	count = func(n, p int) int {
		if len(g.Adj[n]) == 1 {
			return 1
		}
		sum := 0
		for _, c := range neighbors(g.Adj[n]) {
			if c != p {
				sum += count(c, n)
			}
		}
		return sum
	}
	total := count(root, -1)
	var place func(int, int, float64, float64)
	place = func(n, p int, lo, hi float64) {
		angle := (lo + hi) / 2
		length := g.Adj[n][p].Length
		parent := pos[p]
		pos[n] = pt{parent.x + length*math.Cos(angle), parent.y + length*math.Sin(angle)}
		sum := count(n, p)
		start := lo
		for _, c := range neighbors(g.Adj[n]) {
			if c == p {
				continue
			}
			end := start + (hi-lo)*float64(count(c, n))/float64(sum)
			place(c, n, start, end)
			start = end
		}
	}
	lo := -math.Pi
	for _, c := range neighbors(g.Adj[root]) {
		hi := lo + 2*math.Pi*float64(count(c, root))/float64(total)
		place(c, root, lo, hi)
		lo = hi
	}
	minx, maxx, miny, maxy := 0., 0., 0., 0.
	for _, v := range pos {
		minx = math.Min(minx, v.x)
		maxx = math.Max(maxx, v.x)
		miny = math.Min(miny, v.y)
		maxy = math.Max(maxy, v.y)
	}
	height := math.Max(600, float64(total)*20+100)
	scale := math.Min(600/math.Max(maxx-minx, 1e-9), (height-100)/math.Max(maxy-miny, 1e-9))
	xy := func(n int) pt { v := pos[n]; return pt{700 + (v.x-(maxx+minx)/2)*scale, 45 + (v.y-miny)*scale} }
	var b strings.Builder
	fmt.Fprintf(&b, `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1400 %.0f"><rect width="100%%" height="100%%" fill="white"/><g font-family="system-ui" font-size="11">`, height)
	for _, n := range nodes {
		for _, c := range neighbors(g.Adj[n]) {
			if n > c {
				continue
			}
			a, z := xy(n), xy(c)
			fmt.Fprintf(&b, `<path d="M%.5f,%.5f L%.5f,%.5f" stroke="#527b89" fill="none"/>`, a.x, a.y, z.x, z.y)
		}
	}
	for _, side := range []bool{false, true} {
		tips := []int{}
		for _, n := range nodes {
			if len(g.Adj[n]) == 1 && (pos[n].x >= pos[root].x) == side {
				tips = append(tips, n)
			}
		}
		sort.Slice(tips, func(i, j int) bool { return xy(tips[i]).y < xy(tips[j]).y })
		ys := []float64{}
		for _, n := range tips {
			y := xy(n).y
			if len(ys) > 0 {
				y = math.Max(y, ys[len(ys)-1]+28)
			}
			ys = append(ys, y)
		}
		if len(ys) > 0 && ys[len(ys)-1] > height-45 {
			ys[len(ys)-1] = height - 45
			for i := len(ys) - 2; i >= 0; i-- {
				ys[i] = math.Min(ys[i], ys[i+1]-28)
			}
		}
		for i, n := range tips {
			x, anchor := 360., "end"
			if side {
				x = 1040
				anchor = "start"
			}
			p := xy(n)
			label := labels[g.Names[n]]
			if label == "" {
				label = g.Names[n]
			}
			parts := strings.SplitN(label, " · ", 2)
			extra := ""
			if len(parts) > 1 {
				extra = parts[1]
			}
			fmt.Fprintf(&b, `<path d="M%.5f,%.5f L%.5f,%.5f" stroke="#c1cdd2" fill="none"/><text data-tree-label="%s" x="%.5f" y="%.5f" text-anchor="%s"><title>%s</title>%s<tspan x="%.5f" dy="12" font-size="10">%s</tspan></text>`, p.x, p.y, x, ys[i], html.EscapeString(g.Names[n]), x, ys[i], anchor, html.EscapeString(label), html.EscapeString(parts[0]), x, html.EscapeString(extra))
		}
	}
	b.WriteString(phylogramScale(scale, 18, height-28))
	b.WriteString("</g></svg>")
	return b.String(), nil
}

type treeLeaf struct {
	TreeID      string   `json:"tree_id"`
	ProteinSHA  string   `json:"protein_sha256"`
	SequenceIDs []string `json:"sequence_ids"`
	Labels      []string `json:"labels"`
	Included    bool     `json:"included"`
}

func maskedTree(alignment map[string]string, ids []string, byID map[string]Sequence, trimMode string) (map[string]string, []treeLeaf, map[string]any, error) {
	aliases := map[string]int{}
	leaves := []treeLeaf{}
	raw := map[string]string{}
	width := -1
	for _, id := range ids {
		seq, ok := alignment[id]
		if !ok {
			return nil, nil, nil, errors.New("Tree row missing from alignment")
		}
		if width >= 0 && len(seq) != width {
			return nil, nil, nil, errors.New("Unequal alignment rows")
		}
		width = len(seq)
		h := digest([]byte(strings.ReplaceAll(seq, "-", "")))
		index, ok := aliases[h]
		if !ok {
			index = len(leaves)
			aliases[h] = index
			tid := fmt.Sprintf("T%05d", index+1)
			leaves = append(leaves, treeLeaf{TreeID: tid, ProteinSHA: h, Included: true})
			raw[tid] = seq
		}
		leaves[index].SequenceIDs = append(leaves[index].SequenceIDs, id)
		leaves[index].Labels = append(leaves[index].Labels, byID[id].Label)
	}
	columns := []int{}
	trim := map[string]any{"method": trimMode}
	if trimMode == "trimal_automated1" && len(raw) >= 4 && width >= 30 {
		for id, row := range raw {
			coded := []byte(row)
			for i, c := range coded {
				if !canonicalResidue(c) {
					coded[i] = '-'
				}
			}
			raw[id] = string(coded)
		}
		var err error
		columns, trim, err = trimALColumns(raw, width)
		if err != nil {
			return nil, nil, nil, err
		}
	} else if trimMode == "trimal_automated1" {
		trim = map[string]any{"method": "trimal_automated1", "status": "not_run", "reason": "Fewer than 4 unique taxa or 30 input columns"}
	} else if trimMode == "none" {
		for i := 1; i <= width; i++ {
			columns = append(columns, i)
		}
	} else {
		for i := 0; i < width; i++ {
			n := 0
			for _, seq := range raw {
				if canonicalResidue(seq[i]) {
					n++
				}
			}
			if float64(n)/float64(len(raw)) >= .5 {
				columns = append(columns, i+1)
			}
		}
	}
	masked := map[string]string{}
	excluded := []string{}
	for i := range leaves {
		leaf := &leaves[i]
		var b strings.Builder
		has := false
		for _, column := range columns {
			c := raw[leaf.TreeID][column-1]
			if !canonicalResidue(c) {
				c = '-'
			} else {
				has = true
			}
			b.WriteByte(c)
		}
		leaf.Included = has
		if has {
			masked[leaf.TreeID] = b.String()
		} else {
			excluded = append(excluded, leaf.TreeID)
		}
	}
	variable := 0
	for i := range columns {
		states := map[byte]bool{}
		for _, seq := range masked {
			if seq[i] != '-' {
				states[seq[i]] = true
			}
		}
		if len(states) > 1 {
			variable++
		}
	}
	mask := map[string]any{"original_columns": width, "retained_columns": len(columns), "retained_columns_1based": columns, "variable_columns": variable, "trimming": trim, "excluded_all_gap_taxa": excluded}
	if trimMode == "occupancy_50" {
		mask["occupancy"] = .5
	}
	return masked, leaves, mask, nil
}
func (a *App) tree(q Request) (Result, error) {
	var r Result
	if len(q.IDs) < 2 || len(q.IDs) > 2000 {
		return r, errors.New("Select 2–2000 distinct sequences")
	}
	if len(q.OutgroupIDs) > 0 && strings.TrimSpace(q.RootBasis) == "" {
		return r, errors.New("Explain the basis for the selected outgroup")
	}
	source, e := a.verifiedResult(q.AlignmentResultID)
	if e != nil {
		return r, e
	}
	if source.Phase != "alignment" || source.Group != q.Group {
		return r, errors.New("Tree input must be the saved alignment for this group")
	}
	ids := append([]string{}, q.IDs...)
	sort.Strings(ids)
	seen := map[string]bool{}
	protein, cds := map[string]string{}, map[string]string{}
	for _, id := range ids {
		if seen[id] {
			return r, errors.New("Duplicate tree row")
		}
		seen[id] = true
		v, ok := source.Protein[id]
		if !ok {
			return r, errors.New("Tree row is outside the saved alignment")
		}
		protein[id] = v
		if v, ok := source.CDS[id]; ok {
			cds[id] = v
		}
	}
	outSet := map[string]bool{}
	for _, id := range q.OutgroupIDs {
		if outSet[id] || !seen[id] {
			return r, errors.New("Unknown or duplicate outgroup")
		}
		outSet[id] = true
	}
	if len(outSet) == len(ids) {
		return r, errors.New("Outgroup cannot contain the entire alignment")
	}
	trim := q.TrimMode
	if trim == "" {
		trim = "none"
	}
	if trim != "none" && trim != "occupancy_50" && trim != "trimal_automated1" {
		return r, errors.New("Unknown alignment trimming method")
	}
	masked, leaves, mask, e := maskedTree(protein, ids, a.byID, trim)
	if e != nil {
		return r, e
	}
	outer := []string{}
	labels := map[string]string{}
	for _, leaf := range leaves {
		labels[leaf.TreeID] = strings.Join(leaf.Labels, "; ")
		n := 0
		for _, id := range leaf.SequenceIDs {
			if outSet[id] {
				n++
			}
		}
		if n > 0 {
			if !leaf.Included {
				return r, errors.New("No aligned sites remain for the selected outgroup")
			}
			if n != len(leaf.SequenceIDs) {
				return r, errors.New("Identical proteins cannot be split between outgroup and ingroup")
			}
			outer = append(outer, leaf.TreeID)
		}
	}
	coreRaw, _ := json.Marshal([]any{"portable_tree_v2", buildHash, digest(a.raw), source.ID, ids, mask})
	core := digest(coreRaw)
	signature, _ := json.Marshal([]any{core, q.OutgroupIDs, q.RootBasis})
	id := digest(signature)[:20]
	if _, e = os.Stat(filepath.Join(a.root, "review_results", id, "result.json")); e == nil {
		return a.verifiedResult(id)
	}
	r = Result{ID: id, ReviewID: a.data.ReviewID, InputSHA: digest(signature), IDs: ids, Group: q.Group, Protein: protein, CDS: cds, MAFFT: source.MAFFT, Created: time.Now().UTC().Format(time.RFC3339), Phase: "tree", AlignmentID: source.ID, Provenance: map[string]any{"method": "portable_review23_v2", "helper_source_sha256": buildHash, "alignment_result_id": source.ID, "tree_core_signature": core}}
	phy := map[string]any{"status": "not_computed", "method": "Protein ML; ModelFinder BIC; UFBoot 1000; BNNI; seed 42", "leaves": leaves, "mask": mask, "outgroup_sequence_ids": q.OutgroupIDs, "rooting": "unrooted", "root_basis": "", "bootstrap_replicates": 1000, "core_signature": core}
	files := map[string][]byte{"tree.aligned.faa": fasta(masked, func() []string {
		v := []string{}
		for k := range masked {
			v = append(v, k)
		}
		sort.Strings(v)
		return v
	}())}
	if len(masked) < 4 || mask["retained_columns"].(int) < 30 || mask["variable_columns"].(int) < 5 {
		if len(outer) > 0 {
			return r, errors.New("Insufficient taxa or aligned sites to use a selected outgroup. Clear it to save the diagnostic.")
		}
		phy["reason"] = "Insufficient taxa or variable aligned sites: requires 4 unique taxa, 30 retained columns and 5 variable columns"
	} else {
		// Reorient an already verified identical ML core without recomputing it.
		entries, _ := os.ReadDir(filepath.Join(a.root, "review_results"))
		reused := ""
		for _, entry := range entries {
			if !entry.IsDir() || !validResultID(entry.Name()) {
				continue
			}
			old, err := a.verifiedResult(entry.Name())
			if err == nil && old.Phase == "tree" && old.Phylogeny["core_signature"] == core && old.Phylogeny["status"] == "computed" {
				raw, _ := json.Marshal(old.Phylogeny)
				json.Unmarshal(raw, &phy)
				reused = entry.Name()
				break
			}
		}
		newick := ""
		if reused != "" {
			newick, _ = phy["unrooted_newick"].(string)
			phy["reused_result_id"] = reused
		} else {
			if e = a.ensureIQTREE(); e != nil {
				return r, e
			}
			work, e := markedTemp("gcf-review-align-")
			if e != nil {
				return r, e
			}
			defer os.RemoveAll(work)
			input := filepath.Join(work, "tree.aligned.faa")
			if e = os.WriteFile(input, files["tree.aligned.faa"], 0600); e != nil {
				return r, e
			}
			prefix := filepath.Join(work, "ml")
			argv := []string{"-s", input, "-st", "AA", "-keep-ident", "-m", "MFP", "-B", "1000", "-bnni", "-seed", "42", "-T", "2", "--prefix", prefix}
			c := command(a.iqtree, argv...)
			var stdout, stderr bytes.Buffer
			c.Stdout = &stdout
			c.Stderr = &stderr
			started := time.Now()
			if e = c.Start(); e != nil {
				return r, e
			}
			owner, _ := json.Marshal(map[string]any{"schema": "gcf-review-temp-v1", "pid": c.Process.Pid})
			os.WriteFile(filepath.Join(work, ".gcf-owner.json"), owner, 0600)
			lease := ""
			if a.iqtreeTemporary != "" {
				lease = filepath.Join(a.iqtreeTemporary, "leases", strconv.Itoa(c.Process.Pid))
				os.WriteFile(lease, []byte(a.token), 0600)
			}
			e = c.Wait()
			if lease != "" {
				os.Remove(lease)
			}
			if e != nil {
				return r, fmt.Errorf("IQ-TREE failed: %v\n%s\n%s", e, stdout.String(), stderr.String())
			}
			raw, e := os.ReadFile(prefix + ".treefile")
			if e != nil {
				return r, e
			}
			newick = strings.TrimSpace(string(raw))
			binary, _ := os.ReadFile(a.iqtree)
			phy["tool_version"] = a.iqtreeVersion
			phy["tool_sha256"] = digest(binary)
			phy["elapsed_seconds"] = time.Since(started).Seconds()
			phy["commands"] = []any{append([]string{"iqtree"}, []string{"-s", "tree.aligned.faa", "-st", "AA", "-keep-ident", "-m", "MFP", "-B", "1000", "-bnni", "-seed", "42", "-T", "2", "--prefix", "ml"}...)}
			files["iqtree.stdout.log"] = stdout.Bytes()
			files["iqtree.stderr.log"] = stderr.Bytes()
			for _, suffix := range []string{".iqtree", ".log", ".treefile", ".contree"} {
				if raw, e := os.ReadFile(prefix + suffix); e == nil {
					files["ml"+suffix] = raw
				}
			}
		}
		parsed, e := parseTree(newick)
		if e != nil {
			return r, e
		}
		actual := map[string]bool{}
		var visit func(*treeNode) error
		visit = func(n *treeNode) error {
			if len(n.Children) == 0 {
				if actual[n.Name] {
					return errors.New("IQ-TREE returned a duplicate leaf")
				}
				actual[n.Name] = true
			}
			for _, child := range n.Children {
				if err := visit(child); err != nil {
					return err
				}
			}
			return nil
		}
		if err := visit(parsed); err != nil {
			return r, err
		}
		expected := map[string]bool{}
		for tid := range masked {
			expected[tid] = true
		}
		if !sameSet(actual, expected) {
			return r, errors.New("IQ-TREE leaf set differs from exact masked input")
		}
		unrooted, e := unrootedSVG(newick, labels)
		if e != nil {
			return r, e
		}
		for _, key := range []string{"rooted_newick", "rooted_svg", "rooted_phylogram_svg", "root_metadata"} {
			delete(phy, key)
		}
		phy["status"] = "computed"
		phy["unrooted_newick"] = newick
		phy["unrooted_svg"] = unrooted
		phy["unrooted_phylogram_svg"] = unrooted
		phy["rooting"] = "unrooted"
		phy["root_basis"] = ""
		phy["outgroup_sequence_ids"] = q.OutgroupIDs
		if len(outer) > 0 {
			root, e := rootTree(newick, outer)
			if e != nil {
				return r, e
			}
			phy["rooted_newick"] = emitTree(root)
			view := treeSVG(root, labels)
			phy["rooted_svg"] = view
			phy["rooted_phylogram_svg"] = view
			phy["rooting"] = "explicit_user_outgroup"
			phy["root_basis"] = q.RootBasis
			phy["root_metadata"] = map[string]any{"method": "midpoint_of_explicit_outgroup_edge", "outgroup_tree_ids": outer}
		}
	}
	phy["tree_input_fasta"] = string(files["tree.aligned.faa"])
	r.Phylogeny = phy
	files["proteins.aligned.faa"] = fasta(r.Protein, r.IDs)
	files["cds.codon.aligned.fna"] = fasta(r.CDS, r.IDs)
	raw, _ := json.MarshalIndent(r, "", "  ")
	files["result.json"] = raw
	return r, a.saveResult(r.ID, files)
}

// Re-render a saved result from its Newick on read. The sealed result stays intact.
func scaledTreeViews(r *Result) error {
	if r.Phase != "tree" || r.Phylogeny["status"] != "computed" {
		return nil
	}
	labels := map[string]string{}
	rawLeaves, err := json.Marshal(r.Phylogeny["leaves"])
	if err != nil {
		return err
	}
	var leaves []treeLeaf
	if err := json.Unmarshal(rawLeaves, &leaves); err != nil {
		return err
	}
	for _, leaf := range leaves {
		labels[leaf.TreeID] = strings.Join(leaf.Labels, "; ")
	}
	for _, mode := range []string{"unrooted", "rooted"} {
		key := mode + "_phylogram_svg"
		if svg, ok := r.Phylogeny[key].(string); ok && strings.Contains(svg, "tree-distance-scale") {
			continue
		}
		newick, _ := r.Phylogeny[mode+"_newick"].(string)
		if newick == "" {
			continue
		}
		var svg string
		if mode == "rooted" {
			root, e := parseTree(newick)
			if e != nil {
				return e
			}
			svg = treeSVG(root, labels)
		} else {
			var e error
			svg, e = unrootedSVG(newick, labels)
			if e != nil {
				return e
			}
		}
		r.Phylogeny[key] = svg
		r.Phylogeny[mode+"_svg"] = svg
	}
	return nil
}

func (a *App) saveResult(id string, files map[string][]byte) error {
	parent := filepath.Join(a.root, "review_results")
	if e := os.MkdirAll(parent, 0700); e != nil {
		return e
	}
	tmp, e := os.MkdirTemp(parent, ".pending-")
	if e != nil {
		return e
	}
	defer os.RemoveAll(tmp)
	seal := map[string]string{}
	for name, raw := range files {
		if e = copyFile(tmp, name, 0600, bytes.NewReader(raw)); e != nil {
			return e
		}
		seal[name] = digest(raw)
	}
	raw, _ := json.MarshalIndent(seal, "", "  ")
	if e = os.WriteFile(filepath.Join(tmp, "checksums.json"), raw, 0600); e != nil {
		return e
	}
	return os.Rename(tmp, filepath.Join(parent, id))
}
