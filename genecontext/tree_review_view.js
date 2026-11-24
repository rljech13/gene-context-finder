/* View and review only. All placement measurements and signatures come from Python. */
(function () {
  'use strict';
  const $ = id => document.getElementById(id);
  let manifest, token, reviews, panelId, groupId, previewResult, epoch = 0, dirty = false;
  let pending = false, exactMask = null;
  const ns = 'http://www.w3.org/2000/svg';
  const clone = x => JSON.parse(JSON.stringify(x));
  const panel = () => manifest.panels[panelId];
  const group = () => panel().groups[groupId];
  const request = () => reviews.panels[panelId].groups[groupId];
  const defaultMethod = () => ({model: 'MFP', bootstrap: 1000, seed: 42, threads: 2});
  const values = select => Array.from(select.selectedOptions, option => option.value);
  const number = value => value == null ? 'NA' : Number(value).toFixed(3);

  async function api(path, body) {
    const result = await fetch(path, body === undefined ? {cache: 'no-store'} : {
      method: 'POST', headers: {'Content-Type': 'application/json', 'X-GCF-Token': token}, body: JSON.stringify(body)
    });
    const data = await result.json();
    if (!result.ok) throw new Error(data.error || 'Local review request failed');
    return data;
  }
  function notice(text, kind = '') { $('notice').textContent = text; $('notice').className = 'status ' + kind; }
  function changed() { dirty = true; $('saved').textContent = 'Unsaved changes'; }
  function option(select, value, label) {
    const node = document.createElement('option'); node.value = value; node.textContent = label; select.appendChild(node);
  }
  function treeText(source) {
    if (source.newick) return source.newick;
    const saved = group().trees[source.id]; return typeof saved === 'string' ? saved : saved.newick;
  }
  function label(id) {
    const leaf = (group().leaves || []).find(l => l.tree_id === id);
    return leaf ? (leaf.labels || [id]).join(' / ') : id;
  }
  function outgroupChoices() {
    let ids;
    try { ids = tips(parseNewick(treeText(request().source))); } catch (_) { ids = []; }
    for (const name of ['buildOutgroup', 'reviewOutgroup']) {
      const select = $(name); select.replaceChildren();
      for (const id of ids) {
        option(select, id, label(id));
        const leaf = (group().leaves || []).find(l => l.tree_id === id);
        select.lastChild.selected = (request().outgroup_ids || []).includes(id) ||
          (leaf && (leaf.sequence_ids || []).some(alias => (request().outgroup_ids || []).includes(alias)));
      }
    }
  }
  function loadGroup() {
    epoch++; previewResult = null;
    const g = group(), r = request();
    $('members').replaceChildren();
    for (const id of Object.keys(g.rows).sort()) {
      const line = document.createElement('label'), input = document.createElement('input');
      input.type = 'checkbox'; input.value = id; input.checked = (r.selected_ids || []).includes(id);
      input.addEventListener('change', () => {
        request().selected_ids = Array.from($('members').querySelectorAll('input:checked'), n => n.value);
        invalidate('Sequence membership changed. Rebuild the tree before confirming it.');
      });
      line.append(input, document.createTextNode(' ' + label(id))); $('members').appendChild(line);
    }
    $('source').replaceChildren();
    for (const id of Object.keys(g.trees)) option($('source'), id, id);
    if (!(r.source.id in g.trees)) option($('source'), r.source.id, 'Rebuilt tree · ' + r.source.id.slice(-10));
    $('source').value = r.source.id;
    exactMask = Array.isArray(r.mask) ? clone(r.mask) : null;
    $('mask').value = exactMask ? 'exact' : (r.mask || {}).mode || 'all';
    $('mask').querySelector('option[value="exact"]').disabled = !exactMask;
    const m = r.tree_method || (r.source.provenance || {}).tree_method || g.tree_method || defaultMethod();
    for (const key of ['model', 'bootstrap', 'seed', 'threads']) $(key).value = m[key];
    $('buildBasis').value = $('reviewBasis').value = r.root_basis || '';
    $('comment').value = r.comment || '';
    outgroupChoices(); assess();
  }
  function loadPanel(wanted) {
    panelId = $('panel').value; $('group').replaceChildren();
    for (const [id, item] of Object.entries(panel().groups)) {
      if (reviews.panels[panelId].groups[id]) option($('group'), id, item.label || id);
    }
    if (wanted && reviews.panels[panelId].groups[wanted]) $('group').value = wanted;
    groupId = $('group').value;
    if (!groupId) { notice('This panel has no saved tree to review.', 'warning'); return; }
    loadGroup();
  }
  function invalidate(message) {
    const r = request(); r.approved = false; r.confirmed_signature = null;
    $('approved').checked = false; $('approved').disabled = true;
    $('gate').textContent = 'T is not included. The changed configuration needs a fresh exact review.';
    previewResult = null; epoch++; changed();
    notice(message || 'Configuration changed; previous confirmation cleared.'); assess();
  }
  async function assess() {
    const generation = ++epoch; $('approved').disabled = true;
    try {
      const item = {panel: panelId, seal: panel().seal, groups: {[groupId]: clone(request())}};
      const result = await api('/api/assess', item);
      if (generation !== epoch) return;
      previewResult = result.groups[groupId];
      reviews.panels[panelId].groups[groupId] = clone(previewResult.request);
      $('approved').checked = previewResult.approved;
      $('approved').disabled = previewResult.requires_rebuild;
      $('orientation').textContent = previewResult.outgroup_ids.length ?
        'Rooted display · outgroup: ' + previewResult.outgroup_ids.map(label).join(' / ') :
        'Explicitly unrooted · display orientation is arbitrary';
      $('signature').textContent = previewResult.signature;
      $('provenance').textContent = JSON.stringify(previewResult.provenance, null, 2);
      $('limitation').textContent = previewResult.placement.limitation;
      $('gate').textContent = previewResult.requires_rebuild ? 'T not included: membership, mask, or tree method changed; rebuild first.' :
        previewResult.approved ? 'Exact tree/outgroup reviewed. Assessed T may enter the score; unavailable model placements remain NA.' :
          'Proposed T only. No tree points enter the score until this exact tree/outgroup is reviewed.';
      renderModels(result.models); drawTree(previewResult.display_newick);
      notice(previewResult.requires_rebuild ? 'The displayed tree does not match the current construction inputs.' :
        previewResult.approved ? 'Tree and outgroup confirmed' + (dirty ? ' · save changes to persist.' : '.') :
          'Automatic placement ready · awaiting exact tree/outgroup review.', previewResult.requires_rebuild ? 'warning' : '');
    } catch (error) {
      if (generation !== epoch) return;
      request().approved = false; request().confirmed_signature = null;
      $('approved').checked = false; $('approved').disabled = true;
      $('gate').textContent = 'T is not included. Resolve the displayed validation error.';
      notice(error.message, 'error');
    }
  }
  function renderModels(results) {
    $('models').replaceChildren();
    for (const model of group().models || []) {
      const row = results[model.model_id]; if (!row) continue;
      const p = row.placement || {}, counts = p.counts || {};
      const tr = document.createElement('tr');
      const cells = [model.label || model.model_id, number(row.proposed_points), number(row.points),
        (counts.assessed || 0) + '/' + (counts.reference_resolved || 0) +
          (p.coverage == null ? '' : ' (' + (100 * p.coverage).toFixed(1) + '%)') +
          (p.sampling && p.sampling !== 'exhaustive' ? ' · sampled' : ''),
        (counts.concordant || 0) + ' / ' + (counts.discordant || 0),
        [p.status, row.reason, p.ambiguous_peer_species?.length ? 'Ambiguous peers: ' + p.ambiguous_peer_species.join(', ') : ''].filter(Boolean).join(' · ')];
      for (const value of cells) { const td = document.createElement('td'); td.textContent = value; tr.appendChild(td); }
      $('models').appendChild(tr);
    }
  }

  // Presentation parser only: no quartet, split selection or score logic in JS.
  function parseNewick(text) {
    if (text.length > 4 * 1024 * 1024) throw new Error('Tree exceeds display limit');
    const tokens = text.replace(/\[[^\]]*\]/g, '').match(/'(?:[^']|'')*'|[(),:;]|[^\s(),:;]+/g) || [];
    let at = 0; const used = new Set();
    function node(depth) {
      if (depth > 500) throw new Error('Tree exceeds display nesting limit');
      const n = {children: [], name: '', length: null};
      if (tokens[at] === '(') {
        at++; n.children.push(node(depth + 1));
        while (tokens[at] === ',') { at++; n.children.push(node(depth + 1)); }
        if (tokens[at++] !== ')' || n.children.length < 2) throw new Error('Invalid tree branching');
      }
      if (tokens[at] && ![':', ',', ')', ';'].includes(tokens[at])) {
        n.name = tokens[at++];
        if (n.name.startsWith("'")) n.name = n.name.slice(1, -1).replace(/''/g, "'");
      }
      if (!n.children.length) {
        if (!n.name || used.has(n.name)) throw new Error('Invalid duplicate tree leaf'); used.add(n.name);
      }
      if (tokens[at] === ':') {
        at++; n.length = Number(tokens[at++]);
        if (!Number.isFinite(n.length) || n.length < 0) throw new Error('Invalid branch length');
      }
      return n;
    }
    const result = node(0);
    if (tokens[at++] !== ';' || at !== tokens.length) throw new Error('Expected one Newick tree');
    return result;
  }
  function tips(n) { return n.children.length ? n.children.flatMap(tips) : [n.name]; }
  function drawTree(text) {
    const root = parseNewick(text), leaves = []; let maximum = 0;
    function place(n, distance) {
      n.x = distance; maximum = Math.max(maximum, distance);
      if (!n.children.length) { n.y = 24 + leaves.length * 32; leaves.push(n); }
      else { for (const c of n.children) place(c, distance + ($('lengths').checked && c.length != null ? c.length : 1)); n.y = n.children.reduce((sum, c) => sum + c.y, 0) / n.children.length; }
    }
    place(root, 0);
    const width = 1150, height = Math.max(200, leaves.length * 32 + 65), scale = 520 / Math.max(maximum, .00001);
    const svg = document.createElementNS(ns, 'svg'); svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
    svg.setAttribute('role', 'img'); svg.setAttribute('aria-label', 'Reviewed gene tree; labels and branch support');
    svg.style.width = $('zoom').value + '%';
    function element(tag, attrs, textValue) {
      const el = document.createElementNS(ns, tag); for (const [key, value] of Object.entries(attrs)) el.setAttribute(key, value);
      if (textValue !== undefined) el.textContent = textValue; svg.appendChild(el); return el;
    }
    function draw(n) {
      const x = 18 + n.x * scale;
      if (n.children.length) {
        const ys = n.children.map(c => c.y);
        element('path', {d: `M${x},${Math.min(...ys)} V${Math.max(...ys)}`, stroke: '#5c8081', fill: 'none'});
        for (const c of n.children) {
          element('path', {d: `M${x},${c.y} H${18 + c.x * scale}`, stroke: '#5c8081', fill: 'none'}); draw(c);
        }
        if (n.name) element('text', {x: x + 3, y: n.y - 5, 'font-size': 10, fill: '#496e70'}, n.name);
      } else {
        const t = element('text', {x: x + 6, y: n.y + 4, 'font-size': 12, fill: '#183b3c'}, label(n.name));
        const title = document.createElementNS(ns, 'title'); title.textContent = n.name + ' · ' + label(n.name); t.appendChild(title);
      }
    }
    draw(root); element('text', {x: 18, y: height - 14, 'font-size': 11, fill: '#556969'},
      $('lengths').checked ? 'Recorded branch lengths where available; missing lengths use display spacing.' : 'Topology display; branch lengths not to scale.');
    $('tree').replaceChildren(svg);
  }
  async function save(exportFile = false) {
    if (pending) return; pending = true; $('activity').textContent = 'Saving exact review…';
    try {
      const data = await api('/api/save', clone(reviews.panels[panelId])); reviews = data.reviews;
      dirty = false; $('saved').textContent = 'Saved in package .review folder';
      $('activity').textContent = 'Saved'; loadGroup();
      if (exportFile) {
        const url = URL.createObjectURL(new Blob([JSON.stringify(reviews, null, 2)], {type: 'application/json'}));
        const link = document.createElement('a'); link.href = url; link.download = 'gcf-tree-review-v2.json'; link.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
      }
    } catch (e) { notice(e.message, 'error'); $('activity').textContent = 'Save failed'; }
    finally { pending = false; }
  }
  function wire() {
    $('panel').addEventListener('change', () => loadPanel());
    $('group').addEventListener('change', () => { groupId = $('group').value; loadGroup(); });
    for (const name of ['buildOutgroup', 'reviewOutgroup']) $(name).addEventListener('change', () => {
      request().outgroup_ids = values($(name)); outgroupChoices(); invalidate('Outgroup changed; inspect and confirm the exact new orientation.');
    });
    for (const name of ['buildBasis', 'reviewBasis']) $(name).addEventListener('input', () => {
      request().root_basis = $(name).value; $(name === 'buildBasis' ? 'reviewBasis' : 'buildBasis').value = $(name).value;
      invalidate('Outgroup basis changed; previous confirmation cleared.');
    });
    $('unrooted').addEventListener('click', () => {
      request().outgroup_ids = []; request().root_basis = ''; $('buildBasis').value = $('reviewBasis').value = ''; outgroupChoices(); invalidate('Explicitly unrooted review selected.');
    });
    $('comment').addEventListener('input', () => { request().comment = $('comment').value; changed(); });
    $('mask').addEventListener('change', () => { request().mask = $('mask').value === 'exact' ? clone(exactMask) : {mode: $('mask').value}; invalidate('Mask changed; rebuild to obtain a matching tree.'); });
    $('source').addEventListener('change', () => {
      const id = $('source').value;
      if (id in group().trees) {
        const value = group().trees[id], newick = typeof value === 'string' ? value : value.newick;
        request().source = {id}; request().selected_ids = tips(parseNewick(newick)).sort(); request().mask = group().mask || {mode: 'all'};
        request().tree_method = group().tree_method || defaultMethod(); request().outgroup_ids = []; request().root_basis = '';
      }
      request().approved = false; request().confirmed_signature = null; changed(); loadGroup();
    });
    for (const name of ['model', 'bootstrap', 'seed', 'threads']) $(name).addEventListener('change', () => {
      request().tree_method = {model: $('model').value, bootstrap: Number($('bootstrap').value), seed: Number($('seed').value), threads: Number($('threads').value)};
      invalidate('Tree-construction method changed; rebuild first.');
    });
    $('approved').addEventListener('change', () => {
      request().approved = $('approved').checked; request().confirmed_signature = request().approved && previewResult ? previewResult.signature : null;
      changed(); assess();
    });
    $('rebuild').addEventListener('click', async () => {
      if (pending) return; pending = true; $('rebuild').disabled = true;
      request().approved = false; request().confirmed_signature = null; $('approved').checked = false; $('approved').disabled = true; epoch++; changed();
      notice('Building an IQ-TREE result for the selected alignment…');
      try {
        const r = request(), data = await api('/api/rebuild', {panel: panelId, seal: panel().seal, group: groupId,
          inputSha256: r.inputSha256, selected_ids: r.selected_ids, mask: r.mask, tree_method: r.tree_method || defaultMethod(), allow_download: $('download').checked});
        Object.assign(r, data); loadGroup();
      } catch (e) { notice(e.message, 'error'); }
      finally { pending = false; $('rebuild').disabled = false; }
    });
    $('save').addEventListener('click', () => save()); $('export').addEventListener('click', () => save(true));
    $('import').addEventListener('change', async () => {
      try {
        const file = $('import').files[0]; if (!file || file.size > 12 * 1024 * 1024) throw new Error('Select a review JSON under 12 MiB');
        const input = JSON.parse(await file.text());
        if (input.schema !== 'gcf_tree_review_v2') throw new Error('Only v2 automatic reviews can be imported; legacy manual T is not migrated as approved');
        const result = await api('/api/save', input); reviews = result.reviews; dirty = false; $('saved').textContent = 'Imported and revalidated'; loadPanel(groupId);
      } catch (e) { notice(e.message, 'error'); }
    });
    for (const name of ['zoom', 'lengths']) $(name).addEventListener('input', () => { if (previewResult) drawTree(previewResult.display_newick); });
    window.addEventListener('beforeunload', event => { if (dirty) { event.preventDefault(); event.returnValue = ''; } });
  }
  async function start() {
    if (location.protocol === 'file:') {
      notice('Start the local server first: python START_REVIEW.py. This standalone HTML cannot run or save tree review by itself.', 'warning'); return;
    }
    try {
      const init = await api('/api/manifest'); manifest = init.manifest; token = init.token;
      const state = await api('/api/state'); reviews = state.reviews;
      for (const [id, p] of Object.entries(manifest.panels)) option($('panel'), id, p.label || id);
      const params = new URLSearchParams(location.search); if (params.get('panel') in manifest.panels) $('panel').value = params.get('panel');
      wire(); loadPanel(params.get('group'));
      if (state.warnings.length) notice('Saved review needs revision: ' + state.warnings.join('; '), 'warning');
    } catch (e) { notice(e.message + ' · Run python START_REVIEW.py in the extracted package.', 'error'); }
  }
  if (typeof module !== 'undefined' && module.exports) module.exports = {parseNewick, tips};
  if (typeof document !== 'undefined') start();
}());
