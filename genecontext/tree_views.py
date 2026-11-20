"""Two displays of one ML topology; rooting never reruns or renames taxa."""
from html import escape
import math
from copy import deepcopy
from .stage4_phylogeny import parse_newick, tip_names, tree_svg, phylogram_scale


def graph(tree):
    adjacent, names = {}, {}
    def visit(node, parent=None):
        i=len(adjacent);adjacent[i]={};names[i]=node['name']
        if parent is not None:
            edge=(max(0.,node['length']), node['name'] if node['children'] else '')
            adjacent[i][parent]=edge;adjacent[parent][i]=edge
        for child in node['children']:visit(child,i)
        return i
    root=visit(tree)
    if len(adjacent[root])==2:
        a,b=adjacent[root];la,sa=adjacent[a].pop(root);lb,sb=adjacent[b].pop(root)
        adjacent[a][b]=adjacent[b][a]=(la+lb,sa or sb)
        del adjacent[root];del names[root]
    return adjacent,names


def component_tips(adj,names,start,blocked):
    out=set();stack=[(start,blocked)]
    while stack:
        node,parent=stack.pop()
        if len(adj[node])==1:out.add(names[node])
        stack.extend((n,node) for n in adj[node] if n!=parent)
    return out


def rooted_tree(newick,outgroups):
    """Place a root halfway along the exact outgroup/ingroup edge.

    An arbitrary MRCA is never used to force non-monophyletic outgroups.
    """
    adj,names=graph(parse_newick(newick));wanted=set(outgroups)
    all_tips={names[n] for n in adj if len(adj[n])==1}
    if not wanted or not wanted<all_tips:raise ValueError('Outgroup must be a nonempty strict subset of tree leaves')
    split=None
    for a in adj:
        for b in adj[a]:
            if component_tips(adj,names,a,b)==wanted:split=(a,b);break
        if split:break
    if split is None:raise ValueError('Selected outgroup is not a separate clade in this topology; no forced rooting was made')
    a,b=split;length,support=adj[a][b]
    def descend(node,parent,length,support):
        children=[descend(n,node,*adj[node][n]) for n in adj[node] if n!=parent]
        return dict(children=children,name=support if children else names[node],length=length)
    root=dict(children=[descend(a,b,length/2,''),descend(b,a,length/2,'')],name='',length=0.)
    return root,dict(method='midpoint_of_explicit_outgroup_edge',outgroup_tree_ids=sorted(wanted),
                     root_edge_length=length,root_edge_support=support)


def newick_text(tree):
    def emit(n):
        text='('+','.join(emit(c) for c in n['children'])+')'+n['name'] if n['children'] else n['name']
        return text+(':'+format(n['length'],'.12g') if n is not tree else '')
    return emit(tree)+';'


def unrooted_svg(newick,leaves,scale_lengths=False):
    adj,names=graph(parse_newick(newick));labels={l['tree_id']:'; '.join(l['labels']) for l in leaves}
    root=next((n for n in adj if len(adj[n])>=3),next(iter(adj)))
    positions={root:(0.,0.)}
    def count(n,parent):return 1 if len(adj[n])==1 else sum(count(c,n) for c in adj[n] if c!=parent)
    total=sum(count(n,root) for n in adj[root])
    def place(n,parent,start,end):
        angle=(start+end)/2;length=adj[n][parent][0] if scale_lengths else 1.
        x,y=positions[parent];positions[n]=(x+length*math.cos(angle),y+length*math.sin(angle))
        children=[c for c in adj[n] if c!=parent];sizes=[count(c,n) for c in children];lo=start
        for c,size in zip(children,sizes):
            hi=lo+(end-start)*size/sum(sizes);place(c,n,lo,hi);lo=hi
    lo=-math.pi
    for n in adj[root]:
        hi=lo+2*math.pi*count(n,root)/total;place(n,root,lo,hi);lo=hi
    xs=[p[0] for p in positions.values()];ys=[p[1] for p in positions.values()]
    sides={n:positions[n][0]>=positions[root][0] for n in adj if len(adj[n])==1}
    width=1400;height=max(600,max(sum(v==side for v in sides.values()) for side in [True,False])*28+100);scale=min(600/max(max(xs)-min(xs),1e-9),(height-100)/max(max(ys)-min(ys),1e-9))
    def xy(n):
        x,y=positions[n];return (700+(x-(max(xs)+min(xs))/2)*scale,45+(y-min(ys))*scale)
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="Unrooted maximum-likelihood gene tree"><rect width="100%" height="100%" fill="white"/><g font-family="system-ui" font-size="11">']
    for a in adj:
        for b,(length,support) in adj[a].items():
            if a>b:continue
            x,y=xy(a);u,v=xy(b);out.append(f'<path d="M{x},{y} L{u},{v}" stroke="#527b89" fill="none"/>')
            if support and math.hypot(x-u,y-v)>35:out.append(f'<text x="{(x+u)/2+3}" y="{(y+v)/2-3}" fill="#56707b" font-size="9">{escape(support)}</text>')
    for side in [False,True]:
        tips=sorted((n for n in sides if sides[n]==side),key=lambda n:xy(n)[1]);ys_label=[]
        for n in tips:ys_label.append(max(xy(n)[1],ys_label[-1]+28 if ys_label else 30))
        if ys_label and ys_label[-1]>height-45:
            ys_label[-1]=height-45
            for i in range(len(ys_label)-2,-1,-1):ys_label[i]=min(ys_label[i],ys_label[i+1]-28)
        for n,label_y in zip(tips,ys_label):
            x,y=xy(n);label=labels.get(names[n],names[n]);label_x=1040 if side else 360;anchor='start' if side else 'end'
            parts=label.split(' · ',1);first=parts[0];second=parts[1] if len(parts)>1 else ''
            second=second if len(second)<48 else second[:47]+'…'
            out.append(f'<path d="M{x},{y} L{label_x+(-8 if side else 8)},{label_y}" stroke="#c1cdd2" fill="none"/>')
            out.append(f'<text data-tree-label="{names[n]}" x="{label_x}" y="{label_y}" text-anchor="{anchor}"><title>{escape(label)}</title>{escape(first)}<tspan x="{label_x}" dy="12" font-size="10">{escape(second)}</tspan></text>')
    if scale_lengths:
        out.append(phylogram_scale(scale,18,height-28))
    else:
        out.append(f'<text x="18" y="{height-14}">Unrooted · topology layout; branch lengths not to scale</text>')
    out.append('</g></svg>')
    return ''.join(out)


def views(newick,leaves,outgroup_ids=(),root_basis=''):
    result=dict(unrooted_newick=newick,unrooted_svg=unrooted_svg(newick,leaves),unrooted_phylogram_svg=unrooted_svg(newick,leaves,True),rooting='unrooted',
                root_basis='',outgroup_sequence_ids=list(outgroup_ids))
    if not outgroup_ids:return result
    if not root_basis.strip():raise ValueError('Explain the biological basis for the selected outgroup')
    ids=set(outgroup_ids);allids={s for l in leaves if l['included'] for s in l['sequence_ids']}
    if ids-allids:raise ValueError('Outgroup sequence missing after masking')
    if any(ids.intersection(l['sequence_ids']) and set(l['sequence_ids'])-ids for l in leaves if l['included']):
        raise ValueError('Identical proteins cannot be split between outgroup and ingroup; select all exact aliases or a distinct homolog')
    outer=[l['tree_id'] for l in leaves if l['included'] and ids.intersection(l['sequence_ids'])]
    root,metadata=rooted_tree(newick,outer)
    text=newick_text(root)
    topology=deepcopy(root)
    def unit(n):
        for child in n['children']:child['length']=1.;unit(child)
    unit(topology)
    svg=tree_svg(newick_text(topology),leaves,scale_bar=False)
    result.update(rooted_newick=text,rooted_svg=svg,rooted_phylogram_svg=tree_svg(text,leaves),rooting='explicit_user_outgroup',
                  root_basis=root_basis,root_metadata=metadata)
    return result
