#!/usr/bin/python

import json
import sys
from Gfa import *

ingraph = sys.argv[1]
insnarls = sys.argv[2]
outgraph = sys.argv[3]
outchains = sys.argv[4]

print("read graph")

graph = Graph()
graph.load(ingraph)
graph.remove_nonexistent_edges()

snarls = []

print("read ultrabubbles")

with open(insnarls) as f:
	for l in f:
		parts = l.strip().split('\t')
		if len(parts) < 4: continue
		snarls.append(((parts[0]), parts[1] == "+", (parts[2]), parts[3] == "+"))

trav_edges = set()
parent = {}
rank = {}

def find(k1):
	global parent
	assert k1 in parent
	chain = []
	while parent[k1] != k1:
		chain.append(k1)
		k1 = parent[k1]
	for k in chain:
		parent[k] = k1
	return k1

def merge(k1, k2):
	global parent
	p1 = find(k1)
	p2 = find(k2)
	if rank[p1] < rank[p2]:
		(p1, p2) = (p2, p1)
	parent[p2] = p1
	if rank[p2] == rank[p1]: rank[p1] += 1

for n in graph.nodes:
	parent[n] = n
	rank[n] = 0

print("merge bubble nodes")

for snarl in snarls:
	prev = None
	start_node = (snarl[0], snarl[1])
	end_node = (snarl[2], snarl[3])
	nodes = set()
	visit_stack = [start_node]
	trav_here = set()
	while len(visit_stack) > 0:
		pos = visit_stack.pop()
		if pos[0] in nodes: continue
		nodes.add(pos[0])
		if pos == end_node: continue
		if pos not in graph.edges: print(pos)
		assert pos in graph.edges
		for edge in graph.edges[pos]:
			target = edge[0]
			if (pos, target) in trav_here:
				assert (reverse(target), reverse(pos)) in trav_here
				continue
			assert (reverse(target), reverse(pos)) not in trav_here
			trav_here.add((pos, target))
			trav_here.add((reverse(target), reverse(pos)))
			visit_stack.append(target)
	trav_edges = trav_edges.union(trav_here)
	assert len(nodes) > 0
	for n in nodes:
		assert n in graph.nodes
		merge(start_node[0], n)

chains = Graph()

print("merge chains")

for n in graph.nodes:
	node = graph.nodes[n]
	p = find(n)
	graph.nodes[n].chain = str(p)
	if p not in chains.nodes:
		newnode = Node()
		newnode.nodeid = p
		newnode.nodeseq = "*"
		newnode.chain = p
		chains.nodes[p] = newnode
	chains.nodes[p].length += node.length
	chains.nodes[p].readcount += node.readcount

chain_start = {}
chain_end = {}

def add_end(n, pos):
	global chain_start
	global chain_end
	if n not in chain_start:
		assert n not in chain_end 
		chain_start[n] = reverse(pos)
		return
	assert n not in chain_end
	chain_end[n] = pos

print("add ends")

for n in graph.nodes:
	fwpos = (n, True)
	bwpos = (n, False)
	if len(graph.edges[fwpos]) == 0:
		add_end(find(n), fwpos)
	else:
		edge_exists = None
		for e in graph.edges[fwpos]:
			if edge_exists == None:
				edge_exists = (fwpos, e[0]) in trav_edges
			assert edge_exists == ((fwpos, e[0]) in trav_edges)
		if not edge_exists: add_end(find(n), fwpos)
	if len(graph.edges[bwpos]) == 0:
		add_end(find(n), bwpos)
	else:
		edge_exists = None
		for e in graph.edges[bwpos]:
			if edge_exists == None:
				edge_exists = (bwpos, e[0]) in trav_edges
			assert edge_exists == ((bwpos, e[0]) in trav_edges)
		if not edge_exists: add_end(find(n), bwpos)

print("add circular ends")

for c in chains.nodes:
	if (c in chain_start) and (c in chain_end): continue
	assert c not in chain_start
	assert c not in chain_end
	# circular chain
	chains.edges[(c, True)] = set()
	chains.edges[(c, True)].add(((c, True), (0, None)))
	chain_start[c] = (-1, True)
	chain_end[c] = (-1, True)

for c in chains.nodes:
	assert c in chain_start
	assert c in chain_end

print("add chain edges")

for e in graph.edges:
	if e == chain_end[find(e[0])]:
		fromchain = (find(e[0]), True)
	elif reverse(e) == chain_start[find(e[0])]:
		fromchain = (find(e[0]), False)
	else: continue
	for target in graph.edges[e]:
		topos = target[0]
		if topos == chain_start[find(topos[0])]:
			tochain = (find(topos[0]), True)
		elif reverse(topos) == chain_end[find(topos[0])]:
			tochain = (find(topos[0]), False)
		else: assert False 
		if fromchain not in chains.edges: chains.edges[fromchain] = set()
		chains.edges[fromchain].add((tochain, target[1]))

print("write")

graph.write(outgraph)
chains.write(outchains)
