#!/usr/bin/python

import json
import sys
from Gfa import *

ingraph = sys.argv[1]
insnarlstraversals = sys.argv[2]
outgraph = sys.argv[3]
outchains = sys.argv[4]

graph = Graph()
graph.load(ingraph)

with open(insnarlstraversals) as f:
	travs = json.load(f)

trav_edges = set()
parent = {}

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
	parent[p2] = p1

for n in graph.nodes:
	parent[n] = n

for trav in travs:
	prev = None
	nodes = []
	for n in trav["visit"]:
		if "node_id" in n:
			current = (n["node_id"], not ("backward" in n))
			nodes.append(n["node_id"])
			assert current[0] in graph.nodes
			if prev:
				trav_edges.add((prev, current))
				trav_edges.add((reverse(current), reverse(prev)))
			prev = current
	for n in nodes:
		assert nodes[0] in graph.nodes
		assert n in graph.nodes
		merge(nodes[0], n)

chains = Graph()

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

for c in chains.nodes:
	assert c in chain_start
	assert c in chain_end

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

graph.write(outgraph)
chains.write(outchains)
