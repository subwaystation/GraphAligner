#!/usr/bin/python

import sys
import json
from Gfa import *

in_assembly = sys.argv[1]
in_bubbles = sys.argv[2]
out_collapsed = sys.argv[3]

graph = Graph()
graph.load(in_assembly)
graph.remove_nonexistent_edges()

def getfirst(s):
	assert(len(s) >= 1)
	for item in s:
		return item

def get_bubble_nodes(node, end_node):
	global graph
	check = [node]
	bubble_nodes = set()
	while len(check) > 0:
		pos = check.pop()
		if pos[0] in bubble_nodes: continue
		bubble_nodes.add(pos[0])
		if pos == end_node: continue
		for edge in graph.edges[pos]:
			check.append(edge[0])
	return bubble_nodes

def get_bubble_edges(node, end_node):
	global graph
	check = [node]
	bubble_edges = set()
	while len(check) > 0:
		pos = check.pop()
		for edge in graph.edges[pos]:
			if (pos, edge[0]) in bubble_edges: continue
			bubble_edges.add((pos, edge[0]))
			if edge[0] != end_node: check.append(edge[0])
	return bubble_edges

def collapse(start_node, end_node):
	global graph
	bubble_nodes = get_bubble_nodes(start_node, end_node)
	bubble_edges = get_bubble_edges(start_node, end_node)
	path = [start_node]
	while path[-1] != end_node:
		path.append(getfirst(graph.edges[path[-1]])[0])
	assert path[0] == start_node
	assert path[-1] == end_node
	for n in path:
		assert n[0] in bubble_nodes
		bubble_nodes.remove(n[0])
	for i in range(1, len(path)):
		assert (path[i-1], path[i]) in bubble_edges
		bubble_edges.remove((path[i-1], path[i]))
	for n in bubble_nodes:
		assert n in graph.nodes
		del graph.nodes[n]
	for e in bubble_edges:
		assert e[0] in graph.edges
		assert e[1] in graph.edges
		for e2 in graph.edges[e[0]]:
			if e2[0] == e[1]:
				graph.edges[e[0]].remove(e2)
				break
		else:
			assert False
		for e2 in graph.edges[reverse(e[1])]:
			if e2[0] == reverse(e[0]):
				graph.edges[reverse(e[1])].remove(e2)
				break
		else:
			assert False
	for n in path:
		assert n[0] in graph.nodes
	for i in range(1, len(path)):
		assert path[i-1] in graph.edges
		for e in graph.edges[path[i-1]]:
			if e[0] == path[i]:
				break
		else:
			assert False

snarls = []

with open(in_bubbles) as f:
	for l in f:
		parts = l.strip().split('\t')
		if len(parts) < 4: continue
		snarls.append(((parts[0]), parts[1] == "+", (parts[2]), parts[3] == "+"))

for snarl in snarls:
	start_node = (snarl[0], snarl[1])
	end_node = (snarl[2], snarl[3])
	if start_node[0] not in graph.nodes or end_node[0] not in graph.nodes: continue
	collapse(start_node, end_node)

graph.remove_nonexistent_edges()
graph.write(out_collapsed)
