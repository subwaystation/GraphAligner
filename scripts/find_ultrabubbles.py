#!/usr/bin/python

# Detecting Superbubbles in Assembly Graphs, Onodera et al 2013

import sys
from Gfa import *

ingraph = sys.argv[1]

graph = Graph()
graph.load(ingraph)
graph.remove_nonexistent_edges()

def getset(s):
	assert(len(s) == 1)
	for item in s:
		return item

# Fig. 5
def detect_bubble(s):
	if s not in graph.edges: return None
	if len(graph.edges[s]) < 2: return None
	S = [s]
	visited = set()
	seen = set()
	seen.add(s)
	while len(S) > 0:
		v = S.pop()
		assert v in seen
		seen.remove(v)
		assert v not in visited
		visited.add(v)
		if v not in graph.edges: return None
		if len(graph.edges[v]) == 0: return None
		for edge in graph.edges[v]:
			u = edge[0]
			if u[0] == v[0]: return None
			if reverse(u) in visited: return None
			if u == s: return None
			assert u not in visited
			seen.add(u)
			assert reverse(u) in graph.edges
			assert len(graph.edges[reverse(u)]) >= 1
			has_nonvisited_parent = False
			for parent_edge in graph.edges[reverse(u)]:
				parent = reverse(parent_edge[0])
				if parent not in visited: has_nonvisited_parent = True
			if not has_nonvisited_parent: S.append(u)
		if len(S) == 1 and len(seen) == 1 and S[0] == getset(seen):
			t = S.pop()
			for edge in graph.edges[t]:
				if edge[0] == t: return None
			return (s, t)
	return None


for n in graph.nodes:
	bubble = detect_bubble((n, True))
	if bubble:
		leftkey = '{"node_id":"' + str(bubble[0][0]) + '"'
		if not bubble[0][1]: leftkey += ',"backward":true'
		leftkey += '}'
		rightkey = '{"node_id":"' + str(bubble[1][0]) + '"'
		if not bubble[1][1]: rightkey += ',"backward":true'
		rightkey += '}'
		print('{"type":"ULTRABUBBLE","start":' + leftkey + ',"end":' + rightkey + '}')
	bubble = detect_bubble((n, False))
	if bubble:
		leftkey = '{"node_id":"' + str(bubble[0][0]) + '"'
		if not bubble[0][1]: leftkey += ',"backward":true'
		leftkey += '}'
		rightkey = '{"node_id":"' + str(bubble[1][0]) + '"'
		if not bubble[1][1]: rightkey += ',"backward":true'
		rightkey += '}'
		print('{"type":"ULTRABUBBLE","start":' + leftkey + ',"end":' + rightkey + '}')
