#!/usr/bin/python

import sys
import copy
from Gfa import *

infile = sys.argv[1]
outfile = sys.argv[2]

graph = Graph()
graph.load(infile)

added_nodes = []
removed_nodes = []
added_edges = []
next_id = 1

for n in graph.nodes:
	next_id = max(next_id, int(n))
next_id += 1

for n in graph.nodes:
	fwpos = (n, True)
	bwpos = (n, False)
	if len(graph.edges[fwpos]) != 2 or len(graph.edges[bwpos]) != 2: continue
	other_edge = None
	weird = False
	for fw in graph.edges[fwpos]:
		for bw in graph.edges[bwpos]:
			if fw[0] == reverse(bw[0]):
				if other_edge != None: weird = True
				other_edge = fw
	if weird: continue
	if other_edge == None: continue
	other_node = other_edge[0][0]
	fw_node = None
	bw_node = None
	cov_sum = 0
	cov_count = 0
	for fw in graph.edges[fwpos]:
		if fw[0][0] == other_node: continue
		fw_node = fw[0]
		cov_sum += fw[1][1]
		cov_count += 1
	for bw in graph.edges[bwpos]:
		if bw[0][0] == other_node: continue
		bw_node = bw[0]
		cov_sum += bw[1][1]
		cov_count += 1
	if fw_node == None: continue
	if bw_node == None: continue
	if fw_node == bw_node: continue
	if fw_node[0] == n: continue
	if bw_node[0] == n: continue
	if cov_count != 2: continue
	avg_cov = float(cov_sum) / float(cov_count)
	this_votes = int(round(graph.nodes[n].frequency / avg_cov))
	if other_node != n:
		other_votes = int(round(graph.nodes[other_node].frequency / avg_cov))
		if other_votes != this_votes - 1: 
			continue
		removed_nodes.append(n)
		removed_nodes.append(other_node)
		for i in range(0, this_votes):
			if i == 0:
				added_edges += [(next_id, False, target[0][0], target[0][1], target[1]) for target in graph.edges[(n, False)]]
			if i > 0:
				added_edges.append((next_id-1, True, next_id, other_edge[0][1], other_edge[1]))
				added_edges.append((next_id, other_edge[0][1], next_id+1, True, other_edge[1]))
				added_nodes.append(copy.copy(graph.nodes[other_node]))
				added_nodes[-1].nodeid = next_id
				added_nodes[-1].readcount = added_nodes[-1].readcount / other_votes
				next_id += 1
			added_nodes.append(copy.copy(graph.nodes[n]))
			added_nodes[-1].nodeid = next_id
			added_nodes[-1].readcount = added_nodes[-1].readcount / this_votes
			if i == this_votes - 1:
				added_edges += [(next_id, True, target[0][0], target[0][1], target[1]) for target in graph.edges[(n, True)]]
			next_id += 1
	else:
		print(n)
		removed_nodes.append(n)
		for i in range(0, this_votes):
			if i == 0:
				added_edges += [(next_id, False, target[0][0], target[0][1], target[1]) for target in graph.edges[(n, False)]]
			if i > 0:
				added_edges.append((next_id-1, True, next_id, True, other_edge[1]))
			added_nodes.append(copy.copy(graph.nodes[n]))
			added_nodes[-1].nodeid = next_id
			added_nodes[-1].readcount = added_nodes[-1].readcount / this_votes
			if i == this_votes - 1:
				added_edges += [(next_id, True, target[0][0], target[0][1], target[1]) for target in graph.edges[(n, True)]]
			next_id += 1

for node in added_nodes:
	assert node.nodeid not in graph.nodes
	graph.nodes[node.nodeid] = node

for node in removed_nodes:
	assert node in graph.nodes
	del graph.nodes[node]

for edge in added_edges:
	frompos = (edge[0], edge[1])
	topos = (edge[2], edge[3])
	info = edge[4]
	if frompos not in graph.edges: graph.edges[frompos] = set()
	graph.edges[frompos].add((topos, info))

graph.remove_nonexistent_edges()
graph.write(outfile)
