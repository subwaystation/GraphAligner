#!/usr/bin/python

import fileinput

nodes = []
edges = []

nodeid = {}
parent = {}
rank = {}
check_edges = {}

for line in fileinput.input():
	if line[0] == 'S':
		parts = line.split('\t')
		this_id = int(parts[1])
		nodes.append((line, this_id))
		nodeid[this_id] = this_id
		parent[this_id] = this_id
		rank[this_id] = 0
		check_edges[this_id] = []
	if line[0] == 'L':
		parts = line.split('\t')
		edges.append((line, int(parts[1]), int(parts[3])))
		check_edges[int(parts[1])].append(int(parts[3]))
		check_edges[int(parts[3])].append(int(parts[1]))

assert len(nodeid) == len(parent) == len(rank) == len(check_edges)

def find(index):
	global parent
	if parent[index] != index: parent[index] = find(parent[index])
	return parent[index]

def union(first, second):
	global parent
	global rank
	first_root = find(first)
	second_root = find(second)
	if first_root == second_root: return
	if rank[first_root] < rank[second_root]: first_root, second_root = second_root, first_root
	parent[second_root] = first_root
	if rank[second_root] == rank[first_root]: rank[first_root] += 1

while True:
	found_one = False
	for i in nodeid:
		for j in check_edges[i]:
			if find(i) != find(j):
				union(i, j)
				found_one = True
	if not found_one: break

nodes_by_component = {}
edges_by_component = {}

for node in nodes:
	root = find(node[1])
	if root not in nodes_by_component: 
		nodes_by_component[root] = []
		edges_by_component[root] = []
	nodes_by_component[root].append(node[0])

for edge in edges:
	root = find(edge[1])
	assert find(edge[1]) == find(edge[2])
	assert root in edges_by_component
	edges_by_component[root].append(edge[0])

max_by_nodes = None
max_by_edges = None

for component in nodes_by_component:
	if max_by_nodes == None or len(nodes_by_component[component]) > len(nodes_by_component[max_by_nodes]):
		max_by_nodes = component

for component in edges_by_component:
	if max_by_edges == None or len(edges_by_component[component]) > len(edges_by_component[max_by_edges]):
		max_by_edges = component

assert max_by_nodes != None
assert max_by_edges != None

for node in nodes_by_component[max_by_nodes]:
	print(node.strip())
for edge in edges_by_component[max_by_nodes]:
	print(edge.strip())
