#!/usr/bin/python

from Gfa import Graph
import sys

inputgraph = sys.argv[1]
minsize = int(sys.argv[2])
outputgraph = sys.argv[3]


graph = Graph()
graph.load(inputgraph)

checked = set()
removed = set()

def stack_check(n):
	global checked
	global graph
	stack = [n]
	total = 0
	while len(stack) > 0:
		n = stack[-1]
		stack.pop()
		if n in checked: continue
		checked.add(n)
		total += len(graph.nodes[n].nodeseq)
		for e in graph.edges[(n, True)]:
			stack.append(e[0][0])
		for e in graph.edges[(n, False)]:
			stack.append(e[0][0])
	return total

def stack_remove(n):
	global removed
	global graph
	stack = [n]
	while len(stack) > 0:
		n = stack[-1]
		stack.pop()
		if n in removed: continue
		removed.add(n)
		for e in graph.edges[(n, True)]:
			stack.append(e[0][0])
		for e in graph.edges[(n, False)]:
			stack.append(e[0][0])


for n in graph.nodes:
	if n in checked: continue
	size = stack_check(n)
	if size < minsize: stack_remove(n)

for n in removed:
	del graph.nodes[n]

graph.remove_nonexistent_edges()
graph.write(outputgraph)
