#!/usr/bin/python

import fileinput

print("node\tlength\treads\tcoverage\tchain")

for l in fileinput.input():
	if l[0] != 'S': continue
	parts = l.strip().split('\t')
	nodeid = parts[1]
	length = 0
	reads = 0
	coverage = 0.0
	chain = ""
	for part in parts[3:]:
		if part[0:5] == "LN:i:":
			length = part[5:]
		elif part[0:5] == "RC:i:":
			reads = part[5:]
		elif part[0:5] == "bc:Z:":
			chain = part[5:]
	if length > 0: coverage = str(float(reads) / float(length))
	print(nodeid + "\t" + length + "\t" + reads + "\t" + coverage + "\t" + chain)
