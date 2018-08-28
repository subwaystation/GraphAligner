#!/usr/bin/python

import fileinput

sequence = ""

for line in fileinput.input():
	if line[0] == '@' or line[0] == '>':
		if len(sequence) > 0: print(sequence)
		print('>' + line[1:].strip().replace(' ', '_'))
		sequence = ""
	else:
		sequence += line.strip()

if len(sequence) > 0: print(sequence)
