#!/usr/bin/python

import sys

infile = sys.argv[1]
freqcutoff = int(sys.argv[2])
outfile = sys.argv[3]

lencutoff = 110

printing = False

with open(infile) as inf:
	with open(outfile, 'w') as outf:
		for line in inf:
			if line[0] == '>':
				freq = float(line.split(' ')[3].split(':')[2])
				length = int(line.split(' ')[1].split(':')[2])
				printing = freq >= freqcutoff or length >= lencutoff
			if printing: outf.write(line)
