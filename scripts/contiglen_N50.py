#!/usr/bin/python

import fileinput

lens = []
lensum = 0

for line in fileinput.input():
	lens.append(int(line))
	lensum += int(line)

lens.sort(key=lambda x: -x)

partialsum = 0
for i in range(0, len(lens)):
	partialsum += lens[i]
	if partialsum > lensum / 2: 
		print(lens[i])
		break
