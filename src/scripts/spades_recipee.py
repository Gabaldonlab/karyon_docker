#!/bin/python
import sys, numpy, os.path, re
import argparse
from Bio import SeqIO

### Generates a string that call to SPADes or dipSPADes

def call_SPAdes(library_file, path, output, name, commands, no_diploid, memory_limit, nodes, champion, light):
	libstring = ' '
	backstring = ''		
	for i in open(library_file):
		if i.find(champion[-1][champion[-1].rfind("/"):]) > -1: continue
		chunk = i.split()
		if chunk[5] == "1":
			libstring = libstring + "-1 " + os.path.abspath(chunk[0]) + " -2 " + os.path.abspath(chunk[6]) + " "
		elif chunk[5] == "2": continue
		else: backstring = backstring + os.path.abspath(chunk[5]) + " " + os.path.abspath(chunk[0]) + " "
	libstring = libstring + backstring
	
	outputfile = open(output+name+"_karyon.job", 'w')
	if no_diploid == True:
		outputfile.write("python " + path + "bin/spades.py" + libstring + " -t " + str(nodes) + " -m " +  str(memory_limit) + " " +commands + "-o " + output)
	else:
		outputfile.write("python " + path + "bin/dipspades.py" + libstring + " -t " + str(nodes) + " -m " +  str(memory_limit) + " " +commands + "-o " + output)
	outputfile.write("\n")
	outputfile.close()

