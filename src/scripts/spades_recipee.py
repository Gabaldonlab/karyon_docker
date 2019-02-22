#!/bin/python
import sys, numpy, os.path, re
import argparse
from Bio import SeqIO
from trimming_libraries import get_right_path

### Generates a string that call to SPADes or dipSPADes

def call_SPAdes(library_file, path, output, name, commands, no_diploid, memory_limit, nodes, favourite, light, champion):
	libstring = ' '
	backstring = ''		
	for i in open(library_file):
		if favourite != False and i.find(champion[-1][champion[-1].rfind("/"):]) > -1: continue
		chunk = i.split()
		if light == True and i.find(champion[-1][champion[-1].rfind("/"):]) == -1: continue
		elif chunk[5] == "1":
			libstring = libstring + "-1 " + get_right_path(chunk[0]) + " -2 " + get_right_path(chunk[6]) + " "
		elif chunk[5] == "2": continue
		else: backstring = backstring + get_right_path(chunk[5]) + " " + get_right_path(chunk[0]) + " "
		libstring = libstring + backstring
	print libstring, "aaaaaaaaaaaaaaaaaaaaaaaaaa"
	
	outputfile = open(output+name+"_karyon.job", 'w')
	if no_diploid == True:
		outputfile.write("python " + path + "bin/spades.py" + libstring + " -t " + str(nodes) + " -m " +  str(memory_limit) + " " +commands + "-o " + output)
	else:
		outputfile.write("python " + path + "bin/dipspades.py" + libstring + " -t " + str(nodes) + " -m " +  str(memory_limit) + " " +commands + "-o " + output)
	print "los kibis"
	print "python " + path + "bin/spades.py" + libstring + " -t " + str(nodes) + " -m " +  str(memory_limit) + " " +commands + "-o " + output
	outputfile.write("\n")
	outputfile.close()

