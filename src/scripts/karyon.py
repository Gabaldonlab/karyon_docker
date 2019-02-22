#!/usr/bin/python2.7
desc="""Karyon pipeline.
More info at: https://github.com/Gabaldonlab/karyon
"""
epilog="""Author: Miguel Angel Naranjo (email) Barcelona, 01/01/2019"""
import sys, os, re
import argparse
import psutil
import pysam
import string
import random
from spades_recipee import call_SPAdes
from prepare_libraries import preparation
from trimming_libraries import trimming, get_right_path
from varcall_recipee import var_call

def id_generator(size=6, chars=string.ascii_uppercase + string.digits): 
	return ''.join(random.choice(chars) for _ in range(size))

###Parses the config file in order to check the parameters of all the programs.###
def parse_config(config):
	config_dict = {}
	prev = 0
	for line in open(config):
		if line[0] == "#": continue
		elif line[0] == "+":
			prev = line[1:-1]
			config_dict[prev] = ["","",""]
		elif line[0] == "@":
			if config_dict[prev][0] != "": continue
			config_dict[prev][0] = line[1:-1]
		elif line[0] == ">":
			config_dict[prev][1] = config_dict[prev][1] + line[1:-1] + " "
		elif line[0] == "?":
			if config_dict[prev][2] != "": continue
			config_dict[prev][2] = line[1:-1] + " "
	return config_dict

###Selects the main library to use. This is set to accelerate the assembly process and improve the results
def select_champion(fastq):
	parse_dict = {}
	for i in open(fastq):
		chunk = i.split()
		if chunk[5] == "2": continue
		else:
			parse_dict[chunk[0]] = chunk[1:]
	champion=[0,'']
	if args.favourite == False:
		for element in parse_dict:
			if int(parse_dict[element][2]) > champion[0]:
				champion = [int(parse_dict[element][2]), get_right_path(element)]
	else:
		champion = [0, get_right_path(args.favourite)]
	return champion

if __name__ == '__main__':	
	parser = argparse.ArgumentParser(description=desc, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-c', '--configuration', default=False, help="Configuration file. By default will use ./configuration.txt as configuration file")
	parser.add_argument('-l', '--libraries', required=True, nargs='+', help="Fastq libraries to use for assembly and variant calling. Unsuitable libraries for any of the steps will be ignored")
	parser.add_argument('-d', '--output_directory', required=True, help='Directory where all the output files will be generated')
	parser.add_argument('-r', '--reference', default=False, help='Reference assembly. If provided, it will skip the initial assembly step')
	parser.add_argument('-K', '--dirty_kitchen', action='store_true', default=False, help='If this tag is active, the program will not remove all intermediary files in the fodler kitchen after it has finished')
	parser.add_argument('-T', '--no_trimming', action='store_true', default=False, help='If this tag is active, the program will skip the trimming step')
	parser.add_argument('-R', '--no_redundans', action='store_true', default=False, help='If this tag is active, the program will not launch redundans. Remember that redundans is used to perform many downstream analyses. If you skip it, the analyses may not make much sense.')
	parser.add_argument('-V', '--no_varcall', action='store_true', default=False, help="If this tag is active, the program will skip the variant calling step. Many downstream analyses require this and won't be possible if you skip it.")
	parser.add_argument('-P', '--no_plot', action='store_true', default=False, help="If this tag is active, the program will ommit the plots at the end of the skip the variant calling step. Many downstream analyses require this and won't be possible if you skip it.")
	parser.add_argument('-w', '--window_size', default=1000, help='Window size used for some of the analyses. Default is 1000 (1Kb)')
	parser.add_argument('-F', '--favourite', default=False, help='Sets one library as the prefered one for the variant calling analyses. Otherwise, karyon will select the largest library for performing the variant calling protocol.')
	parser.add_argument('-L', '--light', action='store_true', default=False, help='Uses only one library for the assembly process in order to speed the process. If the argument favourite is set, it will use it; otherwise it will select one with the same criteria as the one it uses for the variant calling protocol')
	parser.add_argument('--job_id', default=False, help='Identifier of the intermediate files generated by the different programs. If false, the program will assign a random number')
	parser.add_argument('-o', '--output_name', default=False, help='Prefix name for all the output files. If ommited, it will generate a random string. This random string will be the same as the random identifier for intermediate files.')
	parser.add_argument('-m', '--memory_limit', default=False, help='Memory limit for all the programs set in Gb. By default it will try to use all memory available.')
	parser.add_argument('-M', '--memory_fraction', default=1, help='Percent of total memory to use by all programs. By default it will use all available memory (default=1), but it may be useful to reduce the percent to avoid freezing other tasks of the computer during peaks.')
	parser.add_argument('-n', '--nodes', default=False, help='Number of computation nodes to use. If set a number higher than total, it will use total. If set a number lower than total, it will calculate memory usage based fraction of nodes set respect total existing nodes.')
	parser.add_argument('--max_scaf2plot', default=20, help="Maximum number of scaffolds to plot for scaffold-specific plots")

	args = parser.parse_args()

	###Defines the location of configuration.txt if setting by default###
	config_path = args.configuration
	if not args.configuration:
		selfpath = os.path.dirname(os.path.realpath(sys.argv[0]))
		config_path = selfpath[:selfpath.rfind('/')]
		config_path = selfpath[:selfpath.rfind('/')]+"/configuration.txt"
	
	true_output = os.path.abspath(args.output_directory)
	if true_output[-1] != "/":
		true_output=true_output+"/"

	###Sets RAM usage options###
	total_nodes = n_nodes = psutil.cpu_count()
	if args.nodes and int(args.nodes) < total_nodes:
		n_nodes = int(args.nodes)

	if not args.memory_limit:
		ram_limit = int(psutil.virtual_memory()[0]/1000000000 * float(args.memory_fraction))
		if n_nodes < total_nodes:
			ram_limit = int(psutil.virtual_memory()[0]/1000000000 * float(args.memory_fraction) * (float(n_nodes)/total_nodes))
	else:
		ram_limit = args.memory_limit * int(args.memory_fraction)
	counter = int(args.max_scaf2plot)
	
	###Sets the job ID and the prefix name for the job. If job ID is not user defined, it produces a random 6 character string. If prefix name is not defined, it uses job ID### 
	job_ID = args.job_id if args.job_id else id_generator()
	name = args.output_name if args.output_name else job_ID

	print '###############'
	print 'Config. path: '+str(config_path)
	print "RAM Limit: "+str(ram_limit)+"Gb"
	print "Nodes: "+str(n_nodes)
	print "Job ID: "+str(job_ID)
	print "Job name: "+str(name)
	print '###############'

	config_dict = parse_config(config_path)
	home = config_dict["karyon"][0]
	prepared_libs = home + "kitchen/" + job_ID + "/prepared_libraries.txt"

	###Checks that the output is not a file. If it does not exist, it creates it.###
	if not os.path.isdir(args.output_directory):
		if os.path.isfile == True: 
			message = "Path is a file" #Should raise an exception an exit the program
			sys.stderr.write("\n%s\n\n"%message)
			sys.exit(1)
		else:
			os.mkdir(args.output_directory)
	else :
		os.rmdir(args.output_directory)
		os.mkdir(args.output_directory)
	os.system("mkdir "+ home + "kitchen/"+job_ID)

	###Parses the libraries and checks their parameters for downstream analyses. Also performs trimming.###
	print '###############'
	print 'Preparing libraries'	
	print '###############'
	libs = ''
	for i in args.libraries:
		libs = libs + " " + i
	preparation(libs.split(), 10000, prepared_libs)
	preplibsfile = open(prepared_libs)
	preplibsfile.seek(0)
	for i in preplibsfile: print i, "mandonguilles"

	libs_parsed = ''
	if not args.no_trimming and args.reference == False:
		print '###############'
		print 'Trimmomatic'
		print '###############'
	
		if config_dict['trimmomatic'][1] == '':
			trimmo_commands = ''
		else:
			trimmo_commands = " -c " + config_dict['trimmomatic'][1]
		trimming(prepared_libs, config_dict["trimmomatic"][0], trimmo_commands, home + "kitchen/"+job_ID+"/trimmomatic.job", true_output, False)
		os.system("bash " + home + "kitchen/"+job_ID+"/trimmomatic.job")
	
		for i in os.listdir(args.output_directory):
			if i.find("parsed_") > -1:
				libs_parsed = libs_parsed + " " + true_output + i		
		preparation(libs_parsed.split(), 10000, prepared_libs)

	###Parsing library names, including putting absolute paths #
	libstring = ''
	backstring = ''
	for i in open(prepared_libs):
		chunk = i.split()
		if chunk[5] == "1":
			libstring = libstring + os.path.abspath(chunk[0]) + " " + os.path.abspath(chunk[6]) + " "
		elif chunk[5] == "2": continue
		else: backstring = backstring + os.path.abspath(chunk[0]) + " "
	libstring = libstring + backstring

	champion = select_champion(prepared_libs)

	print '###############'
	print 'Params'
	print '###############'
	print args.window_size
	print true_output+name+".raw.vcf"
	print true_output+"redundans_output/scaffolds.filled.fa"
	print true_output+name+".sorted.bam"
	print true_output+name+".mpileup"
	print champion[-1]
	print config_dict['nQuire'][0]
	print config_dict["KAT"][0]
	print home + "kitchen/"+job_ID+"/"
	print true_output
	print counter
	print '###############'
	
	###Calling spades_recipee.py to generate the assembly job. In the future it should use config file to select the assembly program to use###
	if not args.reference:
		call_SPAdes(prepared_libs, config_dict['SPAdes'][0], true_output, name, config_dict['SPAdes'][1], False, ram_limit, n_nodes, args.favourite, args.light, champion)

	karyonjobfile = open(true_output+name+"_karyon.job", 'a')
	karyonjobfile.write("\n")
	if args.reference == False:
		karyonjobfile.write("python "+config_dict['Redundans'][0]+"redundans.py"+" -f "+true_output+"dipspades/consensus_contigs.fasta -o "+true_output+"redundans_output -i "+libstring+" -t "+str(n_nodes)+" "+config_dict["Redundans"][1])
	else:
		karyonjobfile.write("python "+config_dict['Redundans'][0]+"redundans.py"+" -f "+ args.reference + "-o "+true_output+"redundans_output -i "+libstring+" -t "+str(n_nodes)+" "+config_dict["Redundans"][1])
	karyonjobfile.close()

	#5) Creates job files and calls them
	var_call(prepared_libs, config_dict, true_output, name, args.favourite, home, str(ram_limit), str(n_nodes))
	os.system ("bash "+true_output+name+"_karyon.job")

	#6) We create the plots
	if args.no_redundans == True:
		if args.reference == True:
			reduced_assembly = args.reference
		else:
			reduced_assembly = true_output+"dipspades/consensus_contigs.fasta"
	else: 
		reduced_assembly = true_output+"redundans_output/contigs.fa"

	counter = int(args.max_scaf2plot)
	if args.no_plot == False or args.no_varcall == False:
		from karyonplots import katplot, allplots
		katplot(reduced_assembly, champion[1], config_dict["KAT"][0], true_output)
		allplots(int(args.window_size), 
			true_output+name+".raw.vcf", 
			reduced_assembly, 
			true_output+name+".sorted.bam", 
			true_output+name+".mpileup", 
			os.path.abspath(champion[-1]), 
			config_dict['nQuire'][0], 
			config_dict["KAT"][0], 
			home + "kitchen/"+job_ID+"/", 
			true_output, 
			counter, 
			job_ID, name)

	###We clean the kitchen###
	if args.dirty_kitchen == True:
		existence = open(home + "kitchen/" + job_ID + '/_keep_existing_', 'w')
		existence.close()

	print "Now I'm cleaning the kitchen..."
	if args.dirty_kitchen == True:
		print "...but keeping what you told me..."
	for e in os.listdir(home + "kitchen/"):
		for i in os.listdir(home + "kitchen/"+e):
			if '_keep_existing_' in os.listdir(home + "kitchen/"+e): continue
			else:
				os.remove(home + "kitchen/"+e+"/"+i)
		if '_keep_existing_' in os.listdir(home + "kitchen/"+e): continue
		else: os.rmdir(home + "kitchen/"+e)
	if args.dirty_kitchen == True:
		print "...now kitchen is how you like it!"	
	else:	
		print "...And the kitchen is clean!"



