#!/bin/python
import sys, numpy, os.path, re
import argparse
from Bio import SeqIO


### Selects a single library to use in the mapping step. By default, it selects the library with higher coverage ###
def select_champion(fastq, favourite):
	parse_dict = {}
	for i in open(fastq):
		chunk = i.split()
		if chunk[5] == "2": continue
		else:
			parse_dict[chunk[0]] = chunk[1:]
	champion=[0,'']		
	if favourite == False:
		for element in parse_dict:
			if int(parse_dict[element][2]) > champion[0]:
				champion = [int(parse_dict[element][2]), element]
	else:
		champion = [0,args.favourite]
	return champion, parse_dict

### Converts libraries with phred33 scores to phred33 in order to unify all quality scores ###
def job_description (fastqlist):
	phred64dict = {}
	for element in fastqlist:
		switch = False
		for line in open(element):
			if line[0] == "@": continue
			else:
				if switch == False:
					if line[0] == "+":
						switch = True
					else: continue
				if switch == True:
					if line.find("0") > -1:
						phred64dict[element] = False
						break
				else:
					switch = False
	for element in fastqlist:
		if element not in phred64dict:
			phred64dict[element] = True
	for library in fastq:
		if phred64dict[library] == True:
			SeqIO.convert(library, 'fastq-illumina', library+".converted.fq", 'fastq-sanger')
	return phred64dict


### This function is used to try to guess the complementary library of a given one, mostly based on its name and common naming conventions for forward and reverse pairs ###

def create_hypo_dict(fastq):
	hypo_dict = {}
	for element in fastq:
		for m in re.finditer('2', element):
			hypothetical = element[:m.start()]+"1"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				fastq.remove(element)
				break
		for m in re.finditer('R', element):
			hypothetical = element[:m.start()]+"F"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				fastq.remove(element)
				break
		for m in re.finditer('rev', element):
			hypothetical = element[:m.start()]+"fwd"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				fastq.remove(element)
				break
	return hypo_dict


### Generates a file with all the variant calling protocol ###

def var_call(fastq, config_dict, output, name, favourite, home, memory, nodes, no_reduction):
	outputfile = output+name+"_karyon.job"
	parse_dict = {}
	libstring = ' '
	backstring = ''
	for i in open(fastq):
		chunk = i.split()
		if chunk[5] == "2": continue
		parse_dict[chunk[0]] = chunk[1:]
		if chunk[5] == "1":
			libstring = libstring + os.path.abspath(chunk[0]) + " " + os.path.abspath(chunk[6]) + " "	
		elif chunk[5] == 's' : continue
		else: backstring = backstring + os.path.abspath(chunk[5]) + " " + os.path.abspath(chunk[0]) + " "
	libstring = libstring + backstring
	
	bash_job = open(outputfile, "a")
	locspp = output + name
	pairs = ''
	bash_job.write("\n")
	if no_reduction == False:
		bash_job.write("cp "+output+"redundans_output/scaffolds.filled.fa "+locspp+".fasta\n\n")
	else:
		bash_job.write("cp "+no_reduction+" "+locspp+".fasta\n\n")
	bash_job.write(config_dict["BWA"][0] + "bwa index "+locspp+".fasta\n\n")
	
	champion, parse_dict = select_champion(fastq, favourite)
	bash_job.write(config_dict["GATK"][0] + "gatk CreateSequenceDictionary -R "+locspp+'.fasta -O '+locspp+'.dict\n\n')
	bash_job.write(config_dict["BWA"][0]+"bwa index "+locspp+'.fasta\n\n')
	if parse_dict[champion[1]][4] == '1':
		if config_dict['BWA'][0] == '':
			bash_job.write("python " + home + "scripts/launch_bwa.py -r "+locspp+".fasta -f1 "+os.path.abspath(champion[1])+" -f2 "+os.path.abspath(parse_dict[champion[1]][5])+" -n "+locspp+"\n\n")	
		else:
			bash_job.write("python " + home + "scripts/launch_bwa.py -r "+locspp+".fasta -f1 "+os.path.abspath(champion[1])+" -f2 "+os.path.abspath(parse_dict[champion[1]][5])+" -n "+locspp+" -b "+config_dict['BWA'][0]+"\n\n")	
	if parse_dict[champion[1]][4] == 's':
		if config_dict['BWA'][0] == '':
			bash_job.write("python " + home + "scripts/launch_bwa.py -r "+locspp+".fasta -f1 "+os.path.abspath(champion[1])+" -f2 "+os.path.abspath(parse_dict[champion[1]][5])+" -n "+locspp+"\n\n")	
		else:
			bash_job.write("python " + home + "scripts/launch_bwa.py -r "+locspp+".fasta -f1 "+os.path.abspath(champion[1])+" -n "+locspp+" -b "+parse_dict['BWA'][0]+"\n\n")
	bash_job.write(config_dict["GATK"][0] + "gatk --java-options -Xmx"+memory+"G MarkDuplicates -I "+locspp+'.sorted.bam -O '+locspp+'.marked.sorted.bam  -M '+locspp+'.markedstats.txt ; mv '+locspp+'.marked.sorted.bam '+ locspp+'.sorted.bam\n')	
	bash_job.write(config_dict["samtools"][0]+"samtools index "+locspp+'.sorted.bam\n')
	bash_job.write(config_dict["samtools"][0]+"samtools faidx "+locspp+'.fasta\n')
	bash_job.write(config_dict["GATK"][0] + "gatk --java-options -Xmx"+memory+"G HaplotypeCaller -R "+locspp+'.fasta -I '+locspp+'.sorted.bam -O '+locspp+'.raw.vcf\n\n')
	bash_job.write(config_dict["GATK"][0] + "gatk --java-options -Xmx"+memory+"G VariantFiltration -V "+locspp+'.raw.vcf -O '+ locspp+'.vcf ; rm '+locspp+'.raw.vcf\n')
	bash_job.write(' rm '+locspp+'.bam ; rm '+locspp+'.sam\n')
	bash_job.write(config_dict["GATK"][0] + "gatk --java-options -Xmx"+memory+"G CollectMultipleMetrics -I "+locspp+'.sorted.bam -O '+locspp+'_diagnostics --PROGRAM CollectAlignmentSummaryMetrics --PROGRAM CollectGcBiasMetrics --PROGRAM CollectInsertSizeMetrics\n')
	bash_job.write(config_dict["samtools"][0]+"samtools mpileup -f " + locspp + ".fasta "+locspp+'.sorted.bam > '+locspp+'.mpileup\n')	
