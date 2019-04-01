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
import sys, os, re, subprocess, math
import argparse
import psutil
import pysam
from Bio import SeqIO
import numpy as np
import numpy.random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns
import pandas as pd
import scipy.stats
from decimal import Decimal
from scipy import stats

def launch_nQuire(bam, nQuire, kitchen, bam_file):
	BAMtemp = pysam.AlignmentFile(kitchen+"BAMtemp.bam", 'wb', template=bam_file)
	for i in bam:
		BAMtemp.write(i)
	BAMtemp.close()
	pysam.index(kitchen+"BAMtemp.bam")
	
	os.system(nQuire+" create -b "+kitchen+"BAMtemp.bam -o "+kitchen+"nQuire_temp")
	os.system(nQuire+" lrdmodel "+kitchen+"nQuire_temp.bin > "+kitchen+"nQuire_temp.report")
	nQuire_report = kitchen+"nQuire_temp.report"
		
	free_score, diplo_score, triplo_score, tetra_score = 0.0, 0.0, 0.0, 0.0
	for line in open(nQuire_report):
		if line.find("free") == 0: free_score = float(line.split()[1])
		elif line.find("dipl") == 0: diplo_score = float(line.split()[1])
		elif line.find("trip") == 0: triplo_score = float(line.split()[1])
		elif line.find("tetr") == 0: tetra_score = float(line.split()[1])
		else: continue
	os.remove(kitchen+"BAMtemp.bam")
	os.remove(kitchen+"BAMtemp.bam.bai")
	os.remove(kitchen+"nQuire_temp.bin")
	if free_score < 0.001:
		free_score, diplo_score, triplo_Score, tetra_Score = float('nan'), float('nan'), float('nan'), float('nan')
	return round(diplo_score/free_score, 3) , round(triplo_score/free_score, 3), round(tetra_score/free_score, 3)
	
def ttest_ploidy(number_list):
	y_dip, y_tripA, y_tripB, y_tetraA, y_tetraB = [], [], [], [], []
	for i in range(len(number_list)):
		y_dip.append(0)
		y_tripA.append(0.301)
		y_tripB.append(-0.301)
		y_tetraA.append(0.477)
		y_tetraB.append(-0.477)
	mean_number_list = numpy.nanmean(number_list)
	R2_diploid = (stats.ttest_1samp(number_list, math.log(1.0,2),nan_policy='omit'))
	R2_triploidA, R2_triploidB  = (stats.ttest_1samp(number_list, math.log(2.0,2),nan_policy='omit')), (stats.ttest_1samp(number_list, math.log(2.0,2)*-1,nan_policy='omit'))
	R2_tetraploidA, R2_tetraploidB  = (stats.ttest_1samp(number_list, math.log(3.0,2),nan_policy='omit')), (stats.ttest_1samp(number_list, math.log(1/3.0,2),nan_policy='omit'))
	return R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB#, mean_number_list, numpy.nanstd(number_list)

def window_walker(window_size, step, vcf_file, fasta_file, bam_file, nQuire, kitchen, newpath, counter):
	window_stats = []
	prev_record_name = False
	fasta = SeqIO.parse(fasta_file, "fasta")
	count = 0
	for record in fasta:
		count = count + 1
		if count > counter : break
		if record.name != prev_record_name and prev_record_name != False:
			nQuire_plot(window_stats, window_size, newpath, bam_file)
			window_stats = []
		if prev_record_name == False:
			prev_record_name = record.name
		start, end = 0, len(record)
		log_refalt_list = []
		while start + step <= end:
			window = record.name+":"+str(start)+":"+str(start+window_size)
			VCF = vcf_file.fetch(record.name, start, start + window_size) 
			mean_refaltcov_list = []
			for i in VCF:
				refalt = (str(i).split()[-1].split(':')[1].split(","))
				snp_pos = int(str(i).split()[1])
				ref, alt = float(refalt[0]), float(refalt[1])
				totalcov = numpy.nanmean(bam_file.count_coverage(record.name, snp_pos, snp_pos+1, quality_threshold=0))
				mean_refaltcov_list.append(totalcov)
				if alt == 0 or ref == 0:
					value = float('nan')
				else:
					value = alt/(float(ref)+0.01)
				log_refalt_list.append(math.log(value,2))

			BAMtemp = bam_file.fetch(record.name, start, start + window_size)
			mean_cov = numpy.nanmean(bam_file.count_coverage(record.name, start, start + window_size, quality_threshold=0))

			stdev_cov = numpy.nanstd(bam_file.count_coverage(record.name, start, start + window_size))
			diplo_score, triplo_score, tetra_score = launch_nQuire(BAMtemp, nQuire, kitchen, bam_file)
			R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB = ttest_ploidy(log_refalt_list)

			window_stats.append([window, start+window_size/2, diplo_score, triplo_score, tetra_score, R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB, mean_cov, stdev_cov, mean_refaltcov_list])
			start = start + step
			vcf_file.seek(0)
	return window_stats

def nQuire_plot(value_list, window_size, newpath, bam_file):
	name, x, y1, y2, y3, std_cov, mean_cov, snp_den = '', [], [], [], [], [], [], []
	all_refalt_list, pos_list = [], []
	for i in value_list:
		if i[2] > 0:
			name = i[0].split(":")[0]
			x.append(i[1]) 
			if i[2] > 1: y1.append(1.0)
			else: y1.append(i[2])
			if i[3] > 1: y2.append(1.0)
			else: y2.append(i[3])
			if i[4] > 1: y3.append(1.0)
			else: y3.append(i[4])
			std_cov.append(i[-1])
			mean_cov.append(i[-3])
			snp_den.append(i[-2])
			for e in i[-1]:
				pos_list.append(i[1])
				all_refalt_list.append(e)
		else: continue
	if len(x) > 0:
		sns.set(style="darkgrid")
	
		fig, (p0,p1,p2,p3,p4) = plt.subplots(nrows=5, sharex=True, figsize=(45,15))
		plt.subplot(5,1,1)
		plt.xlim(0,x[-1])
		plt.plot(x, y1, 'ro')
		plt.plot(x, y1, color='#aa0000', linestyle='--')
		plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
		plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

		plt.subplot(5,1,2)
		plt.xlim(0,x[-1])
		plt.plot(x, y2, 'go')
		plt.plot(x, y2, color='#00aa00', linestyle='--')
		plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
		plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

		plt.subplot(5,1,3)
		plt.xlim(0,x[-1])
		plt.plot(x, y3, 'bo')
		plt.plot(x, y3, color='#0000aa', linestyle='--')
		plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
		plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

		plt.subplot(5,1,4)
		plt.xlim(0,x[-1])
		plt.plot(x, mean_cov, color='grey', linestyle='-')
		
		plt.subplot(5,1,5)
		plt.xlim(0,x[-1])
		plt.ylim(0,200)
		xy = np.vstack([pos_list,all_refalt_list])

		z = gaussian_kde(xy)(xy)
		plt.scatter(pos_list, all_refalt_list, c=z, s=30, edgecolor='')
		plt.savefig(newpath+name+".png")
		print newpath+name+".png has been created"
		plt.clf()

