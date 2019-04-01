import sys, os, re
import argparse
from nQuire_plot import *
from karyonplots import *

parser = argparse.ArgumentParser(description='Process some integers.')

if __name__ == '__main__':	
	parser.add_argument('-w', '--window_size', default=1000, help="window size")
	parser.add_argument('-s', '--step', default=1000, help='sliding window size')
	parser.add_argument('-v', '--vcf', required=True, help='vcf file')
	parser.add_argument('-f', '--fasta', required=True, help='reference genome in fasta file')
	parser.add_argument('-b', '--bam', required=True, help='bam file')
	parser.add_argument('-p', '--pileup', required=True, help='pileup file')
	parser.add_argument('--nQuire', default="/home/mnaranjo//users/tg/mnaranjo/scripts_and_stuff/nQuire/nQuire", help='nQuire path')
	args = parser.parse_args()

cwd = os.getcwd()

scaflist = set()
pileup = open(args.pileup)
for i in pileup:
		scaflist.add(i.split()[0])
pileup.seek(0)
VCF = pysam.VariantFile(args.vcf+".gz", 'r')
bam_file = pysam.AlignmentFile(args.bam, 'rb')



window_walker(args.window_size, args.step, VCF, args.fasta, bam_file, args.nQuire, cwd+"/", cwd+"/", 10)
snp_density = extract_vcf_data(open(args.vcf), args.window_size, 20)
var_v_cov(args.vcf, args.pileup, snp_density, args.window_size, cwd+"/")
var_v_cov_per_scaf(args.vcf, args.pileup, snp_density, scaflist, args.window_size, cwd+"/", 10)
