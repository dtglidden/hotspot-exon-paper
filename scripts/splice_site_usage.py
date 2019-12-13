import sys
from os.path import join
from intervaltree import Interval, IntervalTree
import numpy as np
import time
import argparse

site_types = ['3p', '5p']
chroms = set(['chr'+str(i) for i in range(1,23)] + ['chrY', 'chrX'])

def parse_annotation(anno_path):

	print(time.strftime('%m-%d %I:%M:%S%p') + ' - Parsing GENCODE annotation...')

	anno_sites = {ss_type:{chrom:{'+':[], '-':[]} for chrom in chroms} for ss_type in site_types}
	with open(anno_path) as in_file:
		for line in in_file:
			if line[0] != '#':
				chrom, _, entry_type, start, end, _, strand, _, info = line.strip().split('\t')

				if entry_type == 'exon' and chrom in chroms:
					start, end = int(start), int(end)
					if strand == '+':
						three_p_ss, five_p_ss = start, end
					else:
						three_p_ss, five_p_ss = end, start

					anno_sites['3p'][chrom][strand].append(three_p_ss)
					anno_sites['5p'][chrom][strand].append(five_p_ss)

	anno_sites = {ss_type:{chrom:{strand:np.sort(anno_sites[ss_type][chrom][strand]) for strand in ['+', '-']} for chrom in chroms} for ss_type in site_types}
	return anno_sites

def ss_usage(sj_paths, anno_sites):

	print(time.strftime('%m-%d %I:%M:%S%p') + ' - Parsing junction reads...')

	strand_dict = {'1':'+', '2':'-'}
	junctions = {chrom:{'+':IntervalTree(), '-':IntervalTree()} for chrom in chroms}

	for sj_path in sj_paths:
		with open(sj_path) as in_file:

			for line in in_file:

				chrom, start, end, strand, _, _, reads, _, _ = line.strip().split('\t')
				start, end, reads = int(start)-1, int(end)+1, int(reads)

				if chrom not in chroms or strand == '0' or reads == 0:
					continue

				strand = strand_dict[strand]
				if strand == '+':
					five_p_ss, three_p_ss = start, end
				else:
					five_p_ss, three_p_ss = end, start

				if five_p_ss not in anno_sites['5p'][chrom][strand] or three_p_ss not in anno_sites['3p'][chrom][strand]:
					continue

				junc = Interval(start, end, {'reads':reads})
				junctions[chrom][strand].add(junc)

	print(time.strftime('%m-%d %I:%M:%S%p') + ' - Calculating splice site usages...')

	fivep_ss_usage, threep_ss_usage = {chrom:{'+':{}, '-':{}} for chrom in chroms}, {chrom:{'+':{}, '-':{}} for chrom in chroms}
	for chrom in junctions:
		for strand in junctions[chrom]:
			juncs_processed = set()
			for junc_int in junctions[chrom][strand]:
				if (junc_int.begin, junc_int.end) not in juncs_processed:
					if strand == '+':
						fivep_ss, threep_ss = junc_int.begin, junc_int.end
						next_threep_index = np.searchsorted(anno_sites['3p'][chrom][strand], fivep_ss, side='right')
						next_fivep_index = np.searchsorted(anno_sites['5p'][chrom][strand], threep_ss, side='left')-1
						next_fivep_ss = anno_sites['5p'][chrom][strand][next_fivep_index]
						next_threep_ss = anno_sites['3p'][chrom][strand][next_threep_index]
						fivep_juncs = junctions[chrom][strand].overlap(fivep_ss, next_threep_ss)
						threep_juncs = junctions[chrom][strand].overlap(next_fivep_ss, threep_ss)
					else:
						fivep_ss, threep_ss = junc_int.end, junc_int.begin
						next_threep_index = np.searchsorted(anno_sites['3p'][chrom][strand], fivep_ss, side='left')-1
						next_fivep_index = np.searchsorted(anno_sites['5p'][chrom][strand], threep_ss, side='right')
						next_fivep_ss = anno_sites['5p'][chrom][strand][next_fivep_index]
						next_threep_ss = anno_sites['3p'][chrom][strand][next_threep_index]
						fivep_juncs = junctions[chrom][strand].overlap(next_threep_ss, fivep_ss)
						threep_juncs = junctions[chrom][strand].overlap(threep_ss, next_fivep_ss)

					fivep_ss_reads, fivep_reads_total = 0.0, 0.0
					for junc in fivep_juncs:
						if strand == '+':
							site = junc.begin
						else:
							site = junc.end
						if site == fivep_ss:
							fivep_ss_reads += junc.data['reads']
						fivep_reads_total += junc.data['reads']
					fivep_usage = fivep_ss_reads/fivep_reads_total

					threep_ss_reads, threep_reads_total = 0.0, 0.0
					for junc in threep_juncs:
						if strand == '+':
							site = junc.end
						else:
							site = junc.begin
						if site == threep_ss:
							threep_ss_reads += junc.data['reads']
						threep_reads_total += junc.data['reads']
					threep_usage = threep_ss_reads/threep_reads_total

					fivep_ss_usage[chrom][strand][fivep_ss] = (fivep_ss_reads, fivep_reads_total-fivep_ss_reads)
					threep_ss_usage[chrom][strand][threep_ss] = (threep_ss_reads, threep_reads_total-threep_ss_reads)
					juncs_processed.add((junc_int.begin, junc_int.end))

	return fivep_ss_usage, threep_ss_usage

if __name__ == '__main__' :

	parser = argparse.ArgumentParser()
	parser.add_argument('--anno_path', help='Path to GENCODE GTF annotation file')
	parser.add_argument('--sj_files', help='Comma-separated list of splice junction read files from STAR')
	parser.add_argument('--sample', help='Name of the sample to be used in output file names')
	parser.add_argument('--out_dir', help='Where to output the splice site usage files to')
	args = parser.parse_args()
	anno_path, sj_files, sample, out_dir = args.anno_path, args.sj_files.split(','), args.sample, args.out_dir

	anno_sites = parse_annotation(anno_path)

	print(time.strftime('%m-%d %I:%M:%S%p') + ' - Processing usage data for ' + sample + '...')

	fivep_ss_usage, threep_ss_usage = ss_usage(sj_files, anno_sites)

	with open(join(out_dir, sample + '_5p_ss_usage.txt'), 'w') as out_file:
		out_file.write('Chrom\tSite\tStrand\tInclusion_reads\tExclusion_reads\tUsage\n')
		for chrom in fivep_ss_usage:
			for strand in fivep_ss_usage[chrom]:
				for fivep_site in fivep_ss_usage[chrom][strand]:
					inc_reads, exc_reads = fivep_ss_usage[chrom][strand][fivep_site]
					usage = inc_reads/float(inc_reads + exc_reads)
					out_file.write('\t'.join([chrom, str(fivep_site), strand, str(inc_reads), str(exc_reads), str(usage)]) + '\n')

	with open(join(out_dir, sample + '_3p_ss_usage.txt'), 'w') as out_file:
		out_file.write('Chrom\tSite\tStrand\tInclusion_reads\tExclusion_reads\tUsage\n')
		for chrom in threep_ss_usage:
			for strand in threep_ss_usage[chrom]:
				for threep_site in threep_ss_usage[chrom][strand]:
					inc_reads, exc_reads = threep_ss_usage[chrom][strand][threep_site]
					usage = inc_reads/float(inc_reads + exc_reads)
					out_file.write('\t'.join([chrom, str(threep_site), strand, str(inc_reads), str(exc_reads), str(usage)]) + '\n')






