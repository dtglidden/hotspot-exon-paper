import sys
from os.path import join
from numpy import mean
import numpy as np
import argparse

chroms = set(['chr'+str(i) for i in range(1,23)] + ['chrY', 'chrX', 'chrM'])

if __name__ == '__main__' :

	parser = argparse.ArgumentParser()
	parser.add_argument('--anno_path', help='Path to GENCODE GTF annotation file')
	parser.add_argument('--HI_path', help='Path to gene haploinsufficency score file')
	parser.add_argument('--usage_dir', help='Directory with 5\' and 3\' usage data files')
	parser.add_argument('--sample', help='Name of the sample from which the usage data was calculated')
	parser.add_argument('--out_dir', help='Where to output the average usage file to')
	args = parser.parse_args()
	anno_path, HI_path, usage_dir, sample, out_dir = args.anno_path, args.HI_path, args.usage_dir, args.sample, args.out_dir

	#Parse gene haploinsufficency scores
	HI_scores = {}
	with open(HI_path) as in_file:
		for line in in_file:
			gene, score = line.strip().split('\t')
			if gene != 'Gene_name':
				HI_scores[gene] = float(score)

	#Parse splice site usage data
	usages = {ss_type:{chrom:{'+':{}, '-':{}} for chrom in chroms} for ss_type in ['3p', '5p']}
	for ss_type in usages:
		with open(join(usage_dir, '_'.join([sample, ss_type, 'ss_usage.txt']))) as in_file:
			for line in in_file:
				chrom, site, strand, _, _, usage = line.strip().split('\t')
				if chrom != 'Chrom':
					usages[ss_type][chrom][strand][int(site)] = float(usage)

	#Parse the exons of each transcript and map transcript IDs to gene IDs
	tx_exons = {}
	gene_to_tx = {}
	with open(anno_path) as in_file:
		for line in in_file:
			chrom, _, entry_type, start, end, _, strand, _, info = line.strip().split('\t')
			if entry_type == 'exon':
				info_pairs = info.split('; ')[:-1]
				values = set([e.split(' ')[1].strip('\"') for e in info_pairs])
				info_dict = {e.split(' ')[0]:e.split(' ')[1].strip('\"') for e in info_pairs}
				tx_id, gene_id, gene_name = info_dict['transcript_id'].split('.')[0], info_dict['gene_id'].split('.')[0], info_dict['gene_name']
				exon = (int(start), int(end))
				if 'appris_principal_1' in values:
					if gene_id not in gene_to_tx:
						gene_to_tx[gene_id] = {'name':gene_name, 'tx':[tx_id], 'info':(chrom, strand)}
					else:
						gene_to_tx[gene_id]['tx'].append(tx_id)
					if tx_id not in tx_exons:
						tx_exons[tx_id] = [exon]
					else:
						tx_exons[tx_id].append(exon)

	#Identify the longest transcript in the gene
	gene_max_tx = {gene:max(gene_to_tx[gene]['tx'], key=lambda t:len(tx_exons[t])) for gene in gene_to_tx}
	#Calculate average usage for each gene that has an HI score
	gene_usage_by_HI = []
	for gene in gene_to_tx:
		gene_name = gene_to_tx[gene]['name']
		if gene_name in HI_scores:
			HI = HI_scores[gene_name]
			max_tx = gene_max_tx[gene]
			exon_num = len(tx_exons[max_tx])
			avg_usage = []
			chrom, strand = gene_to_tx[gene]['info']

			for i in range(len(tx_exons[max_tx])):
				exon_start, exon_end = tx_exons[max_tx][i]
				if strand == '+':
					fivep_site, threep_site = exon_end, exon_start
				else:
					fivep_site, threep_site = exon_start, exon_end
				if i != 0 and threep_site in usages['3p'][chrom][strand]:
					avg_usage.append(usages['3p'][chrom][strand][threep_site])
				if i != len(tx_exons[max_tx]) - 1 and fivep_site in usages['5p'][chrom][strand]:
					avg_usage.append(usages['5p'][chrom][strand][fivep_site])

			#Only output data for transcripts that have usage for each splice site in every exon
			if len(avg_usage) > max(2.0*len(tx_exons[max_tx])-3, 0):
				gene_usage_by_HI.append((gene, chrom, strand, str(mean(avg_usage)), str(HI)))

	with open(join(out_dir, '{}_gene_usage_HI_primary_tx.txt'.format(sample)), 'w') as out_file:
			out_file.write('Ensembl_ID\tChrom\tStrand\tAvg_usage\tHI_score\n')
			for i in range(len(gene_usage_by_HI)):
				out_file.write('\t'.join(gene_usage_by_HI[i]) + '\n')





















