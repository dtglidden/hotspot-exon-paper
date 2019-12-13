import sys
from os.path import join
from numpy import mean
import argparse

chroms = set(['chr'+str(i) for i in range(1,23)] + ['chrY', 'chrX', 'chrM'])

if __name__ == '__main__' :

	parser = argparse.ArgumentParser()
	parser.add_argument('--anno_path', help='Path to GENCODE GTF annotation file')
	parser.add_argument('--intron_path', help='Path to a BED file of introns')
	parser.add_argument('--usage_dir', help='Directory with 5\' and 3\' usage data files')
	parser.add_argument('--sample', help='Name of the sample from which the usage data was calculated')
	parser.add_argument('--out_dir', help='Where to output the average usage file to')
	args = parser.parse_args()
	anno_path, intron_path, usage_dir, sample, out_dir = args.anno_path, args.intron_path, args.usage_dir, args.sample, args.out_dir

	#Parse introns from BED file of GENCODE annotation obtained from UCSC table browser
	tx_introns = {}
	with open(intron_path) as in_file:
		for line in in_file:
			chrom, start, end, info, _, strand = line.strip().split('\t')
			if chrom in chroms:
				tx_id, intron = info.split('.')[0], (int(start), int(end)+1)
				if tx_id not in tx_introns:
					tx_introns[tx_id] = [intron]
				else:
					tx_introns[tx_id].append(intron)


	#Map gene IDs to the transcript IDs of each intron set
	gene_to_tx = {}
	with open(anno_path) as in_file:
		for line in in_file:
			chrom, _, entry_type, start, end, _, strand, _, info = line.strip().split('\t')
			if entry_type == 'transcript' and chrom in chroms:
				info_pairs = info.split('; ')[:-1]
				values = set([e.split(' ')[1].strip('\"') for e in info_pairs])
				info_dict = {e.split(' ')[0]:e.split(' ')[1].strip('\"') for e in info_pairs}
				gene_id, tx_id = info_dict['gene_id'].split('.')[0], info_dict['transcript_id'].split('.')[0]
				if tx_id in tx_introns and 'appris_principal_1' in values:
					if gene_id not in gene_to_tx:
						gene_to_tx[gene_id] = {'tx':[tx_id], 'info':(chrom, strand)}
					else:
						gene_to_tx[gene_id]['tx'].append(tx_id)


	#Identify the longest transcript in the gene
	gene_longest_tx = {gene:max(gene_to_tx[gene]['tx'], key=lambda t:len(tx_introns[t])) for gene in gene_to_tx}

	#Parse splice site usage data
	usages = {ss_type:{chrom:{'+':{}, '-':{}} for chrom in chroms} for ss_type in ['3p', '5p']}
	for ss_type in usages:
		with open(join(usage_dir, '_'.join([sample, ss_type, 'ss_usage.txt']))) as in_file:
			for line in in_file:
				chrom, site, strand, _, _, usage = line.strip().split('\t')
				if chrom != 'Chrom':
					usages[ss_type][chrom][strand][int(site)] = float(usage)

	#Calculate average usage for the gene's longest transcript
	gene_usage_by_introns = []
	for gene in gene_longest_tx:
		max_tx = gene_longest_tx[gene]
		intron_num = len(tx_introns[max_tx])
		avg_usage = []
		chrom, strand = gene_to_tx[gene]['info']

		for intron_start, intron_end in tx_introns[max_tx]:

			if strand == '+':
				fivep_site, threep_site = intron_start, intron_end
			else:
				fivep_site, threep_site = intron_end, intron_start

			if fivep_site in usages['5p'][chrom][strand] and threep_site in usages['3p'][chrom][strand]:
				avg_usage.append(usages['5p'][chrom][strand][fivep_site])
				avg_usage.append(usages['3p'][chrom][strand][threep_site])

		#Only output data for transcripts that have usage for each splice site in every intron
		if len(avg_usage) >= 2.0*len(tx_introns[max_tx]):
			gene_usage_by_introns.append((gene, chrom, strand, str(mean(avg_usage)), str(intron_num)))


	with open(join(out_dir, '{}_gene_usage_intron_primary_tx.txt'.format(sample)), 'w') as out_file:
		out_file.write('Ensembl_ID\tChrom\tStrand\tAvg_usage\tIntron_num\n')
		for i in range(len(gene_usage_by_introns)):
			out_file.write('\t'.join(gene_usage_by_introns[i]) + '\n')











