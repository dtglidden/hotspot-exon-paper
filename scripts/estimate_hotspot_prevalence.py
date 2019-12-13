import sys
from os.path import join
from numpy import mean, concatenate, searchsorted
import matplotlib.pyplot as plt
import time
import argparse

chroms = set(['chr'+str(i) for i in range(1,23)] + ['chrY', 'chrX'])

def parse_annotation(anno_path):

	id_to_name, name_to_id = {}, {}
	gene_tx = {}
	with open(anno_path) as in_file:
		for line in in_file:
			if line[0] != '#':
				chrom, _, seq_type, start, end, _, strand, _, info = line.strip().split('\t')
				info_pairs = info.split('; ')[:-1]
				values = set([e.split(' ')[1].strip('\"') for e in info_pairs])
				info_dict = {e.split(' ')[0]:e.split(' ')[1].strip('\"') for e in info_pairs}

				if seq_type == 'exon' and info_dict['transcript_type'] == 'protein_coding' and chrom in chroms and 'appris_principal_1' in values:
					gene_id, gene_name = info_dict['gene_id'].split('.')[0], info_dict['gene_name']
					tx_id = info_dict['transcript_id'].split('.')[0]

					if gene_id not in gene_tx:
						gene_tx[gene_id] = {}
						id_to_name[gene_id], name_to_id[gene_name] = gene_name, gene_id
					if tx_id not in gene_tx[gene_id]:
						gene_tx[gene_id][tx_id] = []

					exon = (chrom, strand, int(start), int(end))
					gene_tx[gene_id][tx_id].append(exon)

	return gene_tx, id_to_name, name_to_id

def parse_HI_scores(HI_path):

	HI_scores = {}
	with open(HI_path) as in_file:
		for line in in_file:
			gene, score = line.strip().split('\t')
			if gene != 'Gene_name':
				HI_scores[gene] = float(score)

	return HI_scores

def parse_usage(sample, usage_dir):

	usages = {ss_type:{chrom:{'+':{}, '-':{}} for chrom in chroms} for ss_type in ['3p', '5p']}
	for ss_type in usages:
		with open(join(usage_dir, '_'.join([sample, ss_type, 'ss_usage.txt']))) as in_file:
			for line in in_file:
				chrom, site, strand, _, _, usage = line.strip().split('\t')
				if chrom != 'Chrom':
					usages[ss_type][chrom][strand][int(site)] = float(usage)

	return usages


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('--anno_path', help='Path to GENCODE GTF annotation file')
	parser.add_argument('--HI_path', help='Path to gene haploinsufficiency scores')
	parser.add_argument('--usage_dir', help='Directory with 5\' and 3\' usage files')
	parser.add_argument('--sample', help='Name of sample to process (should match first portion of usage file name)')
	parser.add_argument('--out_dir', help='Where to output the predicted hotspots file')
	args = parser.parse_args()
	anno_path, HI_path, usage_dir, sample, out_dir = args.anno_path, args.HI_path, args.usage_dir, args.sample, args.out_dir

	gene_tx, id_to_name, name_to_id = parse_annotation(anno_path)
	longest_tx = {g:max(gene_tx[g].keys(), key=lambda t:len(gene_tx[g][t])) for g in gene_tx}

	print(time.strftime('%m-%d %I:%M:%S%p') + ' - Done parsing annotation...')

	HI_scores = parse_HI_scores(HI_path)

	print(time.strftime('%m-%d %I:%M:%S%p') + ' - Done parsing HI scores...')

	print(time.strftime('%m-%d %I:%M:%S%p') + ' - Processing {0}...'.format(sample))

	usages = parse_usage(sample, usage_dir)

	gene_usages = {}
	for gene_id in gene_tx:
		if id_to_name[gene_id] in HI_scores:
			gene_usages[gene_id] = []
			for chrom, strand, start, end in gene_tx[gene_id][longest_tx[gene_id]][1:-1]:
				if strand == '+':
					fivep_site, threep_site = end, start
				else:
					fivep_site, threep_site = start, end
				if fivep_site in usages['5p'][chrom][strand] and threep_site in usages['3p'][chrom][strand]:
					gene_usages[gene_id] += [usages['5p'][chrom][strand][fivep_site], usages['3p'][chrom][strand][threep_site]]


	with open(join(out_dir, '{0}_predicted_hotspot_exons_by_usage.txt'.format(sample)), 'w') as out_file:
		exons_low, exons_reg = 0, 0
		out_file.write('Gene_ID\tChrom\tStart\tEnd\tStrand\t5p_usage\t3p_usage\tLow_usage\tNeighborhood_size\n')
		for gene_id in gene_tx:
			gene_name, intron_num = id_to_name[gene_id], len(gene_tx[gene_id][longest_tx[gene_id]])-1
			if gene_name in HI_scores and intron_num > 1:
				gene_HI = HI_scores[gene_name]
				HI_neighbors = [name_to_id[g] for g in HI_scores if g in name_to_id and abs(HI_scores[g] - gene_HI) <= 0.1]
				intron_and_HI_neighbors = [g for g in HI_neighbors if abs(len(gene_tx[g][longest_tx[g]])-intron_num)<0.1*intron_num and len(gene_tx[g][longest_tx[g]])>2]
				neighbor_usages = [gene_usages[g] for g in intron_and_HI_neighbors]
				if len(neighbor_usages) > 0:
					neighbor_usages = sorted(concatenate(neighbor_usages))
					for chrom, strand, start, end in gene_tx[gene_id][longest_tx[gene_id]][1:-1]:
						if strand == '+':
							fivep_site, threep_site = end, start
						else:
							fivep_site, threep_site = start, end
						if fivep_site in usages['5p'][chrom][strand] and threep_site in usages['3p'][chrom][strand]:
							fivep_usage, threep_usage = usages['5p'][chrom][strand][fivep_site], usages['3p'][chrom][strand][threep_site]
							fivep_rank = (searchsorted(neighbor_usages, fivep_usage, side='right')+1)/float(len(neighbor_usages))
							threep_rank = (searchsorted(neighbor_usages, threep_usage, side='right')+1)/float(len(neighbor_usages))
							output = [gene_id, chrom, start, end, strand, fivep_usage, threep_usage]
							if fivep_rank <= 0.1 and threep_rank <= 0.1:
								output += ['True', len(neighbor_usages)]
								exons_low += 1
							else:
								output += ['False', len(neighbor_usages)]
								exons_reg += 1
							out_file.write('\t'.join([str(e) for e in output]) + '\n')


