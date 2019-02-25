import sys
from os.path import join

def ss_usage(sj_paths):

	valid_chroms = set(['chr'+str(i) for i in range(1,23)] + ['chrY', 'chrX'])
	strand_dict = {'1':'+', '2':'-'}
	five_to_three = {chrom:{'+':{}, '-':{}} for chrom in valid_chroms}
	three_to_five = {chrom:{'+':{}, '-':{}} for chrom in valid_chroms}

	for sj_path in sj_paths:
		with open(sj_path) as in_file:

			for line in in_file:

				chrom, start, end, strand, _, _, reads, _, _ = line.strip().split('\t')
				start, end, reads = int(start)-1, int(end)+1, int(reads)

				if chrom not in valid_chroms or strand == '0' or reads == 0:
					continue

				strand = strand_dict[strand]
				if strand == '+':
					five_p_ss, three_p_ss = start, end
				else:
					five_p_ss, three_p_ss = end, start

				if five_p_ss in five_to_three[chrom][strand]:
					if three_p_ss in five_to_three[chrom][strand][five_p_ss]:
						five_to_three[chrom][strand][five_p_ss][three_p_ss] += reads
					else:
						five_to_three[chrom][strand][five_p_ss][three_p_ss] = reads
				else:
					five_to_three[chrom][strand][five_p_ss] = {three_p_ss:reads}

				if three_p_ss in three_to_five[chrom][strand]:
					if five_p_ss in three_to_five[chrom][strand][three_p_ss]:
						three_to_five[chrom][strand][three_p_ss][five_p_ss] += reads
					else:
						three_to_five[chrom][strand][three_p_ss][five_p_ss] = reads
				else:
					three_to_five[chrom][strand][three_p_ss] = {five_p_ss:reads}

	fivep_ss_usage, threep_ss_usage = {chrom:{'+':{}, '-':{}} for chrom in valid_chroms}, {chrom:{'+':{}, '-':{}} for chrom in valid_chroms}
	for chrom in five_to_three:
		for strand in five_to_three[chrom]:
			for five_p_ss in five_to_three[chrom][strand]:
				three_p_sites = five_to_three[chrom][strand][five_p_ss].keys()
				max_three_p = max(three_p_sites, key=lambda tp:five_to_three[chrom][strand][five_p_ss][tp])
				five_p_sites_max_three = three_to_five[chrom][strand][max_three_p].keys()
				five_p_usage = three_to_five[chrom][strand][max_three_p][five_p_ss]/float(sum([three_to_five[chrom][strand][max_three_p][fp] for fp in five_p_sites_max_three]))
				fivep_ss_usage[chrom][strand][five_p_ss] = five_p_usage
	for chrom in three_to_five:
		for strand in three_to_five[chrom]:
			for three_p_ss in three_to_five[chrom][strand]:
				five_p_sites = three_to_five[chrom][strand][three_p_ss].keys()
				max_five_p = max(five_p_sites, key=lambda fp:three_to_five[chrom][strand][three_p_ss][fp])
				three_p_sites_max_five = five_to_three[chrom][strand][max_five_p].keys()
				three_p_usage = five_to_three[chrom][strand][max_five_p][three_p_ss]/float(sum([five_to_three[chrom][strand][max_five_p][tp] for tp in three_p_sites_max_five]))
				threep_ss_usage[chrom][strand][three_p_ss] = three_p_usage

	return fivep_ss_usage, threep_ss_usage

if __name__ == '__main__' :

	sj_info, out_dir = sys.argv[1:]

	cell_line_usage = sj_info.split(':')
	cell_line, sj_files = cell_line_usage[0], cell_line_usage[1].split(',')

	fivep_ss_usage, threep_ss_usage = ss_usage(sj_files)

	with open(join(out_dir, cell_line + '_5p_ss_usage.txt'), 'w') as out_file:
		out_file.write('Chrom\tSite\tStrand\tUsage\n')
		for chrom in fivep_ss_usage:
			for strand in fivep_ss_usage[chrom]:
				for fivep_site in fivep_ss_usage[chrom][strand]:
					fivep_usage = fivep_ss_usage[chrom][strand][fivep_site]
					out_file.write('\t'.join([chrom, str(fivep_site), strand, str(fivep_usage)]) + '\n')

	with open(join(out_dir, cell_line + '_3p_ss_usage.txt'), 'w') as out_file:
		out_file.write('Chrom\tSite\tStrand\tUsage\n')
		for chrom in threep_ss_usage:
			for strand in threep_ss_usage[chrom]:
				for threep_site in threep_ss_usage[chrom][strand]:
					threep_usage = threep_ss_usage[chrom][strand][threep_site]
					out_file.write('\t'.join([chrom, str(threep_site), strand, str(threep_usage)]) + '\n')
