import sys

# Global path
g_dir = sys.argv[1]
ref_path = f'{g_dir}/reference/ref.fa'

# Get vcf data
with open('variants_biallelic.vcf') as vcf_f:
    vcf_l = vcf_f.readlines()

vcf_l = [(x.split()[1], x.split()[3], x.split()[4], x.split()[9].split(':')[0], x.split()[9].split(':')[2]) for x in vcf_l if x[0] != '#' and 'GT:PL:AD:GQ' in x]
positions = [int(x[0]) for x in vcf_l]
refs = [x[1] for x in vcf_l]
alts = [x[2] for x in vcf_l]
genotypes = [x[3] for x in vcf_l]

depths = [x[4] for x in vcf_l]
depths = [x.split(',') for x in depths]

depths_int = []
for depth in depths:
    depth_sum = 0
    for x in depth:
        depth_sum += int(x)

    depths_int.append(depth_sum)

variants = len(positions)

# Compare actual ref to refs in VCF file
with open(ref_path) as ref_f:
    ref_seq = ref_f.readlines()[1]

for pos, ref in zip(positions, refs):
    if ref_seq[pos - 1] != ref:
        print('REF - VCF MISMATCH AT {}'.format(pos))

# Find false positives, no segsite at actual chromosome
out_lines = []
false_positives = 0

with open('../../chr.1.fa') as chr1_f, open('../../chr.2.fa') as chr2_f:
    chr1_seq = chr1_f.readlines()[1]
    chr2_seq = chr2_f.readlines()[1]

for pos, ref, alt, gt in zip(positions, refs, alts, genotypes):
    chr1_site = chr1_seq[pos - 1]
    chr2_site = chr2_seq[pos - 1]
    ref_site = ref_seq[pos - 1]

    # Check for segsite
    if gt == '0/1':
        if chr1_site == chr2_site:
            false_positives += 1
            if chr1_site == ref:
                out_lines.append('FALSE POSITIVE AT {}, BOTH SITES ARE ACTUALLY REF \'{}\'\n'.format(pos, ref))
            if chr1_site == alt:
                out_lines.append('FALSE POSITIVE AT {}, BOTH SITES ARE ACTUALLY ALT \'{}\'\n'.format(pos, alt))
    elif gt == '1/1':
        if chr1_site != chr2_site:
            false_positives += 1
            out_lines.append('FALSE POSITIVE AT {}, SITES ARE DIFFERENT\n'.format(pos))


# Find false negatives, segsite not in vcf file
false_negatives = 0

for pos, s1, s2 in zip(range(len(chr1_seq)), chr1_seq, chr2_seq):
    if s1 != s2 and pos + 1 not in positions:
        false_negatives += 1
        out_lines.append('FALSE NEGATIVE AT {}, SITES ARE \'{}\' AND \'{}\'\n'.format(pos + 1, s1, s2))

# Output stats
if variants > 0:
    out_lines.append('FOUND {} FALSE POSITIVES OUT OF {} FILTERED VARIANTS: {}%\n'.format(
        false_positives, variants, false_positives / variants * 100.0))
    out_lines.append('FOUND {} FALSE NEGATIVES\n'.format(false_negatives))

out_lines.append('AVG DEPTH:\n')

avg_depth = 0
if len(depths_int) > 0:
    avg_depth = sum(depths_int) / float(len(depths_int))

out_lines.append('{}\n'.format(avg_depth))

out_lines.append('FILTERED VARIANTS:\n')
out_lines.append('{}'.format(variants))

with open('vcf_stats.txt', 'w') as out_f:
	for line in out_lines:
		out_f.write(line)
