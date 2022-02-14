import sys

# Global path
g_dir = sys.argv[1]
ref_path = f'{g_dir}/reference/ref.fa'

# Read sequence data
with open('../../chr.1.fa', 'r') as chr1_f, open('../../chr.2.fa', 'r') as chr2_f, open(ref_path, 'r') as ref_f:
    chr1_seq = chr1_f.readlines()[1]
    chr2_seq = chr2_f.readlines()[1]
    ref_seq = ref_f.readlines()[1]

out_lines = []
out_lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\treadset_sort_nodup.bam\n')

for pos in range(len(chr1_seq)):
    ref_allele = ref_seq[pos]
    chr1_allele = chr1_seq[pos]
    chr2_allele = chr2_seq[pos]

    # Check for variant
    if chr1_allele == chr2_allele == ref_allele:
        continue

    genotype = '0/1'
    alt_allele = ''
    if chr1_allele == chr2_allele:
        genotype = '1/1'
        alt_allele = chr1_allele
    elif chr1_allele != ref_allele:
        alt_allele = chr1_allele
    elif chr2_allele != ref_allele:
        alt_allele = chr2_allele
        

    out_lines.append('ref_1\t{}\t.\t{}\t{}\t200\t.\t.\tGT:AD\t{}:{},{}\n'.format(pos + 1, ref_allele, alt_allele, genotype, 5, 5))

with open('variants_biallelic.vcf', 'w') as out:
    for line in out_lines:
        out.write(line)
