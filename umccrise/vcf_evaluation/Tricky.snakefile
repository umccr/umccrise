from os.path import join

samples = [
	'batch1-ensemble-annotated',
	'batch1-ensemble-annotated-bwa',
	'batch1-mutect2-annotated',
	'batch1-mutect2-annotated-bwa',
	'batch1-strelka2-annotated',
	'batch1-strelka2-annotated-bwa',
	'batch1-vardict-annotated',
	'batch1-vardict-annotated-bwa',
]

tricky_bed = '/home/vlad/bcbio/genomes/Hsapiens/GRCh37/coverage/problem_regions/GA4GH/merged.sorted.bed.gz'


rule all:
	input:
		expand('{sample}_bcftools_isec/000{k}.anno.vcf', sample=samples, k=[0,1,2])

# rule fix_truth_vcf:
# 	input:
# 		'{sample}_bcftools_isec/000{k}.vcf'
# 	output:	
# 		'{sample}_bcftools_isec/000{k}.fixed.vcf'
# 	run:
# 		if input[0].endswith('0001.vcf'):
# 			shell('bcftools annotate {input} -h <(echo \'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\') > {output}')
# 		else:
# 			shell('cp {input} {output}')

# rule anno_tricky_bcftools:
# 	input: 
# 		vcf = rules.fix_truth_vcf.output[0]
# 	output:
# 		'{sample}_bcftools_isec/000{k}.fixed.anno.vcf.gz'
# 	run:
# 		cmd = 'bcftools view {input.vcf} -Ob'
# 		for tb in tricky_beds:
# 			key = tb.replace('.bed.gz', '')
# 			fp = join('/home/vlad/bcbio/genomes/Hsapiens/GRCh37/coverage/problem_regions/GA4GH', tb)
# 			cmd += (
# 				f' | bcftools annotate'
# 				f' -h <(echo \'##INFO=<ID=TRICKY_{key},Number=1,Type=String,Description="Tricky regions overlap">\')'
# 				f' -a {fp} -c CHROM,FROM,TO,TRICKY_{key} -Ob')
# 		cmd += ' | bcftools view -Oz -o {output} && tabix {output} -f'
# 		shell(cmd)

TOML = f'''[[annotation]]
file="{tricky_bed}"
names=["TRICKY"]
columns=[4]
ops=["self"]'''

rule anno_tricky:
	input:
		vcf = '{sample}_bcftools_isec/000{k}.vcf'
	output:
		'{sample}_bcftools_isec/000{k}.anno.vcf'
	params:
		toml = TOML.replace('\n', r'\n')
	shell:
		"vcfanno <(echo '{params.toml}') {input} > {output}"
