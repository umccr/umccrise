"""
Structural variants
------------------
Re-do the CNV plots. This will need lots of love (drop gene names, make the scatterplot viable again, etc.).
"""
from cyvcf2 import VCF

vcftobedpe = 'vcfToBedpe'


localrules: prep_sv_vcf, filter_sv_vcf, prep_sv_tsv, ribbon_filter_manta, ribbon_filter_vcfbedtope_starts, ribbon_filter_vcfbedtope_ends, ribbon, structural


def get_manta_path(b):
    return join(batch_by_name[b].tumor.dirpath, f'{batch_by_name[b].name}-sv-prioritize-manta.vcf.gz')
def get_sv_tsv_path(b):
    return join(batch_by_name[b].tumor.dirpath, f'{batch_by_name[b].name}-sv-prioritize.tsv')
if not all(isfile(get_manta_path(b)) for b in batch_by_name.keys()):
    # CWL?
    def get_manta_path(b):
        return join(run.date_dir, batch_by_name[b].tumor.name + '-manta-prioritized.vcf.gz')
    def get_sv_tsv_path(b):
        return join(run.date_dir, batch_by_name[b].tumor.name + '-prioritize.tsv')

    if not all(isfile(get_manta_path(b)) for b in batch_by_name.keys()):
        critical('Could not find manta files for all batches neither under sample folders as '
                 '<tumor>/<batch>-sv-prioritize-manta.vcf.gz (conventional bcbio), nor in the project folder as'
                 'project/<tumor>-manta-prioritized.vcf.gz (CWL bcbio).')


rule prep_sv_vcf:
    """ Keep variants with the FILTER values *only* in PASS, Intergenic, or MissingAnn
        Note: `bcftools view -f .,PASS,Intergenic,MissingAnn` doesn't work because it keep variants 
        with other values in FILTER too. That's why we remove those FILTER values.
    """
    input:
        lambda wc: get_manta_path(wc.batch)
    output:
        vcf = '{batch}/structural/{batch}-sv-prioritize-manta.vcf'
    # group: "sv_vcf"
    shell:
        '(bcftools annotate -x "FILTER/Intergenic,FILTER/MissingAnn" {input} || bcftools annotate -x "FILTER/REJECT" {input}) | '
        'bcftools view -f .,PASS > {output}'

rule filter_sv_vcf:
    """ Keep all with read support above 10x; or allele frequency above 10%, but only if read support is above 5x
    """
    input:
        vcf = rules.prep_sv_vcf.output.vcf
    output:
        vcf = '{batch}/structural/{batch}-sv-prioritize-manta-filter.vcf'
    # group: "sv_vcf"
    run:
        print(f'VCF samples: {VCF(input.vcf).samples}')
        print(f'Bcbio batch tumor name:: {batch_by_name[wildcards.batch].tumor.name}')
        tumor_id = VCF(input.vcf).samples.index(batch_by_name[wildcards.batch].tumor.name)
        print(f'Derived tumor VCF index: {tumor_id}')
        shell('''cat {input.vcf} | \
bcftools filter -e "FORMAT/SR[{tumor_id}:1]<5  & FORMAT/PR[{tumor_id}:1]<5" | \
bcftools filter -e "FORMAT/SR[{tumor_id}:1]<10 & FORMAT/PR[{tumor_id}:1]<10 & (BPI_AF[0] < 0.1 | BPI_AF[1] < 0.1)" > {output.vcf}
''')

#### Bring in the prioritized SV calls from Manta. This should also include a basic plot at some stage.
rule prep_sv_tsv:
    input:
        sv_prio = lambda wc: get_sv_tsv_path(wc.batch),
        vcf = rules.filter_sv_vcf.output.vcf
    output:
        '{batch}/structural/{batch}-sv-prioritize-manta-pass.tsv'
    shell: """
head -n1 {input.sv_prio} > {output} && grep manta {input.sv_prio} | grep -f <(grep -v ^# {input.vcf} | cut -f1,2) | cat >> {output}
    """

#### At least for the most conservative manta calls generate a file for viewing in Ribbon ###

rule ribbon_filter_manta:
    input:
        manta_vcf = rules.filter_sv_vcf.output.vcf
    output:
        'work/{batch}/structural/ribbon/manta.vcf'
    # group: "ribbon"
    shell:
        'bcftools view {input.manta_vcf} > {output}'

rule ribbon_filter_vcfbedtope_starts:
    input:
        bed = rules.ribbon_filter_manta.output[0],
        fai = ref_fa + '.fai'
    output:
        'work/{batch}/structural/ribbon/manta-starts.bed'
    params:
        vcftobedpe = vcftobedpe
    # group: "ribbon"
    shell:
        'cat {input.bed} | {params.vcftobedpe}'
        ' | cut -f 1-3'
        ' | bedtools slop -b 5000 -i stdin -g {input.fai}'
        ' > {output}'

rule ribbon_filter_vcfbedtope_ends:
    input:
        bed = rules.ribbon_filter_manta.output[0],
        fai = ref_fa + '.fai'
    output:
        'work/{batch}/structural/ribbon/manta-ends.bed'
    params:
        vcftobedpe = vcftobedpe
    # group: "ribbon"
    shell:
        'cat {input.bed} | {params.vcftobedpe}'
        ' | cut -f 4-6'
        ' | grep -v \'CHROM\''
        ' | bedtools slop -b 5000 -i stdin -g {input.fai}'
        ' > {output}'

rule ribbon:
    input:
        starts = rules.ribbon_filter_vcfbedtope_starts.output[0],
        ends = rules.ribbon_filter_vcfbedtope_ends.output[0]
    output:
        '{batch}/structural/{batch}-sv-prioritize-manta.ribbon.bed'
    params:
        vcftobedpe = vcftobedpe
    # group: "ribbon"
    shell:
        'cat {input.starts} {input.ends} | bedtools sort -i stdin | bedtools merge -i stdin > {output}'


#### Convert matna VCF to bedpe ####
rule bedpe:
    input:
        manta_vcf = rules.filter_sv_vcf.output.vcf
    output:
        '{batch}/structural/{batch}-sv-prioritize-manta.bedpe'
    params:
        vcftobedpe = vcftobedpe
    shell:
        'bcftools view {input.manta_vcf}'
        ' | {params.vcftobedpe}'
        ' | cut -f 1-7'
        ' > {output}'


#############

rule structural:
    input:
        expand(rules.bedpe.output, batch=batch_by_name.keys()),
        expand(rules.ribbon.output, batch=batch_by_name.keys()),
        expand(rules.prep_sv_tsv.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/structural.done'))
