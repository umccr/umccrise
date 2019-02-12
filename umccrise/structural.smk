"""
Structural variants
------------------
Re-do the CNV plots. This will need lots of love (drop gene names, make the scatterplot viable again, etc.).
"""
import itertools

from cyvcf2 import VCF

from ngs_utils.file_utils import safe_mkdir
from ngs_utils.reference_data import get_key_genes_txt, get_known_fusion_pairs, get_known_fusion_heads, get_known_fusion_tails
from vcf_stuff import count_vars, vcf_contains_field, iter_vcf

vcftobedpe = 'vcfToBedpe'


localrules: structural


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

rule prep_sv_prio_lists:
    input:
        fusion_pairs = get_known_fusion_pairs(),
        fusion_heads = get_known_fusion_heads(),
        fusion_tails = get_known_fusion_tails(),
    output:
        pairs = 'work/fusion_pairs.txt',
        promisc = 'work/fusion_promisc.txt',
    run:
        with open(input.fusion_pairs) as f, open(output.pairs, 'w') as out:
            for l in f:
                l = l.strip().replace('"', '')
                if l:
                    g1, g2 = l.split(',')[0:2]
                    if g1 and g2 and g1 != 'H_gene':
                        out.write(f'{g1},{g2}\n')
        with open(input.fusion_heads) as f1, open(input.fusion_tails) as f2, open(output.promisc, 'w') as out:
            for l in itertools.chain(f1, f2):
                l = l.strip().replace('"', '')
                if l:
                    gene = l.split(',')[0]
                    if gene and gene != 'gene':
                        out.write(f'{gene}\n')

rule sv_prioritize:
    input:
        vcf = lambda wc: get_manta_path(wc.batch),
        known_pairs = rules.prep_sv_prio_lists.output.pairs,
        known_promisc = rules.prep_sv_prio_lists.output.promisc,
        cancer_genes = get_key_genes_txt(),
    output:
        vcf = 'work/{batch}/structural/prioritize/{batch}-manta.vcf.gz'
    run:
        cmd = f'cat {input.vcf}'
        # remove previous annotation
        filts_to_remove = [f'{f}' for f in ['INFO/SIMPLE_ANN', 'INFO/SV_HIGHEST_TIER',
                                            'FILTER/Intergenic', 'FILTER/MissingAnn', 'FILTER/REJECT']
                           if vcf_contains_field(input.vcf, f)]
        if filts_to_remove:
            cmd += f' | bcftools annotate -x "' + ','.join(f'{f}' for f in filts_to_remove) + '"'

        cmd += (f' | simple_sv_annotation.py - '
                f'--known_fusion_pairs {input.known_pairs} --known_fusion_promiscuous {input.known_promisc} '
                f'--gene_list {input.cancer_genes} -o - | bgzip -c > {output.vcf} && tabix -p vcf {output.vcf}')
        shell(cmd)
        before = count_vars(input.vcf)
        after = count_vars(output.vcf)
        assert before == after, (before, after)

# Keep variants with the FILTER values only in PASS, Intergenic, or MissingAnn.
# Unlesss there are too many SV calls (FFPE?) - in this case keep *only* PASS.
rule sv_keep_pass:
    input:
        vcf = rules.sv_prioritize.output.vcf
    output:
        vcf = 'work/{batch}/structural/keep_pass/{batch}-manta.vcf.gz'
    group: "sv_vcf"
    run:
        cmd = f'cat {input.vcf}'
        if count_vars(input.vcf, filter='.,PASS,Intergenic,MissingAnn') <= 1000:
            # Not too many variants, so Intergenic won't clutter the reports much - keeping them
            filts_to_remove = [f'FILTER/{f}' for f in ['Intergenic', 'MissingAnn']
                               if vcf_contains_field(input[0], f, 'FILTER')]
            if filts_to_remove:
                # Want to keep Intergenic variants from sv-prioritize.
                # Note that `bcftools view -f .,PASS,Intergenic,MissingAnn` doesn't work because
                # it keep variants with other values in FILTER also. That's why we remove those
                # FILTER values with `bcftools annotate -x`.
                cmd += f' | bcftools annotate -x "' + ','.join(f'{f}' for f in filts_to_remove) + '"'
        cmd += f' | bcftools view -f.,PASS'
        if vcf_contains_field(input[0], 'ANN', 'INFO'):  # very cluttered ANN field, ran for all tanscripts
            cmd += ' | bcftools annotate -x INFO/ANN'
        cmd += f' -Oz -o {output.vcf}'
        shell(cmd)

rule sv_maybe_keep_prioritize:
    input:
        vcf = rules.sv_keep_pass.output.vcf
    output:
        vcf = 'work/{batch}/structural/maybe_keep_prio/{batch}-manta.vcf.gz'
    group: "sv_vcf"
    run:
        if count_vars(input.vcf, filter='.,PASS') > 1000:
            # Still too cluttered (FFPE?) - removing all NOT_PRIORITISED as well
            def func(rec):
                ann = rec.INFO.get('SIMPLE_ANN')
                if ann:
                    anns = [a for a in ann.split(',') if '|NOT_PRIORITISED|' not in a]
                    if anns:
                        ann = ','.join(anns)
                        rec.INFO['SIMPLE_ANN'] = ann
                        return rec
            iter_vcf(input.vcf, output.vcf, func)
        else:
            shell('cp {input.vcf} {output.vcf}')

# if BPI was disabled in bcbio
rule sv_maybe_bpi:
    input:
        vcf = rules.sv_maybe_keep_prioritize.output.vcf,
        tumor_bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
    output:
        vcf = 'work/{batch}/structural/maybe_bpi/{batch}-manta.vcf'
    group: "sv_vcf"
    log:
        'log/structural/{batch}/{batch}-bpi_stats.txt'
    params:
        xms = 1000,
        xmx = 2800,
        tmp_dir = '{batch}/structural/maybe_bpi/tmp_dir'
    resources:
        mem_mb = 3000
    run:
        if not vcf_contains_field(input.vcf, 'BPI_AF', 'INFO'):
            safe_mkdir(params.tmp_dir)
            shell(
                'break-point-inspector -Xms{params.xms}m -Xmx{params.xmx}m '
                '-Djava.io.tmpdir={params.tmp_dir} '
                '-vcf {input.vcf} '
                '-ref {input.normal_bam} '
                '-tumor {input.tumor_bam} '
                '-output_vcf {output.vcf} '
                '> {log}'
            )
        else:
            shell('cp {input.vcf} {output.vcf}')

# Keep all with read support above 10x; or allele frequency above 10%, but only if read support is above 5x
rule filter_sv_vcf:
    input:
        vcf = rules.sv_maybe_bpi.output.vcf
    output:
        vcf = '{batch}/structural/{batch}-manta.vcf'
    group: "sv_vcf"
    run:
        print(f'VCF samples: {VCF(input.vcf).samples}')
        print(f'Bcbio batch tumor name:: {batch_by_name[wildcards.batch].tumor.name}')
        tumor_id = VCF(input.vcf).samples.index(batch_by_name[wildcards.batch].tumor.name)
        print(f'Derived tumor VCF index: {tumor_id}')
        shell('''
bcftools view -f.,PASS {input.vcf} | \
bcftools filter -e "FORMAT/SR[{tumor_id}:1]<5  & FORMAT/PR[{tumor_id}:1]<5" | \
bcftools filter -e "FORMAT/SR[{tumor_id}:1]<10 & FORMAT/PR[{tumor_id}:1]<10 & (BPI_AF[0] < 0.1 | BPI_AF[1] < 0.1)" \
> {output.vcf}
''')

# Bring in the prioritized SV calls from Manta
rule prep_sv_tsv:
    input:
        sv_prio = lambda wc: get_sv_tsv_path(wc.batch),
        vcf = rules.filter_sv_vcf.output.vcf
    output:
        '{batch}/structural/{batch}-manta.tsv'
    group: "sv_vcf"
    shell: """
head -n1 {input.sv_prio} > {output} && grep manta {input.sv_prio} | grep -f <(grep -v ^# {input.vcf} | cut -f1,2) | cat >> {output}
    """

# At least for the most conservative manta calls, generate a file for viewing in Ribbon
rule ribbon_filter_manta:
    input:
        manta_vcf = rules.filter_sv_vcf.output.vcf
    output:
        'work/{batch}/structural/ribbon/manta.vcf'
    group: "sv_vcf"
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
    group: "sv_vcf"
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
    group: "sv_vcf"
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
        '{batch}/structural/{batch}-manta.ribbon.bed'
    params:
        vcftobedpe = vcftobedpe
    group: "sv_vcf"
    shell:
        'cat {input.starts} {input.ends} | bedtools sort -i stdin | bedtools merge -i stdin > {output}'


#### Convert matna VCF to bedpe ####
rule bedpe:
    input:
        manta_vcf = rules.filter_sv_vcf.output.vcf
    output:
        '{batch}/structural/{batch}-manta.bedpe'
    params:
        vcftobedpe = vcftobedpe
    group: "sv_vcf"
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
