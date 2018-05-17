"""
Structural variants
------------------
Re-do the CNV plots. This will need lots of love (drop gene names, make the scatterplot viable again, etc.).
"""

#######################
######### CNV #########

#### Drop gene labels
rule cnvkit_cleanup:
    input:
        lambda wc: join(batch_by_name[wc.batch].tumor.dirpath, f'{batch_by_name[wc.batch].name}-cnvkit.cns')
    output:
        'work/{batch}/structural/{batch}-cnvkit-nolabels.cns'
    shell:
        'cat {input}'
        ' | grep -v ^GL '
        ' | py -x "\'\\t\'.join((x.split()[:3] + [\'.\'] + x.split()[4:]) if not x.startswith(\'chromosome\') else x.split())"'
        ' > {output}'

#### Plot
rule cnvkit_plot:
    input:
        rules.cnvkit_cleanup.output[0]
    output:
        '{batch}/structural/{batch}-cnvkit-diagram.pdf'
    shell:
        'cnvkit.py diagram -s {input} -o {output}'


#######################
######### SV ##########

rule prep_sv_vcf:
    input:
        manta_vcf = lambda wc: join(batch_by_name[wc.batch].tumor.dirpath, f'{batch_by_name[wc.batch].name}-sv-prioritize-manta.vcf.gz')
    output:
        vcf = '{batch}/structural/{batch}-sv-prioritize-manta-pass.vcf'
    shell:
        'bcftools view -f.,PASS,REJECT {input.manta_vcf} > {output}'

#### Bring in the prioritized SV calls from Manta. This should also include a basic plot at some stage.
rule prep_sv_tsv:
    input:
        sv_prio = lambda wc: join(batch_by_name[wc.batch].tumor.dirpath, f'{batch_by_name[wc.batch].name}-sv-prioritize.tsv'),
        vcf = rules.prep_sv_vcf.output.vcf
    output:
        '{batch}/structural/{batch}-sv-prioritize-manta-pass.tsv'
    shell:
        'head -n1 {input.sv_prio} > {output} && '
        'grep manta {input.sv_prio} | grep -f <(cut -f1,2 {input.vcf}) >> {output}'


#### At least for the most conservative manta calls generate a file for viewing in Ribbon ###

rule ribbon_filter_manta:
    input:
        manta_vcf = rules.prep_sv_vcf.output.vcf
    output:
        'work/{batch}/structural/ribbon/manta-pass.vcf'
    shell:
        'bcftools view {input.manta_vcf} > {output}'

rule ribbon_filter_vcfbedtope_starts:
    input:
        bed = rules.ribbon_filter_manta.output[0],
        fai = ref_fa + '.fai'
    output:
        'work/{batch}/structural/ribbon/manta-pass-strats.bed'
    shell:
        'cat {input.bed} | vcftobedpe'
        ' | cut -f 1-3'
        ' | bedtools slop -b 5000 -i stdin -g {input.fai}'
        ' > {output}'

rule ribbon_filter_vcfbedtope_ends:
    input:
        bed = rules.ribbon_filter_manta.output[0],
        fai = ref_fa + '.fai'
    output:
        'work/{batch}/structural/ribbon/manta-pass-ends.bed'
    shell:
        'cat {input.bed} | vcftobedpe'
        ' | cut -f 4-6'
        ' | grep -v \'CHROM\''
        ' | bedtools slop -b 5000 -i stdin -g {input.fai}'
        ' > {output}'

rule ribbon:
    input:
        starts = rules.ribbon_filter_vcfbedtope_starts.output[0],
        ends = rules.ribbon_filter_vcfbedtope_ends.output[0]
    output:
        '{batch}/structural/{batch}-sv-prioritize-manta-pass.ribbon.bed'
    shell:
        'cat {input.starts} {input.ends} | bedtools sort -i stdin | bedtools merge -i stdin > {output}'


#### Convert matna VCF to bedpe ####
rule bedpe:
    input:
        manta_vcf = rules.prep_sv_vcf.output.vcf
    output:
        '{batch}/structural/{batch}-sv-prioritize-manta-pass.bedpe'
    shell:
        'bcftools view {input.manta_vcf}'
        ' | vcftobedpe'
        ' | cut -f 1-7'
        ' > {output}'


#############

rule structural:
    input:
        expand(rules.bedpe.output, batch=batch_by_name.keys()),
        expand(rules.ribbon.output, batch=batch_by_name.keys()),
        expand(rules.prep_sv_tsv.output, batch=batch_by_name.keys()),
        expand(rules.cnvkit_plot.output, batch=batch_by_name.keys())
    output:
        temp(touch('structural.done'))
