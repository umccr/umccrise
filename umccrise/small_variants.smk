#################
#### Somatic ####
from os.path import isfile, join
from ngs_utils.file_utils import get_ungz_gz
from umccrise import package_path
import toml
from vcf_stuff import iter_vcf
import cyvcf2
import yaml
from ngs_utils.reference_data import get_key_genes_set
import subprocess


localrules: small_variants, prep_anno_toml, prep_giab_bed, somatic_stats_report, germline_stats_report


SUBSET_SOMATIC_IF_MORE_THAN = 500_000


# Call variants with 1%
# Apply to VarDict only: (INFO/QUAL * TUMOR_AF) >= 4
# Call ensemble
# First PoN round: remove PoN_CNT>2'
# Annotate with PCGR (VEP+known cancer databases)
#   Tier 1 - variants of strong clinical significance
#   Tier 2 - variants of potential clinical significance
#   Tier 3 - variants of unknown clinical significance
#   Tier 4 - other coding variants
#   Noncoding variants
# Tier 1-3 - keep all variants
# Tier 4 and noncoding - filter with:
#   Remove gnomad_AF <=0.02
#   Remove PoN_CNT>{0 if issnp else 1}'
#   Remove indels in "bad_promoter" tricky regions
#   Remove DP<25 & AF<5% in tricky regions: gc15, gc70to75, gc75to80, gc80to85, gc85, low_complexity_51to200bp,
#           low_complexity_gt200bp, non-GIAB confident, unless coding in cancer genes

vcf_suffix = '-annotated'
def get_somatic_vcf_path(b, suf):
    return join(run.date_dir, f'{batch_by_name[b].name}-{run.somatic_caller}{suf}.vcf.gz')
def get_germline_vcf_path(b, suf):
    return join(run.date_dir, f'{batch_by_name[b].normal.name}-germline-{run.germline_caller}{suf}.vcf.gz')
if not all(isfile(get_somatic_vcf_path(b, vcf_suffix)) for b in batch_by_name.keys()):
    # CWL?
    vcf_suffix = ''
    if not all(isfile(get_somatic_vcf_path(b, vcf_suffix)) for b in batch_by_name.keys()):
        critical('Could not find somatic variants files for all batches neither as '
                 '<datestamp_dir>/<batch>-{run.somatic_caller}-annotated.vcf.gz (conventional bcbio), '
                 'nor as project/<batch>-{run.somatic_caller}.vcf.gz (CWL bcbio).')


# if total_vars > 500_000:  # to avoid PCGR choking on too many variants
rule somatic_vcf_keygenes:
    input:
        vcf = lambda wc: get_somatic_vcf_path(wc.batch, vcf_suffix)
    output:
        'work/{batch}/small_variants/keygenes_subset/{batch}-somatic-' + run.somatic_caller + '.vcf.gz',
    run:
        genes = get_key_genes_set()
        def func(rec):
            if rec.INFO.get('ANN') is not None and rec.INFO['ANN'].split('|')[3] in genes:
                return rec
        iter_vcf(input.vcf, output.vcf, func)

def cnt_vars(vcf_path, passed=False):
    snps = 0
    indels = 0
    others = 0
    for rec in cyvcf2.VCF(vcf_path):
        if passed and rec.FILTER is not None and rec.FILTER != 'PASS':
            continue
        if rec.is_snp:
            snps += 1
        elif rec.is_indel:
            indels += 1
        else:
            others += 1
    return snps, indels, others

# def count_vars(vcf_fpath):
#     return int(subprocess.check_output(f'bcftools view -H {vcf_fpath} | wc -l', shell=True).strip())

rule somatic_vcf_sort:
    input:
        subset_vcf = rules.somatic_vcf_keygenes.output[0],
        full_vcf = lambda wc: get_somatic_vcf_path(wc.batch, vcf_suffix),
    output:
        'work/{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-SORT.vcf.gz',
    run:
        vcf = input.subset_vcf \
            if sum(cnt_vars(input.full_vcf)) > SUBSET_SOMATIC_IF_MORE_THAN \
            else input.full_vcf
        shell(f'(bcftools view -h {vcf} ; bcftools view -H -f.,PASS {vcf} | sort -k1,1V -k2,2n) | bgzip -c > {output} && '
               'tabix -f -p vcf {output}')

rule somatic_vcf_annotate:
    input:
        vcf = rules.somatic_vcf_sort.output[0],
    output:
        vcf = 'work/{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-ANNO.vcf.gz',
    params:
        genome = run.genome_build,
        work_dir = 'work/{batch}/small_variants',
#        full_out_path = lambda wc, input, output: abspath(output.vcf)
    # group: "small_variants"
    resources:
        mem_mb = 20000
    shell:
        'anno_somatic_vcf {input.vcf} -g {params.genome} -o {output.vcf} -w {params.work_dir}'

rule somatic_extract_tumor_sample:
    input:
        vcf = rules.somatic_vcf_annotate.output.vcf,
    output:
        vcf = 'work/{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-ANNO-onesample.vcf.gz',
    params:
        tumor_sample = lambda wc: batch_by_name[wc.batch].tumor.name,
    # group: "small_variants"
    shell:
        'bcftools view -s {params.tumor_sample} {input.vcf} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'

rule somatic_vcf_filter:
    input:
        vcf = rules.somatic_extract_tumor_sample.output.vcf,
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '.vcf.gz',
    # group: "small_variants"
    shell:
        'filter_somatic_vcf {input.vcf} -o {output.vcf}'

rule somatic_vcf_filter_pass:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS.vcf.gz.tbi',
    # group: "small_variants"
    shell:
        'bcftools view -f.,PASS {input.vcf} -Oz -o {output.vcf} && tabix -f -p vcf {output.vcf}'

rule somatic_vcf_filter_pass_af10:
    input:
        vcf = rules.somatic_vcf_filter_pass.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS-AF10.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS-AF10.vcf.gz.tbi',
    # group: "small_variants"
    shell:
        'bcftools filter -e "INFO/TUMOR_AF<0.1" -Oz {input.vcf} -o {output.vcf} && tabix -f -p vcf -f {output.vcf}'

rule bcftools_stats_somatic:
    input:
        rules.somatic_vcf_filter.output.vcf
    output:
        '{batch}/small_variants/stats/{batch}_bcftools_stats.txt'
    params:
        sname = lambda wc: batch_by_name[wc.batch].tumor.name,
    shell:
        'bcftools stats -s {params.sname} {input} | sed s#{input}#{params.sname}# > {output}'

##################
#### Germline ####
# Annotate any events found in ~200 cancer predisposition gene set.
rule germline_vcf_subset:  # {batch}
    input:
        vcf = lambda wc: get_germline_vcf_path(wc.batch, vcf_suffix)
    output:
        vcf = 'work/{batch}/small_variants/raw_normal-' + run.germline_caller + '-predispose_genes.vcf.gz',
        tbi = 'work/{batch}/small_variants/raw_normal-' + run.germline_caller + '-predispose_genes.vcf.gz.tbi',
    params:
        ungz = lambda wc, output: get_ungz_gz(output[0])[0]
    # group: "small_variants"
    run:
        pcgr_toml_fpath = join(package_path(), 'pcgr', 'cpsr.toml')
        genes = [g for g in toml.load(pcgr_toml_fpath)['cancer_predisposition_genes']]
        def func(rec):
            if rec.INFO.get('ANN') is not None and rec.INFO['ANN'].split('|')[3] in genes:
                return rec
        iter_vcf(input.vcf, output.vcf, func)

# Preparations: annotate TUMOR_X and NORMAL_X fields, remove non-standard chromosomes and mitochondria, remove non-PASSed calls.
# Suites for PCGR, but for all other processing steps too
rule germline_vcf_prep:
    input:
        vcf = rules.germline_vcf_subset.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-normal-' + run.germline_caller + '-predispose_genes.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-normal-' + run.germline_caller + '-predispose_genes.vcf.gz.tbi',
    # group: "small_variants"
    shell:
        'pcgr_prep {input.vcf} |'
        ' bcftools view -f.,PASS -Oz -o {output.vcf}'
        ' && tabix -f -p vcf {output.vcf}'

# rule bcftools_stats_germline:
#     input:
#         vcf = lambda wc: get_germline_vcf_path(wc.batch, vcf_suffix)
#     output:
#         '{batch}/small_variants/stats/{batch}_bcftools_stats_germline.txt'
#     params:
#         sname = lambda wc: batch_by_name[wc.batch].normal.name,
#     shell:
#         'bcftools stats -s {params.sname} {input} | sed s#{input}#{params.sname}# > {output}'

rule somatic_stats_report:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf,
        full_vcf = lambda wc: get_somatic_vcf_path(wc.batch, vcf_suffix)
    output:
        'work/{batch}/small_variants/somatic_stats.yml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name,
    run:
        snps, indels, others = cnt_vars(input.vcf, passed=True)
        total_snps, total_indels, total_others = cnt_vars(input.vcf)

        # pon_cnt = 0
        # gnomad_cnt = 0
        # pon_total_cnt = 0
            # if 'PoN' in (rec.FILTER or ''):
            #     pon_filt += 1
            # if rec.INFO.get('PoN', 0) > 0:
            #     pon_cnt += 1
        # pass_pct = pass_cnt / total_cnt
        # pon_pct = pon_cnt / total_cnt
        # pon_filt_pct = pon_filt / total_cnt

        ori_snps, ori_indels, ori_others = cnt_vars(input.full_vcf)
        data = dict(
            snps=snps,
            indels=indels,
            others=others if others else None,
            filt_snps=total_snps - snps,
            filt_indels=total_indels - indels,
            filt_others=((total_others - others) if others and total_others else None),
        )
        if ori_snps + ori_indels + ori_others > SUBSET_SOMATIC_IF_MORE_THAN:
            data.update(dict(
                ori_snps=ori_snps,
                ori_indels=ori_indels,
                ori_others=ori_others,
            ))

        with open(output[0], 'w') as out:
            data = {
                'data': {
                    params.sample: data
                }
            }
            yaml.dump(data, out, default_flow_style=False)

            # out.write('\t'.join(['Sample Name', 'somatic', 'filtered', 'snps', 'indels'] + [['other'] if others else [])) + '\n')
            # out.write('\t'.join(map(str, [params.sample, pass_cnt, total_cnt - pass_cnt, ])) + '\n')
            # yaml.dump(data, out, default_flow_style=False)

        # data = dict(
        #     id='umccrise',
        #     plot_type='generalstats',
        #     # pconfig=[
        #     #     dict(
        #     #         pon=dict(
        #     #             title='PoN',
        #     #             description='Variants found in UMCCR panel of normals',
        #     #             id='pon',
        #     #             min=0,
        #     #             max=100,
        #     #             scale='RdYlGn',
        #     #             suffix='%',
        #     #         ),
        #     #     ),
        #     #     dict(
        #     #         passed=dict(
        #     #             title='PASS',
        #     #             description='Variants passed umccrise somatic filters',
        #     #             id='passed',
        #     #             min=0,
        #     #             max=100,
        #     #             scale='RdYlGn',
        #     #             suffix='%',
        #     #         ),
        #     #     ),
        #     #     dict(
        #     #         pon_filt=dict(
        #     #             title='PoN filt',
        #     #             description='Variants filtered due to hits in UMCCR panel of normals '
        #     #                         '(filtering threshold - 2 hits, but keeping all hotspots and TIER 1-3)',
        #     #             id='pon_filt',
        #     #             min=0,
        #     #             max=100,
        #     #             scale='RdYlGn',
        #     #             suffix='%',
        #     #             hidden=True,
        #     #         )
        #     #     )
        #     # ],
        #     data={
        #         params.sample: {
        #             'PoN': pon_pct,
        #             'PASS': pass_cnt,
        #             'PoN_filt': pon_filt_pct,
        #             # TODO: add germline hits (gnomad cnt?)
        #         }
        #     },
        # )


rule germline_stats_report:
    input:
        vcf = lambda wc: get_germline_vcf_path(wc.batch, vcf_suffix)
    output:
        'work/{batch}/small_variants/germline_stats.yml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].normal.name
    run:
        pass_cnt = 0
        vcf = cyvcf2.VCF(input.vcf)
        for rec in vcf:
            if rec.FILTER is None or rec.FILTER == 'PASS':
                pass_cnt += 1

        with open(output[0], 'w') as out:
            data = {
                'data': {
                    params.sample: dict(
                        germline = pass_cnt
                    )
                }
            }
            yaml.dump(data, out, default_flow_style=False)


#############

rule small_variants:
    input:
        expand(rules.somatic_vcf_filter.output.vcf, batch=batch_by_name.keys()),
        expand(rules.somatic_vcf_filter_pass.output.vcf, batch=batch_by_name.keys()),
        expand(rules.somatic_vcf_filter_pass_af10.output.vcf, batch=batch_by_name.keys()),
        expand(rules.germline_vcf_prep.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/small_variants.done'))




