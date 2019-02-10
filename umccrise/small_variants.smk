#################
#### Somatic ####
from os.path import isfile, join
from ngs_utils.file_utils import get_ungz_gz
from umccrise import package_path
import toml
from vcf_stuff import iter_vcf
import cyvcf2
import yaml


localrules: small_variants


vcf_suffix = '-annotated'
def get_somatic_vcf_path(b, suf):
    return join(run.date_dir, f'{batch_by_name[b].name}-{run.somatic_caller}{suf}.vcf.gz')
def get_germline_vcf_path(b, suf):
    return join(run.date_dir, f'{batch_by_name[b].normal.name}-germline-{run.germline_caller}{suf}.vcf.gz')
if not all(isfile(get_somatic_vcf_path(b, vcf_suffix)) for b in batch_by_name.keys()):
    # CWL?
    vcf_suffix = ''
    if not all(isfile(get_somatic_vcf_path(b, vcf_suffix)) for b in batch_by_name.keys()):
        critical(f'Could not find somatic variants files for all batches neither as '
                 f'{run.date_dir}/<batch>-{run.somatic_caller}-annotated.vcf.gz (conventional bcbio), '
                 f'nor as project/<batch>-{run.somatic_caller}.vcf.gz (CWL bcbio).')

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

rule somatic_vcf_pass_sort:
    input:
        vcf = lambda wc: get_somatic_vcf_path(wc.batch, vcf_suffix),
    output:
        vcf = 'work/{batch}/small_variants/input/{batch}-somatic-' + run.somatic_caller + '.vcf.gz',
    group: "somatic_anno"
    shell:
        '(bcftools view -h {input.vcf} ; bcftools view -H -f.,PASS {input.vcf} | sort -k1,1V -k2,2n) | '
        'bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'

rule somatic_vcf_annotate:
    input:
        vcf = rules.somatic_vcf_pass_sort.output[0],
    output:
        vcf = 'work/{batch}/small_variants/annotate/{batch}-somatic-' + run.somatic_caller + '.vcf.gz',
        subset_to_cancer = 'work/{batch}/small_variants/somatic_anno/subset_to_cancer_genes.flag',
    params:
        genome = run.genome_build,
        genomes_dir_opt = f" --genomes-dir {config.get('genomes_dir')}" if config.get('genomes_dir') else "",
        work_dir = 'work/{batch}/small_variants',
    resources:
        mem_mb = 20000
    group: "somatic_anno"
    shell:
        'anno_somatic_vcf {input.vcf} -g {params.genome} -o {output.vcf} -w {params.work_dir}{params.genomes_dir_opt}'

rule run_sage:
    input:
        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
        coding_bed = get_ref_file(run.genome_build, 'coding_regions'),
        ref_fa = ref_fa,
        hotspots = get_ref_file(run.genome_build, key='hmf_hotspot'),
    output:
        vcf = 'work/{batch}/small_variants/sage/{batch}-somatic-' + run.somatic_caller + '.vcf.gz'
    params:
        jar = join(package_path(), 'jars', 'sage-1.0-jar-with-dependencies.jar'),
        rundir = 'work/{batch}/purple',
        outdir = 'work/{batch}/purple',
        normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.name,
        xms = 2000,
        xmx = min(40000, 3000*threads_per_batch),
    shell:
        'java -Xms{params.xms}m -Xmx{params.xmx}m -cp {params.jar} com.hartwig.hmftools.sage.SageHotspotApplication '
        '-tumor {params.tumor_sname} -tumor_bam {input.tumor_bam} '
        '-reference {params.normal_sname} -reference_bam {input.normal_bam} '
        '-known_hotspots <(gunzip -c {input.hotspots}) '
        '-coding_regions {input.coding_bed} '
        '-ref_genome {input.ref_fa} '
        '-out {output.vcf} '
        '&& tabix -p vcf {output.vcf}'

rule annotate_from_sage:
    input:
        vcf = rules.somatic_vcf_annotate.output.vcf,
        sage_vcf = rules.run_sage.output.vcf,
    output:
        vcf = 'work/{batch}/small_variants/annotate/{batch}-somatic-' + run.somatic_caller + '-sage.vcf.gz'
    group: "somatic_filt"
    run:
        with open(input.sage_vcf + '_toml', 'w') as f:
            f.write(f"""[[annotation]]
file="{input.sage_vcf}"
fields = ["FILTER", "AF", "HOTSPOT"]
names = ["SAGE_FILTER", "SAGE_AF", "SAGE_HOTSPOT"]
ops = ["self", "self", "self"]
""")
        shell(f'vcfanno {input.sage_vcf}_toml {input.vcf} | bgzip -c > {output.vcf} && tabix -p vcf {output.vcf}')

rule somatic_vcf_filter:
    input:
        vcf = rules.annotate_from_sage.output.vcf,
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '.vcf.gz',
        # vcf = 'work/{batch}/small_variants/filter/{batch}-somatic-' + run.somatic_caller + '.vcf.gz',
    group: "somatic_filt"
    shell:
        'filter_somatic_vcf {input.vcf} -o {output.vcf}'

# rule somatic_extract_tumor_sample:
#     input:
#         vcf = rules.somatic_vcf_filter.output.vcf,
#     output:
#         vcf = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '.vcf.gz',
#     params:
#         tumor_sample = lambda wc: batch_by_name[wc.batch].tumor.name,
#     # group: "small_variants"
#     shell:
#         'bcftools view -s {params.tumor_sample} {input.vcf} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'

rule somatic_vcf_filter_pass:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS.vcf.gz.tbi',
    group: "somatic_filt"
    shell:
        'bcftools view -f.,PASS {input.vcf} -Oz -o {output.vcf} && tabix -f -p vcf {output.vcf}'

rule somatic_vcf_filter_pass_af10:
    input:
        vcf = rules.somatic_vcf_filter_pass.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS-AF10.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-somatic-' + run.somatic_caller + '-PASS-AF10.vcf.gz.tbi',
    group: "somatic_filt"
    shell:
        'bcftools filter -e "INFO/TUMOR_AF<0.1" -Oz {input.vcf} -o {output.vcf} && tabix -f -p vcf -f {output.vcf}'

rule bcftools_stats_somatic:
    input:
        rules.somatic_vcf_filter_pass.output.vcf
    output:
        '{batch}/small_variants/stats/{batch}_bcftools_stats.txt'
    group: "somatic_filt"
    params:
        sname = lambda wc: batch_by_name[wc.batch].tumor.name,
    shell:
        'bcftools stats -s {params.sname} {input} | sed s#{input}#{params.sname}# > {output}'

##################
#### Germline ####
rule germline_vcf_pass:
    input:
        vcf = lambda wc: get_germline_vcf_path(wc.batch, vcf_suffix),
    output:
        vcf = 'work/{batch}/small_variants/germline/{batch}-normal-' + run.germline_caller + '-PASS.vcf.gz',
    group: "germline_snv"
    shell:
        'bcftools view {input.vcf} -f.,PASS -Oz -o {output.vcf} && tabix -f -p vcf {output.vcf}'

# Annotate any events found in ~200 cancer predisposition gene set.
rule germline_vcf_subset:
    input:
        vcf = rules.germline_vcf_pass.output.vcf,
    output:
        vcf = 'work/{batch}/small_variants/germline/{batch}-normal-' + run.germline_caller + '-predispose_genes.vcf.gz',
        tbi = 'work/{batch}/small_variants/germline/{batch}-normal-' + run.germline_caller + '-predispose_genes.vcf.gz.tbi',
    params:
        ungz = lambda wc, output: get_ungz_gz(output[0])[0]
    group: "germline_snv"
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
    group: "germline_snv"
    shell:
        'pcgr_prep {input.vcf} | bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'

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
        full_vcf = rules.somatic_vcf_pass_sort.output.vcf,
        subset_to_cancer = rules.somatic_vcf_annotate.output.subset_to_cancer,
    output:
        'work/{batch}/small_variants/somatic_stats.yml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name,
    group: "somatic_filt"
    run:
        snps, indels, others = cnt_vars(input.vcf, passed=True)
        all_snps, all_indels, all_others = cnt_vars(input.full_vcf)
        ori = None
        if open(input.subset_to_cancer).read().strip() == 'YES':
            # In this case, "all" correspond to pre-cancer-gene subset variants,
            # Thus we can't compare it with cancer-gene-subset PASSing variants
            ori = all_snps + all_indels + all_others
            all_snps, all_indels, all_others = cnt_vars(input.vcf, passed=False)

        data = dict(
            snps=snps,
            indels=indels,
            others=others if others else None,
            filt_snps=(all_snps - snps) / all_snps * 100.0,
            filt_indels=(all_indels - indels) / all_indels * 100.0,
            filt_others=(((all_others - others) / all_others * 100.0) if others and all_others else None),
        )
        if ori is not None:
            data.update(dict(
                ori = ori
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
        vcf = rules.germline_vcf_pass.output.vcf,
    output:
        'work/{batch}/small_variants/germline_stats.yml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].normal.name
    group: "germline_snv"
    run:
        pass_cnt = 0
        vcf = cyvcf2.VCF(input.vcf)
        for rec in vcf:
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




