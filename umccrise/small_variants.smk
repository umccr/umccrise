#################
#### Somatic ####
import os
import subprocess
from os.path import isfile, join, dirname
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.reference_data import get_predispose_genes_bed
from ngs_utils.logger import critical
from umccrise import package_path
import toml
from ngs_utils.vcf_utils import iter_vcf
import cyvcf2
import yaml
from hpc_utils import hpc


localrules: somatic, germline, small_variants


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

# rule somatic_vcf_reheader  # change RGIDs to tumor and normal names?

rule somatic_vcf_pass_sort:
    input:
        vcf = lambda wc: batch_by_name[wc.batch].somatic_vcf,
    output:
        vcf = 'work/{batch}/small_variants/pass_sort/{batch}-somatic.vcf.gz',
        tbi = 'work/{batch}/small_variants/pass_sort/{batch}-somatic.vcf.gz.tbi',
    group: "somatic_anno"
    shell:
        '(bcftools view -h {input.vcf} ; bcftools view -H -f.,PASS {input.vcf} | sort -k1,1V -k2,2n) | '
        'bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'

rule somatic_vcf_select_noalt:
    input:
        vcf = rules.somatic_vcf_pass_sort.output.vcf,
        noalts_bed = hpc.get_ref_file(run.genome_build, 'noalt_bed'),
    output:
        vcf = 'work/{batch}/small_variants/noalt/{batch}-somatic.vcf.gz',
        tbi = 'work/{batch}/small_variants/noalt/{batch}-somatic.vcf.gz.tbi',
    group: "somatic_anno"
    shell:
        'bcftools view -R {input.noalts_bed} {input.vcf} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'

rule sage:
    input:
        vcf = rules.somatic_vcf_select_noalt.output.vcf,
        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
    output:
        vcf = 'work/{batch}/small_variants/sage/{batch}-somatic.vcf.gz',
        tbi = 'work/{batch}/small_variants/sage/{batch}-somatic.vcf.gz.tbi',
        sage_vcf = '{batch}/small_variants/sage/{batch}-sage.vcf.gz',
        sage_tbi = '{batch}/small_variants/sage/{batch}-sage.vcf.gz.tbi',
    params:
        genome = run.genome_build,
        genomes_dir = hpc.genomes_dir,
        work_dir = 'work/{batch}/small_variants/sage',
        unlock_opt = ' --unlock' if config.get('unlock', 'no') == 'yes' else '',
        tumor_sample = lambda wc: batch_by_name[wc.batch].tumor.rgid,
        normal_sample = lambda wc: batch_by_name[wc.batch].normal.rgid,
    resources:
        mem_mb = 20000
    group: "somatic_anno"
    shell:
        'sage -t {input.tumor_bam} -n {input.normal_bam} -v {input.vcf} -o {output.vcf} -s {output.sage_vcf} '
        '-tn {params.tumor_sample} -nn {params.normal_sample} '
        '-w {params.work_dir} -g {params.genome} --genomes-dir {params.genomes_dir} '
        '{params.unlock_opt}'

rule somatic_vcf_annotate:
    input:
        vcf = rules.sage.output.vcf,
    output:
        vcf = 'work/{batch}/small_variants/annotate/{batch}-somatic.vcf.gz',
        subset_highly_mutated_stats = 'work/{batch}/small_variants/somatic_anno/subset_highly_mutated_stats.yaml',
    params:
        genome = run.genome_build,
        genomes_dir = hpc.genomes_dir,
        work_dir = 'work/{batch}/small_variants',
        unlock_opt = ' --unlock' if config.get('unlock', 'no') == 'yes' else '',
        tumor_sample = lambda wc: batch_by_name[wc.batch].tumor.rgid,
        normal_sample = lambda wc: batch_by_name[wc.batch].normal.rgid,
    resources:
        mem_mb = 20000
    group: "somatic_anno"
    shell:
        'anno_somatic_vcf {input.vcf} -o {output.vcf} '
        '-tn {params.tumor_sample} -nn {params.normal_sample} '
        '-w {params.work_dir} -g {params.genome} --genomes-dir {params.genomes_dir} '
        '{params.unlock_opt}'

rule somatic_vcf_filter:
    input:
        vcf = rules.somatic_vcf_annotate.output.vcf,
    output:
        vcf = '{batch}/small_variants/{batch}-somatic.vcf.gz',
    group: "somatic_filt"
    params:
        tumor_sample = lambda wc: batch_by_name[wc.batch].tumor.rgid,
        normal_sample = lambda wc: batch_by_name[wc.batch].normal.rgid,
    shell:
        'filter_somatic_vcf {input.vcf} -o {output.vcf}'
         ' -tn {params.tumor_sample} -nn {params.normal_sample}'

rule somatic_vcf_filter_pass:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic.PASS.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-somatic.PASS.vcf.gz.tbi',
    group: "somatic_filt"
    shell:
        'bcftools view -f.,PASS {input.vcf} -Oz -o {output.vcf} && tabix -f -p vcf {output.vcf}'

rule bcftools_stats_somatic:
    input:
        rules.somatic_vcf_filter_pass.output.vcf
    output:
        '{batch}/small_variants/stats/{batch}_bcftools_stats.txt'
    group: "somatic_filt"
    params:
        sname = lambda wc: batch_by_name[wc.batch].tumor.rgid,
    shell:
        'bcftools stats -s {params.sname} {input} | sed s#{input}#{params.sname}# > {output}'

##################
#### Germline ####
include_germline = all(b.germline_vcf for b in batch_by_name.values())
if include_germline:
    rule germline_vcf_pass:
        input:
            vcf = lambda wc: batch_by_name[wc.batch].germline_vcf,
        output:
            vcf = 'work/{batch}/small_variants/germline/{batch}-germline.PASS.vcf.gz',
        group: "germline_snv"
        shell:
            'bcftools view {input.vcf} -f.,PASS -Oz -o {output.vcf} && tabix -f -p vcf {output.vcf}'

    # Subset to ~200 cancer predisposition gene set.
    rule germline_predispose_subset:
        input:
            vcf = rules.germline_vcf_pass.output.vcf,
            predispose_genes_bed = get_predispose_genes_bed(run.genome_build, coding_only=True),
        output:
            vcf = 'work/{batch}/small_variants/germline/{batch}-germline.predispose_genes.vcf.gz',
            tbi = 'work/{batch}/small_variants/germline/{batch}-germline.predispose_genes.vcf.gz.tbi',
        params:
            ungz = lambda wc, output: get_ungz_gz(output[0])[0]
        group: "germline_snv"
        run:
            shell('bcftools view -T {input.predispose_genes_bed} {input.vcf} -Oz -o {output.vcf}'
                  ' && tabix -p vcf {output.vcf}')
            cnt = int(subprocess.check_output(f'bcftools view -H {output.vcf} | wc -l', shell=True).strip())
            if cnt == 0:
                critical(f'Found zero germline variants within predisposition genes in {output.vcf}')

    # Preparations: annotate TUMOR_X and NORMAL_X fields
    # Used for PCGR, but for all other processing steps too
    # rule germline_predispose_subset_vcf_prep:
    #     input:
    #         vcf = rules.germline_predispose_subset.output.vcf
    #     output:
    #         vcf = 'work/{batch}/small_variants/germline/{batch}-normal-' + run.germline_caller + '-predispose_genes.vcf.gz',
    #         tbi = 'work/{batch}/small_variants/germline/{batch}-normal-' + run.germline_caller + '-predispose_genes.vcf.gz.tbi',
    #     group: "germline_snv"
    #     shell:
    #         'pcgr_prep {input.vcf} | bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'

    # # TODO: merge with filtered out somatic varians.
    # #  1. PoN, normal fail     - low freq (< 1/3 of purity) - remove; otherwise add to germline
    # #  2. PoN, normal ok       - low freq (< 1/3 of purity) - remove; any normal support - add; otherwise remove
    # #  3. gnomAD, normal fail  - add to germline
    # #  4. gnomAD, normal ok    - add to germline if has any normal support
    rule germline_leakage:
        input:
            vcf = rules.somatic_vcf_filter.output.vcf
        output:
            vcf = 'work/{batch}/small_variants/germline/{batch}-tumor-germline-leakage.vcf.gz',
        group: "germline_snv"
        params:
            toremove = 'INFO/AC,INFO/AF,INFO/TUMOR_AF,INFO/TUMOR_VD,INFO_TUMOR_MQ,INFO/TUMOR_DP,FORMAT/AD,FORMAT/ADJAF,FORMAT/AF,FORMAT/VD,FILTER',
            tumor_sample = lambda wc: batch_by_name[wc.batch].tumor.rgid,
            normal_sample = lambda wc: batch_by_name[wc.batch].normal.rgid,
        shell:
            'bcftools filter -i "Germline=1" {input.vcf} | '
            'bcftools annotate -x "{params.toremove}" | ' 
            'bcftools view -s {params.tumor_sample} | ' 
            'sed \'s/{params.tumor_sample}/{params.normal_sample}/\' | ' 
            'bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'

    # Subset to ~200 cancer predisposition gene set.
    rule germline_leakage_predispose_subset:
        input:
            vcf = rules.germline_leakage.output.vcf,
        output:
            vcf = 'work/{batch}/small_variants/germline/{batch}-tumor-germline-leakage-predispose_genes.vcf.gz',
            tbi = 'work/{batch}/small_variants/germline/{batch}-tumor-germline-leakage-predispose_genes.vcf.gz.tbi',
        params:
            ungz = lambda wc, output: get_ungz_gz(output[0])[0]
        group: "germline_snv"
        run:
            pcgr_toml_fpath = join(package_path(), 'pcgr', 'cpsr.toml')
            genes = [g for g in toml.load(pcgr_toml_fpath)['cancer_predisposition_genes']]
            def func(rec, vcf):
                if rec.INFO.get('PCGR_SYMBOL') is not None and rec.INFO['PCGR_SYMBOL'] in genes:
                    return rec
            iter_vcf(input.vcf, output.vcf, func)

    rule germline_merge_with_leakage:
        input:
            vcfs = [
                rules.germline_predispose_subset.output.vcf,
                rules.germline_leakage_predispose_subset.output.vcf,
            ]
        output:
            vcf = 'work/{batch}/small_variants/{batch}-germline.predispose_genes.vcf.gz',
        group: "germline_snv"
        shell:
            'bcftools concat -a {input.vcfs} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'


rule somatic_stats_report:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf,
        full_vcf = rules.somatic_vcf_select_noalt.output.vcf,
        subset_highly_mutated_stats = rules.somatic_vcf_annotate.output.subset_highly_mutated_stats,
    output:
        'work/{batch}/small_variants/{batch}_somatic_stats.yml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name,
    group: "somatic_filt"
    run:
        snps, indels, others = cnt_vars(input.vcf, passed=True)
        all_snps, all_indels, all_others = cnt_vars(input.full_vcf)
        subset_genes = None

        with open(input.subset_highly_mutated_stats) as inp:
            stats = yaml.load(inp)
            total_vars = stats['total_vars']
            vars_no_gnomad = stats.get('vars_no_gnomad')
            vars_cancer_genes = stats.get('vars_cancer_genes')
            if vars_cancer_genes is not None:
                # In this case, "all_snps", "all_indels", "all_others" correspond to pre-cancer-gene subset variants
                # Thus we can't compare it with out final cancer-gene-subset PASSing variants.
                subset_genes = (total_vars - vars_cancer_genes) / total_vars * 100.0
                all_snps, all_indels, all_others = cnt_vars(input.vcf, passed=False)

        all_vars = all_snps + all_indels + all_others
        vars = snps + indels + others

        data = dict(
            snps=snps,
            indels=indels,
            others=others if others else None,
            filt_vars=(all_vars - vars) / all_vars * 100.0,
            filt_snps=(all_snps - snps) / all_snps * 100.0,
            filt_indels=(all_indels - indels) / all_indels * 100.0,
            filt_others=(((all_others - others) / all_others * 100.0) if others and all_others else None),
        )
        if subset_genes is not None:
            data.update(dict(
                subset_genes=subset_genes
            ))

        with open(output[0], 'w') as out:
            data = {
                'data': {
                    params.sample: data
                }
            }
            yaml.dump(data, out, default_flow_style=False)

if include_germline:
    rule germline_stats_report:
        input:
            pass_vcf = rules.germline_vcf_pass.output.vcf,
            pass_predispose_vcf = rules.germline_merge_with_leakage.output.vcf \
                                  if 'somatic' in stages \
                                  else rules.germline_predispose_subset.output.vcf,
        output:
            'work/{batch}/small_variants/{batch}_germline_stats.yml',
        params:
            sample = lambda wc: batch_by_name[wc.batch].normal.name
        group: "germline_snv"
        run:
            pass_cnt = int(subprocess.check_output(
                f'bcftools view -H {input.pass_vcf} | wc -l', shell=True).strip())
            predispose_pass_cnt = int(subprocess.check_output(
                f'bcftools view -H {input.pass_predispose_vcf} | wc -l', shell=True).strip())

            with open(output[0], 'w') as out:
                data = {
                    'data': {
                        params.sample: dict(
                            germline = pass_cnt,
                            germline_predispose = predispose_pass_cnt
                        )
                    }
                }
                yaml.dump(data, out, default_flow_style=False)

    # Produces the same stats as bcbio file in QC, but have to rerun because of CWL version,
    # which doesn't suffix the file with _germline, so MultiQC can't relate it to the germline
    # stats section.
    rule bcftools_stats_germline:
        input:
            rules.germline_merge_with_leakage.output.vcf \
            if 'somatic' in stages else \
            rules.germline_predispose_subset.output.vcf,
        output:
            '{batch}/small_variants/stats/{batch}_bcftools_stats_germline.txt'
        group: "germline_snv"
        params:
            sname = lambda wc: batch_by_name[wc.batch].normal.rgid,
        shell:
            'bcftools stats -s {params.sname} {input} | sed s#{input}#{params.sname}# > {output}'

    rule copy_germline:
        input:
            germline_vcfs = expand(rules.germline_merge_with_leakage.output.vcf
                                   if 'somatic' in stages else \
                                   rules.germline_predispose_subset.output.vcf, batch=batch_by_name.keys())
        output:
            temp(touch('work/copy_germline.done'))
        run:
            for bn, batch, germline_vcf in zip(batch_by_name.keys(), batch_by_name.values(), input.germline_vcfs):
                assert batch.name in germline_vcf
                shell(f'mkdir -p {join(bn, "small_variants")}')
                renamed_germline_vcf = join(bn, 'small_variants',
                         f'{batch}__{batch.normal.name}-germline.predispose_genes.vcf.gz')
                shell(f'cp {germline_vcf} {renamed_germline_vcf}')
                shell(f'tabix -p vcf {renamed_germline_vcf}')


#############

rule germline:
    input:
        germline = 'work/copy_germline.done'
    output:
        temp(touch('log/germline.done'))


rule somatic:
    input:
        somatic_vcfs = expand(rules.somatic_vcf_filter.output.vcf, batch=batch_by_name.keys()),
        somatic_pass_vcfs = expand(rules.somatic_vcf_filter_pass.output.vcf, batch=batch_by_name.keys()),
    output:
        temp(touch('log/somatic.done'))


rule small_variants:
    input:
        'log/germline.done',
        'log/somatic.done',
    output:
        temp(touch('log/small_variants.done'))




