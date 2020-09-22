#################
#### Somatic ####
import subprocess
from os.path import isfile, join, dirname
import toml
import cyvcf2
import yaml
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.reference_data import get_predispose_genes_bed
from ngs_utils.logger import critical
from ngs_utils.vcf_utils import iter_vcf
from reference_data import api as refdata
from umccrise import package_path, cnt_vars


localrules: germline, germline_batch


rule germline_vcf_pass:
    input:
        vcf = lambda wc: batch_by_name[wc.batch].germline_vcf or [],
    output:
        vcf = 'work/{batch}/small_variants/germline/{batch}-germline-PASS.vcf.gz',
    group: "germline_snv"
    shell:
        'bcftools view {input.vcf} -f.,PASS -Oz -o {output.vcf} && tabix -f -p vcf {output.vcf}'

# Subset to ~200 cancer predisposition gene set.
rule germline_predispose_subset:
    input:
        vcf = rules.germline_vcf_pass.output.vcf,
        predispose_genes_bed = get_predispose_genes_bed(run.genome_build, coding_only=False),
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
        tumor_sample = lambda wc: batch_by_name[wc.batch].tumors[0].rgid,
        normal_sample = lambda wc: batch_by_name[wc.batch].normals[0].rgid,
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
        predispose_genes_bed = get_predispose_genes_bed(run.genome_build, coding_only=False),
    output:
        vcf = 'work/{batch}/small_variants/germline/{batch}-tumor-germline-leakage-predispose_genes.vcf.gz',
        tbi = 'work/{batch}/small_variants/germline/{batch}-tumor-germline-leakage-predispose_genes.vcf.gz.tbi',
    params:
        ungz = lambda wc, output: get_ungz_gz(output[0])[0]
    group: "germline_snv"
    shell:
        'bcftools view -T {input.predispose_genes_bed} {input.vcf} -Oz -o {output.vcf}'
        ' && tabix -p vcf {output.vcf}'

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

rule germline_stats_report:
    input:
        pass_vcf = rules.germline_vcf_pass.output.vcf,
        pass_predispose_vcf = rules.germline_merge_with_leakage.output.vcf \
                              if 'somatic' in stages \
                              else rules.germline_predispose_subset.output.vcf,
    output:
        'work/{batch}/small_variants/{batch}_germline_stats.yml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].normals[0].name
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
        sname = lambda wc: batch_by_name[wc.batch].normals[0].rgid,
    shell:
        'bcftools stats -s {params.sname} {input} | sed s#{input}#{params.sname}# > {output}'


#############

rule germline_batch:
    input:
        lambda wc: [] if not batch_by_name[wc.batch].germline_vcf \
                  else rules.germline_merge_with_leakage.output.vcf \
                       if 'somatic' in stages else \
                       rules.germline_predispose_subset.output.vcf
    output:
        temp(touch('log/germline_{batch}.done'))
    run:
        for bn, batch, germline_vcf in zip(batch_by_name.keys(),
                                           batch_by_name.values(),
                                           input):
            shell(f'mkdir -p {join(bn, "small_variants")}')
            renamed_germline_vcf = join(bn, 'small_variants',
                     f'{batch}__{batch.normals[0].name}-germline.predispose_genes.vcf.gz')
            shell(f'cp {germline_vcf} {renamed_germline_vcf}')
            shell(f'tabix -p vcf {renamed_germline_vcf}')

# rule copy_germline:
#     input:
#         germline_vcfs = expand(rules.germline_merge_with_leakage.output.vcf
#                                if 'somatic' in stages else \
#                                rules.germline_predispose_subset.output.vcf, batch=batch_by_name.keys())
#     output:
#         temp(touch('log/germline.done'))
#     run:
#         for bn, batch, germline_vcf in zip(batch_by_name.keys(), batch_by_name.values(), input.germline_vcfs):
#             assert batch.name in germline_vcf
#             shell(f'mkdir -p {join(bn, "small_variants")}')
#             renamed_germline_vcf = join(bn, 'small_variants',
#                      f'{batch}__{batch.normal.name}-germline.predispose_genes.vcf.gz')
#             shell(f'cp {germline_vcf} {renamed_germline_vcf}')
#             shell(f'tabix -p vcf {renamed_germline_vcf}')

rule germline:
    input:
        germline = expand(rules.germline_batch.output, batch=batch_by_name.keys())
    output:
        temp(touch('log/germline.done'))




