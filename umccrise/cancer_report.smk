import shutil
import yaml
from os.path import join, abspath, dirname
from ngs_utils.logger import warn
from ngs_utils.reference_data import get_key_genes, get_key_genes_bed
from ngs_utils.file_utils import safe_mkdir
from ngs_utils.vcf_utils import count_vars
from reference_data import api as refdata

from umccrise import package_path


localrules: cancer_report, conda_list


# Subset final SNVs to GiaB confident intervals
rule subset_to_giab:
    input:
        vcf = '{batch}/small_variants/{batch}-somatic-PASS.vcf.gz',
    params:
        regions = refdata.get_ref_file(run.genome_build, key=['hmf_giab_conf'])
    output:
        'work/{batch}/cancer_report/afs/somatic-confident.vcf.gz'
    group: "rmd_prep"
    shell:
        'bcftools view {input.vcf} -T <(gunzip -c {params.regions}) -Oz -o {output}'

# Keep INFO/TUMOR_AF and split multiallelics (to avoid R parsing issues)
rule split_multiallelic:
    input:
        vcf = 'work/{batch}/cancer_report/afs/somatic-confident.vcf.gz',
        ref_fa = refdata.get_ref_file(run.genome_build, key='fa')
    output:
        'work/{batch}/cancer_report/afs/somatic-confident-singleallelic.vcf.gz'
    group: "rmd_prep"
    shell:
        'bcftools annotate -x ^INFO/TUMOR_AF {input.vcf} -Ob | '
        'bcftools norm -m \'-\' -f {input.ref_fa} -Ob | '
        'bcftools sort -Oz -o {output} '
        '&& tabix -p vcf {output}'

# Extract just INFO/TUMOR_AF
rule afs:
    input:
        'work/{batch}/cancer_report/afs/somatic-confident-singleallelic.vcf.gz'
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumors[0].rgid
    output:
        'work/{batch}/cancer_report/afs/af_tumor.txt'
    group: "rmd_prep"
    shell:
        'bcftools view {input} -s {params.tumor_name} -Ou | '
        '(printf "af\n"; bcftools query -f "%INFO/TUMOR_AF\\n") > {output} '
        '&& test -e {output}'

# Subset to just key genes
rule afs_keygenes:
    input:
        vcf = 'work/{batch}/cancer_report/afs/somatic-confident-singleallelic.vcf.gz',
        bed = get_key_genes_bed(run.genome_build, coding_only=True),
    params:
        tumor_name = lambda wc: batch_by_name[wc.batch].tumors[0].rgid
    output:
        'work/{batch}/cancer_report/afs/af_tumor_keygenes.txt'
    group: "rmd_prep"
    shell:
        'bcftools view -f .,PASS {input.vcf} -s {params.tumor_name} -Ov'
        ' | bedtools intersect -a stdin -b {input.bed} -header'
        ' | (printf "chrom\\tpos\\tid\\tref\\talt\\taf\\n" ; bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/TUMOR_AF\\n")'
        ' > {output} && test -e {output}'

# List all conda pkgs per env
rule conda_list:
    output:
        txt = 'work/{batch}/conda_pkg_list.txt'
    params:
        env = ['""', '_pcgr', '_hmf', '_cancer_report'],
    group: "rmd_prep"
    shell:
        "for e in {params.env}; do conda list -p {env_path}$e "
        "| grep -v ^# "
        "| awk -v var=env$e '{{ print var, $0 }}' >> {output} ; done"

rule somatic_snv_summary:
    input:
        raw = lambda wc: batch_by_name[wc.batch].somatic_vcf,
        raw_pass = 'work/{batch}/small_variants/pass_sort/{batch}-somatic.vcf.gz',
        noalt = 'work/{batch}/small_variants/noalt/{batch}-somatic.vcf.gz',
        sage = 'work/{batch}/small_variants/sage1/{batch}-somatic.vcf.gz',
        subset_highly_mutated_stats = 'work/{batch}/small_variants/somatic_anno/subset_highly_mutated_stats.yaml', # this has sage + gnomad counts
        anno = 'work/{batch}/small_variants/annotate/{batch}-somatic.vcf.gz',
        filt = '{batch}/small_variants/{batch}-somatic.vcf.gz',
        filt_pass = '{batch}/small_variants/{batch}-somatic-PASS.vcf.gz',
        giab = 'work/{batch}/cancer_report/afs/somatic-confident.vcf.gz'
    output:
        yaml = 'work/{batch}/cancer_report/somatic_snv_summary.yaml'
    run:
        stats = dict(
          raw = count_vars(input.raw),
          raw_pass = count_vars(input.raw_pass),
          noalt = count_vars(input.noalt),
          sage = count_vars(input.sage),
          anno = count_vars(input.anno),
          filt = count_vars(input.filt),
          filt_pass = count_vars(input.filt_pass),
          giab = count_vars(input.giab)
        )

        with open(output.yaml, 'w') as out:
            yaml.dump(stats, out, default_flow_style=False)



rule run_cancer_report:
    input:
        key_genes            = get_key_genes(),
        af_global            = rules.afs.output[0],
        af_keygenes          = rules.afs_keygenes.output[0],
        somatic_snv_summary  = rules.somatic_snv_summary.output.yaml,
        somatic_snv_vcf      = '{batch}/small_variants/{batch}-somatic-PASS.vcf.gz',
        somatic_sv_tsv       = lambda wc: rules.prep_sv_tsv.output[0]
                               if (batch_by_name[wc.batch].sv_vcf and 'structural' in stages) else [],
        somatic_sv_vcf       = lambda wc: '{batch}/structural/{batch}-manta.vcf.gz'
                               if (batch_by_name[wc.batch].sv_vcf and 'structural' in stages) else [],
        purple_som_snv_vcf   = 'work/{batch}/purple/{batch}.purple.somatic.vcf.gz',
        purple_som_cnv       = 'work/{batch}/purple/{batch}.purple.cnv.somatic.tsv',
        purple_som_gene_cnv  = 'work/{batch}/purple/{batch}.purple.cnv.gene.tsv',
        purple_germ_cnv      = 'work/{batch}/purple/{batch}.purple.cnv.germline.tsv',
        purple_purity        = 'work/{batch}/purple/{batch}.purple.purity.tsv',
        purple_qc            = 'work/{batch}/purple/{batch}.purple.qc',

        purple_circos_png    = 'work/{batch}/purple/plot/{batch}.circos.png',
        purple_input_png     = 'work/{batch}/purple/plot/{batch}.input.png',
        purple_cn_png        = 'work/{batch}/purple/plot/{batch}.copynumber.png',
        purple_ma_png        = 'work/{batch}/purple/plot/{batch}.map.png',
        purple_purity_png    = 'work/{batch}/purple/plot/{batch}.purity.range.png',
        purple_segment_png   = 'work/{batch}/purple/plot/{batch}.segment.png',
        purple_clonality_png = 'work/{batch}/purple/plot/{batch}.somatic.clonality.png',
        purple_ploidy_png    = 'work/{batch}/purple/plot/{batch}.somatic.png',
        purple_baf_png       = 'work/{batch}/purple/circos_baf/{batch}.circos_baf.png',

        wait_for_integration_sites = get_integration_sites_tsv_fn
             if 'oncoviruses' in stages else [],
        oncoviral_present_viruses = 'work/{batch}/oncoviruses/present_viruses.txt'
             if 'oncoviruses' in stages else [],

        conda_list           = rules.conda_list.output.txt,

    params:
        tumor_name          = lambda wc: batch_by_name[wc.batch].tumors[0].rgid,
        result_outdir       = lambda wc, output: join(os.getcwd(), dirname(output[0]), 'cancer_report_tables'),
        output_file         = lambda wc, output: join(os.getcwd(), output[0]),
        af_global           = lambda wc, input: abspath(input.af_global),
        af_keygenes         = lambda wc, input: abspath(input.af_keygenes),
        somatic_snv_summary = lambda wc, input: abspath(input.somatic_snv_summary),
        somatic_snv_vcf     = lambda wc, input: abspath(input.somatic_snv_vcf),
        somatic_sv_tsv      = lambda wc, input: abspath(input.somatic_sv_tsv) if input.somatic_sv_tsv else 'NA',
        somatic_sv_vcf      = lambda wc, input: abspath(input.somatic_sv_vcf) if input.somatic_sv_vcf else 'NA',
        purple_som_snv_vcf  = lambda wc, input: abspath(input.purple_som_snv_vcf),
        purple_som_cnv      = lambda wc, input: abspath(input.purple_som_cnv),
        purple_som_gene_cnv = lambda wc, input: abspath(input.purple_som_gene_cnv),
        purple_germ_cnv     = lambda wc, input: abspath(input.purple_germ_cnv),
        purple_purity       = lambda wc, input: abspath(input.purple_purity),
        purple_qc           = lambda wc, input: abspath(input.purple_qc),
        conda_list          = lambda wc, input: abspath(input.conda_list),
        img_dir_abs         = lambda wc, output: abspath(output.img_dir),
    output:
        report_html = '{batch}/{batch}_cancer_report.html',
        img_dir = directory('work/{batch}/cancer_report/img'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000
    run:
        ov_cmdl = ''
        if 'oncoviruses' in stages:
            oncoviral_present_viruses = abspath(input.oncoviral_present_viruses)
            oncoviral_breakpoints_tsv = ''
            with open(input.oncoviral_present_viruses) as f:
                if [v for v in f.read().strip().split(',') if v]:
                    oncoviral_breakpoints_tsv = \
                        abspath(f'work/{wildcards.batch}/oncoviruses/oncoviral_breakpoints.tsv')
            ov_cmdl = (
                f"--oncoviral_present_viruses '{oncoviral_present_viruses}' "
                f"--oncoviral_breakpoints_tsv '{oncoviral_breakpoints_tsv}'"
            )

        # copy BAF circos + PURPLE plots to img dir
        shell('mkdir -p {output.img_dir}')
        shell('cp {input.purple_baf_png} {output.img_dir}')
        shell('cp $(dirname {input.purple_circos_png})/* {output.img_dir}')

        shell(conda_cmd.format('cancer_report') + """
gpgr.R canrep \
  --af_global '{params.af_global}' \
  --af_keygenes '{params.af_keygenes}' \
  --batch_name '{wildcards.batch}' \
  --conda_list '{params.conda_list}' \
  --img_dir '{params.img_dir_abs}' \
  --key_genes '{input.key_genes}' \
  --somatic_snv_summary '{params.somatic_snv_summary}' \
  --somatic_snv_vcf '{params.somatic_snv_vcf}' \
  --somatic_sv_tsv '{params.somatic_sv_tsv}' \
  --somatic_sv_vcf '{params.somatic_sv_vcf}' \
  --purple_som_gene_cnv '{params.purple_som_gene_cnv}' \
  --purple_som_cnv '{params.purple_som_cnv}' \
  --purple_germ_cnv '{params.purple_germ_cnv}' \
  --purple_purity '{params.purple_purity}' \
  --purple_qc '{params.purple_qc}' \
  --purple_som_snv_vcf '{params.purple_som_snv_vcf}' \
  {ov_cmdl} \
  --out_file '{params.output_file}' \
  --result_outdir '{params.result_outdir}' \
  --tumor_name '{params.tumor_name}'
""")


#############

rule cancer_report:
    input:
        expand(rules.run_cancer_report.output[0], batch=batch_by_name.keys())
    output:
        temp(touch('log/cancer_report.done'))

