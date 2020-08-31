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
from umccrise import package_path


localrules: sage2


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


rule run_sage:
    input:
        tumor_bams  = lambda wc: [s.bam for s in batch_by_name[wc.batch].tumors],
        normal_bams = lambda wc: [s.bam for s in batch_by_name[wc.batch].normals] +\
                                 [s.bam for s in batch_by_name[wc.batch].rna_samples],
        ref_fa        = refdata.get_ref_file(run.genome_build, key='fa'),
        hotspots_vcf  = refdata.get_ref_file(run.genome_build, key='hotspots'),
        coding_bed    = refdata.get_ref_file(run.genome_build, key='coding_regions'),
        high_conf_bed = refdata.get_ref_file(run.genome_build, key='hmf_giab_conf'),
    output:
        vcf = 'work/{batch}/sage/1_run/{batch}.vcf.gz',
        tbi = 'work/{batch}/sage/1_run/{batch}.vcf.gz.tbi',
    params:
        jar = join(package_path(), 'sage-2.2.jar'),
        xms = 4000,
        xmx = 28000,
        genome = run.genome_build,
        tumor_names  = lambda wc: [s.name for s in batch_by_name[wc.batch].tumors],
        normal_names = lambda wc: [s.name for s in batch_by_name[wc.batch].normals] +\
                                  [s.name for s in batch_by_name[wc.batch].rna_samples],
    resources:
        mem_mb = lambda wildcards, attempt: 30000 * attempt
    threads: threads_per_batch,
    group: 'sage'
    run:
        shell(conda_cmd.format('hmf') +
            f'java -Xms{params.xms}m -Xmx{params.xmx}m -cp {params.jar} '
            f'com.hartwig.hmftools.sage.SageApplication '
            f'-threads {threads} '
            f'-tumor {",".join(params.tumor_names)} -tumor_bam {",".join(input.tumor_bams)} '
            f'-reference {",".join(params.normal_names)} -reference_bam {",".join(input.normal_bams)} '
            f'-hotspots {input.hotspots_vcf} '
            f'-panel_bed {input.coding_bed} '
            f'-high_confidence_bed {input.high_conf_bed} '
            f'-assembly {params.genome} '
            f'-ref_genome {input.ref_fa} '
            f'-out {output.vcf} '
            f'&& tabix -f -p vcf {output.vcf}'
        )


rule sage_pass:
    input:
        vcf = 'work/{batch}/sage/1_run/{batch}.vcf.gz',
    output:
        vcf = 'work/{batch}/sage/2_pass/{batch}.vcf.gz',
        tbi = 'work/{batch}/sage/2_pass/{batch}.vcf.gz.tbi',
    group: 'sage'
    shell:
        'bcftools view -f.,PASS {input.vcf} -Oz -o {output.vcf} '
        '&& tabix -p vcf {output.vcf}'


rule sage_pon:
    input:
        vcf = 'work/{batch}/sage/2_pass/{batch}.vcf.gz',
        hmf_pon = refdata.get_ref_file(run.genome_build, key='hmf_pon'),
    output:
        vcf = 'work/{batch}/sage/3_pon/{batch}.vcf.gz',
    group: 'sage'
    shell:
        # PON_MAX is the maximum allelic variant depth for a PoN variant
        # PON_COUNT is the number of hits in the PoN
        'bcftools annotate -a {input.hmf_pon} '
        '-c PON_COUNT,PON_MAX {input.vcf} -Oz -o {output.vcf}'


rule sage_somatic_anno:
    input:
        vcf = 'work/{batch}/sage/3_pon/{batch}.vcf.gz',
    output:
        vcf = 'work/{batch}/sage/4_somatic_anno/{batch}.vcf.gz',
    params:
        genome = run.genome_build,
        genomes_dir = refdata.genomes_dir,
        work_dir = 'work/{batch}/sage/4_somatic_anno',
        unlock_opt = ' --unlock' if config.get('unlock', 'no') == 'yes' else '',
        tumor_names  = lambda wc: ','.join(s.name for s in batch_by_name[wc.batch].tumors),
    resources:
        mem_mb = 20000
    group: 'somatic_anno'
    shell:
        'anno_somatic_vcf {input.vcf} -o {output.vcf} '
        '-tn {params.tumor_names} '
        '-w {params.work_dir} -g {params.genome} --genomes-dir {params.genomes_dir} '
        '{params.unlock_opt}'


rule sage_somatic_filt:
    input:
        vcf = 'work/{batch}/sage/4_somatic_anno/{batch}.vcf.gz',
    output:
        vcf = '{batch}/sage/{batch}.vcf.gz',
        tbi = '{batch}/sage/{batch}.vcf.gz.tbi',
    params:
        min_af = 0.01,
        tumor_names  = lambda wc: ','.join(s.name for s in batch_by_name[wc.batch].tumors),
    group: 'somatic_anno'
    shell:
        'filter_somatic_vcf -tn {params.tumor_names} {input.vcf} '
        '-o {output.vcf} --min-af {params.min_af}'


rule sage_somatic_filt_pass:
    input:
        vcf = '{batch}/sage/{batch}.vcf.gz',
    output:
        vcf = '{batch}/sage/{batch}.pass.vcf.gz',
        tbi = '{batch}/sage/{batch}.pass.vcf.gz.tbi',
    group: 'somatic_anno'
    shell:
        'bcftools view -f.,PASS {input.vcf} -Oz -o {output.vcf} '
        '&& tabix -p vcf {output.vcf}'


rule sage_somatic_stats_report:
    input:
        full_vcf = '{batch}/sage/{batch}.vcf.gz',
        vcf = '{batch}/sage/{batch}.pass.vcf.gz',
    output:
        'work/{batch}/sage/{batch}_somatic_stats.yml',
    group: 'somatic_anno'
    run:
        snps, indels, others = cnt_vars(input.vcf, passed=True)
        all_snps, all_indels, all_others = cnt_vars(input.full_vcf)
        subset_genes = None

        all_vars = all_snps + all_indels + all_others
        vars = snps + indels + others

        data = dict(
            snps=snps,
            indels=indels,
            others=others if others else None,
            filt_vars=(all_vars - vars) / all_vars * 100.0,
            filt_snps=(all_snps - snps) / all_snps * 100.0,
            filt_indels=(all_indels - indels) / all_indels * 100.0,
            filt_others=(((all_others - others) / all_others * 100.0)
                         if others and all_others else None),
        )
        if subset_genes is not None:
            data.update(dict(
                subset_genes=subset_genes
            ))

        with open(output[0], 'w') as out:
            data = {
                'data': {
                    wildcards.batch: data
                }
            }
            yaml.dump(data, out, default_flow_style=False)


rule sage_vcf2maf:
    input:
        vcf = '{batch}/sage/{batch}.pass.vcf.gz',
        fa = refdata.get_ref_file(genome=run.genome_build, key='fa')
    output:
        maf = '{batch}/sage/{batch}.pass.maf',
    params:
        tname  = lambda wc: [s.name for s in batch_by_name[wc.batch].tumors][0],
        nname = lambda wc: [s.name for s in batch_by_name[wc.batch].normals][0],  # +\
#                           [s.name for s in batch_by_name[wc.batch].rna_samples],
        ncbi_build = {'hg38': 'GRCh38', 'GRCh37': 'GRCh37'}.get(run.genome_build),
        uncompressed_tmp_vcf = '{batch}/sage/{batch}.pass.vcf.tmp'
    shell:
        'gunzip -c {input.vcf} > {params.uncompressed_tmp_vcf} '
        '&& vcf2maf.pl --inhibit-vep ' 
        '--input-vcf {params.uncompressed_tmp_vcf} '
        '--output-maf {output.maf} '
        '--ref-fasta {input.fa} '
        '--filter-vcf 0 '
        '--tumor-id {params.tname} '
        '--normal-id {params.nname} '
        '--ncbi-build {params.ncbi_build} '
        '&& rm {params.uncompressed_tmp_vcf}'


rule sage2:
    input:
        vcf = expand('{batch}/sage/{batch}.pass.vcf.gz', batch=batch_by_name.keys()),
        maf = expand('{batch}/sage/{batch}.pass.maf', batch=batch_by_name.keys()),
    output:
        temp(touch('log/sage2.done'))




