#################
#### Somatic ####
import subprocess
from os.path import isfile, join, dirname
import toml
import cyvcf2
import yaml
import shutil
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.reference_data import get_predispose_genes_bed
from ngs_utils.logger import critical
from ngs_utils.vcf_utils import iter_vcf
from reference_data import api as refdata
from umccrise import package_path, cnt_vars


localrules: somatic, germline, germline_batch, pierian


# rule somatic_vcf_reheader  # change RGIDs to tumor and normal names?

rule run_sage:
    input:
        tumor_bams  = lambda wc: [s.bam          for s in batch_by_name[wc.batch].tumors],
        tumor_bais  = lambda wc: [s.bam + '.bai' for s in batch_by_name[wc.batch].tumors],
        normal_bams = lambda wc: [s.bam          for s in batch_by_name[wc.batch].normals],
        normal_bais = lambda wc: [s.bam + '.bai' for s in batch_by_name[wc.batch].normals],
        rna_bams    = lambda wc: [s.bam          for s in batch_by_name[wc.batch].rna_samples],
        rna_bais    = lambda wc: [s.bam + '.bai' for s in batch_by_name[wc.batch].rna_samples],

        ref_fa        = refdata.get_ref_file(run.genome_build, key='fa'),
        hotspots_vcf  = refdata.get_ref_file(run.genome_build, key='hotspots'),
        coding_bed    = refdata.get_ref_file(run.genome_build, key='coding_regions'),
        high_conf_bed = refdata.get_ref_file(run.genome_build, key='hmf_giab_conf'),
    output:
        vcf = 'work/{batch}/small_variants/sage2/{batch}.vcf.gz',
        tbi = 'work/{batch}/small_variants/sage2/{batch}.vcf.gz.tbi',
    params:
        jar = join(package_path(), 'sage-2.2.jar'),
        xms = 4000,
        xmx = 31000,
        genome = run.genome_build,
    resources:
        mem_mb = 50000
    threads: threads_per_batch,
    group: 'sage'
    run:
        tumor_bams   = input.tumor_bams
        normal_bams  = input.normal_bams
        tumor_names  = [s.rgid for s in batch_by_name[wildcards.batch].tumors]
        normal_names = [s.rgid for s in batch_by_name[wildcards.batch].normals]

        rna_bam = ''
        rna_name = ''
        if input.rna_bams:
            rna_bam = input.rna_bams[0]
            rna_name = wildcards.batch + '_rna'
            normal_names.append(rna_name)
            normal_bams.append(rna_bam)

        shell(conda_cmd.format('hmf') +
            f'java -Xms{params.xms}m -Xmx{params.xmx}m -cp {params.jar} '
            f'com.hartwig.hmftools.sage.SageApplication '
            f'-threads {threads} '
            f'-tumor {",".join(tumor_names)} -tumor_bam {",".join(tumor_bams)} '
            f'-reference {",".join(normal_names)} -reference_bam {",".join(normal_bams)} '
            f'-hotspots {input.hotspots_vcf} '
            f'-panel_bed {input.coding_bed} '
            f'-high_confidence_bed {input.high_conf_bed} '
            f'-assembly {params.genome} '
            f'-ref_genome {input.ref_fa} '
            f'-out {output.vcf} '
            f'&& test -s {output.vcf} '
            f'&& tabix -f -p vcf {output.vcf}'
        )
        verify_file(output.vcf)

rule somatic_vcf_pass_sort:
    input:
        vcf = lambda wc: batch_by_name[wc.batch].somatic_vcf or
                         f'work/{wc.batch}/small_variants/sage2/{wc.batch}.vcf.gz',
    output:
        vcf = 'work/{batch}/small_variants/pass_sort/{batch}-somatic.vcf.gz',
        tbi = 'work/{batch}/small_variants/pass_sort/{batch}-somatic.vcf.gz.tbi',
    group: "somatic_anno"
    shell:
        'test -s {input.vcf} && '
        '(bcftools view -h {input.vcf} ; '
        'bcftools view -H -f.,PASS {input.vcf} | sort -k1,1V -k2,2n) | '
        'bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'

rule somatic_vcf_select_noalt:
    input:
        vcf = rules.somatic_vcf_pass_sort.output.vcf,
        noalts_bed = refdata.get_ref_file(run.genome_build, 'noalt_bed'),
    output:
        vcf = 'work/{batch}/small_variants/noalt/{batch}-somatic.vcf.gz',
        tbi = 'work/{batch}/small_variants/noalt/{batch}-somatic.vcf.gz.tbi',
    group: "somatic_anno"
    shell:
        'bcftools view -R {input.noalts_bed} {input.vcf} -Oz -o {output.vcf} '
        '&& tabix -p vcf {output.vcf}'

rule somatic_vcf_sage1:
    input:
        vcf = rules.somatic_vcf_select_noalt.output.vcf,
        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumors[0].bam,
        normal_bam = lambda wc: batch_by_name[wc.batch].normals[0].bam,
    output:
        vcf = 'work/{batch}/small_variants/sage1/{batch}-somatic.vcf.gz',
        tbi = 'work/{batch}/small_variants/sage1/{batch}-somatic.vcf.gz.tbi',
        sage_vcf = '{batch}/small_variants/sage1/{batch}-sage.vcf.gz',
        sage_tbi = '{batch}/small_variants/sage1/{batch}-sage.vcf.gz.tbi',
    params:
        genome = run.genome_build,
        genomes_dir = refdata.genomes_dir,
        work_dir = 'work/{batch}/small_variants/sage1',
        unlock_opt = ' --unlock' if config.get('unlock', 'no') == 'yes' else '',
        tumor_sample = lambda wc: batch_by_name[wc.batch].tumors[0].rgid,
        normal_sample = lambda wc: batch_by_name[wc.batch].normals[0].rgid,
    resources:
        mem_mb = 20000
    group: "somatic_anno"
    shell:
        'sage '
        '-t {input.tumor_bam} '
        '-n {input.normal_bam} '
        '-v {input.vcf} '
        '-o {output.vcf} '
        '-s {output.sage_vcf} '
        '-tn {params.tumor_sample} '
        '-nn {params.normal_sample} '
        '-w {params.work_dir} '
        '-g {params.genome} '
        '--genomes-dir {params.genomes_dir} '
        '{params.unlock_opt}'

rule somatic_vcf_annotate:
    input:
        vcf = lambda wc: f'work/{wc.batch}/small_variants/sage1/{wc.batch}-somatic.vcf.gz' \
                if batch_by_name[wc.batch].somatic_vcf else
                         f'work/{wc.batch}/small_variants/noalt/{wc.batch}-somatic.vcf.gz',
    output:
        vcf = 'work/{batch}/small_variants/annotate/{batch}-somatic.vcf.gz',
        subset_highly_mutated_stats =
          'work/{batch}/small_variants/somatic_anno/subset_highly_mutated_stats.yaml',
    params:
        genome = run.genome_build,
        genomes_dir = refdata.genomes_dir,
        work_dir = 'work/{batch}/small_variants',
        unlock_opt = ' --unlock' if config.get('unlock', 'no') == 'yes' else '',
        tumor_names   = lambda wc: ','.join(s.rgid for s in batch_by_name[wc.batch].tumors),
        normal_names  = lambda wc: ','.join(s.rgid for s in batch_by_name[wc.batch].normals),
    resources:
        mem_mb = 20000
    group: "somatic_anno"
    run:
        t_samples = [s.rgid for s in batch_by_name[wildcards.batch].tumors]
        assert t_samples
        n_samples = [s.rgid for s in batch_by_name[wildcards.batch].normals]
        tn_opt = f'-tn {",".join(t_samples)}'
        nn_opt = f'-nn {",".join(n_samples)}' if n_samples else ''
        rn_opt = f'-rn {wildcards.batch}_rna' if batch_by_name[wildcards.batch].rna_samples else ''
        shell(
            f'anno_somatic_vcf {input.vcf} -o {output.vcf} '
            f'{tn_opt} {nn_opt} {rn_opt} '
            f'-w {params.work_dir} -g {params.genome} --genomes-dir {params.genomes_dir} '
            f'{params.unlock_opt}'
        )

rule somatic_vcf_filter:
    input:
        vcf = rules.somatic_vcf_annotate.output.vcf,
    output:
        vcf = '{batch}/small_variants/{batch}-somatic.vcf.gz',
    group: "somatic_filt"
    params:
        min_af = config.get('min_af', 0.1),
    run:
        t_samples = [s.rgid for s in batch_by_name[wildcards.batch].tumors]
        assert t_samples
        n_samples = [s.rgid for s in batch_by_name[wildcards.batch].normals]
        tn_opt = f'-tn {",".join(t_samples)}'
        nn_opt = f'-nn {",".join(n_samples)}' if n_samples else ''
        rn_opt = f'-rn {wildcards.batch}_rna' if batch_by_name[wildcards.batch].rna_samples else ''
        shell(
            f'filter_somatic_vcf {input.vcf} -o {output.vcf} '
            f'{tn_opt} {nn_opt} {rn_opt} '
            f'--min-af {params.min_af}'
        )

rule somatic_vcf_filter_pass:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-PASS.vcf.gz',
        tbi = '{batch}/small_variants/{batch}-somatic-PASS.vcf.gz.tbi',
    group: "somatic_filt"
    shell:
        'bcftools view -f.,PASS {input.vcf} -Oz -o {output.vcf} '
        '&& tabix -f -p vcf {output.vcf}'

rule bcftools_stats_somatic:
    input:
        rules.somatic_vcf_filter_pass.output.vcf
    output:
        '{batch}/small_variants/stats/{batch}_bcftools_stats.txt'
    group: "somatic_filt"
    params:
        sname = lambda wc: batch_by_name[wc.batch].tumors[0].rgid,
    shell:
        'bcftools stats -s {params.sname} {input} | sed s#{input}#{params.sname}# > {output}'

rule somatic_stats_report:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf,
        full_vcf = rules.somatic_vcf_select_noalt.output.vcf,
        subset_highly_mutated_stats = 'work/{batch}/small_variants/somatic_anno/subset_highly_mutated_stats.yaml',
    output:
        'work/{batch}/small_variants/{batch}_somatic_stats.yml',
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumors[0].name,
    group: "somatic_filt"
    run:
        snps, indels, others = cnt_vars(input.vcf, passed=True)
        all_snps, all_indels, all_others = cnt_vars(input.full_vcf)
        subset_genes = None

        with open(input.subset_highly_mutated_stats) as inp:
            stats = yaml.safe_load(inp)
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
                'id': 'umccrise',
                'data': {
                    params.sample: data
                }
            }
            yaml.dump(data, out, default_flow_style=False)

rule somatic_vcf2maf:
    input:
        vcf = rules.somatic_vcf_filter_pass.output.vcf,
        fa = refdata.get_ref_file(genome=run.genome_build, key='fa'),
    output:
        maf = '{batch}/small_variants/{batch}-somatic-PASS.maf',
    params:
        # Use with caution as vcf2maf doesn't support multi-tumor and multi-normal VCFs,
        # so fields parsed from FORMAT will be for the first tumor and normal only
        tname = lambda wc: [s.rgid for s in batch_by_name[wc.batch].tumors][0],
        nname = lambda wc: [s.rgid for s in batch_by_name[wc.batch].normals][0],
        ncbi_build = {'hg38': 'GRCh38', 'GRCh37': 'GRCh37'}.get(run.genome_build),
        uncompressed_tmp_vcf = 'work/{batch}/small_variants/{batch}-somatic.vcf.tmp',
    shell:
        'gunzip -c {input.vcf} > {params.uncompressed_tmp_vcf} '
        '&& vcf2maf.pl --inhibit-vep ' 
        '--input-vcf {params.uncompressed_tmp_vcf} '
        '--output-maf {output.maf} '
        '--ref-fasta {input.fa} '
        '--filter-vcf 0 '
        '--tumor-id {params.tname} '
        '--normal-id {params.nname} '
        '--ncbi-build {params.ncbi_build} 2> >(grep -v "Use of uninitialized value" >&2)'
        '&& rm {params.uncompressed_tmp_vcf}'


rule pierian:
    input:
        snv = rules.somatic_vcf_filter_pass.output.vcf,
        sv = '{batch}/structural/{batch}-manta.vcf.gz',
        svtbi = '{batch}/structural/{batch}-manta.vcf.gz.tbi',
        cnv = '{batch}/purple/{batch}.purple.cnv.somatic.tsv',
    output:
        snv_renamed = '{batch}/pierian/{batch}.somatic-PASS-single.grch38.vcf.gz',
        sv_renamed = '{batch}/pierian/{batch}-manta.single.vcf.gz',
        cnv_renamed = '{batch}/pierian/{batch}.purple.cnv',
    run:
        # SV: copy sv since no processing required
        shutil.copy(f'{input.sv}', f'{output.sv_renamed}')
        shutil.copy(f'{input.svtbi}', f'{output.sv_renamed}.tbi')

        # CNV: handle PURPLE renamed column
        shell("sed 's/AlleleCopyNumber/AllelePloidy/g' {input.cnv} > {output.cnv_renamed}")

        # SNV: grab tumor column from snv vcf
        t_name = batch_by_name[wildcards.batch].tumors[0].rgid
        vcf_samples = cyvcf2.VCF(input.snv).samples
        assert t_name in vcf_samples, f"Tumor name {t_name} not in VCF {input.snv}, available: {vcf_samples}"
        shell("bcftools view {input.snv} -s {t_name} -Oz -o {output.snv_renamed} && tabix -p vcf {output.snv_renamed}")


#############

rule somatic:
    input:
        somatic_vcfs = expand(rules.somatic_vcf_filter.output.vcf, batch=batch_by_name.keys()),
        somatic_mafs = expand(rules.somatic_vcf2maf.output.maf, batch=batch_by_name.keys()),
        somatic_pass_vcfs = expand(rules.somatic_vcf_filter_pass.output.vcf, batch=batch_by_name.keys()),
        pierian = expand(rules.pierian.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/somatic.done'))


rule maf:
    input:
        somatic_mafs = expand(rules.somatic_vcf2maf.output.maf, batch=batch_by_name.keys()),
    output:
        temp(touch('log/maf.done'))



