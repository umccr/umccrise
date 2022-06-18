"""
Structural variants
"""

import glob
from os.path import join, dirname, basename, isfile
from cyvcf2 import VCF
from ngs_utils.file_utils import safe_mkdir, verify_dir, get_ungz_gz, verify_file
from ngs_utils.logger import critical
from ngs_utils.vcf_utils import count_vars, vcf_contains_field, iter_vcf
from reference_data import api as refdata


vcftobedpe = 'vcfToBedpe'

MAX_SVS = 50000


localrules: structural, structural_batch, copy_sv_vcf_ffpe_mode


# Keep PASSed variants from final Manta, and reset (some) INFO & FILTER annotations
rule sv_keep_pass:
    input:
        vcf = lambda wc: batch_by_name[wc.batch].sv_vcf,
    output:
        vcf = 'work/{batch}/structural/keep_pass/{batch}-manta.vcf.gz',
        tbi = 'work/{batch}/structural/keep_pass/{batch}-manta.vcf.gz.tbi',
    group: "sv_vcf"
    run:
        cmd = f'cat {input.vcf}'
        # remove previous annotation
        filts_to_remove = [f'{f}' for f in ['INFO/SIMPLE_ANN', 'INFO/SV_HIGHEST_TIER',
                                            'FILTER/Intergenic', 'FILTER/MissingAnn', 'FILTER/REJECT']
                           if vcf_contains_field(input.vcf, f)]
        if filts_to_remove:
            cmd += f' | bcftools annotate -x "' + ','.join(f'{f}' for f in filts_to_remove) + '"'
        cmd += ' | bcftools view -f.,PASS -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'
        shell(cmd)

# Run snpEff if no INFO/ANN
rule sv_snpeff_maybe:
    input:
        vcf = rules.sv_keep_pass.output.vcf,
        tbi = rules.sv_keep_pass.output.tbi,
    output:
        vcf  = 'work/{batch}/structural/snpeff/{batch}-sv-snpeff.vcf.gz',
        tbi  = 'work/{batch}/structural/snpeff/{batch}-sv-snpeff.vcf.gz.tbi',
    params:
        genome = run.genome_build,
        tmp_dir = 'work/{batch}/structural/snpeff/tmp',
        csv  = 'work/{batch}/structural/snpeff/{batch}-sv-snpeff-stats.csv',
        html = 'work/{batch}/structural/snpeff/{batch}-sv-snpeff-stats.html',
    resources:
        mem_mb = min(8000, max(30000, 3000*threads_per_batch + 1000))
    group: "sv_vcf",
    threads:
        threads_per_batch
    run:
        if vcf_contains_field(input.vcf, 'INFO/ANN'):
            logger.info(f'Manta {input.vcf} is already annotated with SnpEff, reusing')
            shell(f"cp {input.vcf} {output.vcf}")
            shell(f"cp {input.tbi} {output.tbi}")
        else:
            snpeff_db = refdata.get_ref_file(genome=params.genome, key='snpeff')
            snpeff_db_dir = dirname(snpeff_db)
            snpeff_db_name = 'GRCh38.86'  # it takes the genome build from `envs/umccrise/share/snpeff-4.3.1t-3/snpEff.config`
                                          # so it doesn't matter if the subdir is named GRCh38.92
            mem_jvm_gb = min(7, max(30, 3*threads))
            jvm_opts = f'-Xms750m -Xmx{mem_jvm_gb}g'
            java_args = f'-Djava.io.tmpdir={params.tmp_dir}'

            shell('snpEff {jvm_opts} {java_args} '
                  '-dataDir {snpeff_db_dir} {snpeff_db_name} '
                  '-hgvs -cancer -i vcf -o vcf '
                  '-csvStats {params.csv} -s {params.html} '
                  '{input.vcf} '
                  '| bgzip --threads {threads} -c > {output.vcf} && tabix -p vcf {output.vcf}')
            verify_file(output.vcf, is_critical=True)

# Run VEP on final Manta
rule sv_vep:
    input:
        vcf = lambda wc: batch_by_name[wc.batch].sv_vcf,
    output:
        vcf = 'work/{batch}/structural/vep/{batch}-sv-vep.vcf.gz',
        tbi = 'work/{batch}/structural/vep/{batch}-sv-vep.vcf.gz.tbi',
    params:
        genome = run.genome_build,
        pcgr_dir = refdata.get_ref_file(genome=run.genome_build, key='pcgr_data'),
    group: "sv_vcf",
    threads:
        threads_per_batch
    run:
        vep_genome = 'GRCh38' if ('38' in params.genome) else 'GRCh37'
        pcgr_genome = vep_genome.lower()
        vep_dir = join(params.pcgr_dir, pcgr_genome, '.vep')
        verify_dir(vep_dir, is_critical=True)
        ref_glob = join(vep_dir, 'homo_sapiens', f'*_{vep_genome}', f'Homo_sapiens.{vep_genome}.dna.primary_assembly.fa.gz')
        ref_found = glob.glob(ref_glob)
        if len(ref_found) == 0:
            critical(f'Can\'t find VEP assembly fasta file at {ref_glob}')
        ref_fa = ref_found[0]
        vep_version = basename(dirname(ref_fa)).split('_')[0]  # homo_sapiens/{vep_version}_{vep_genome}/Homo_sapiens.{vep_genome}.dna.primary_assembly.fa.gz

        # see https://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html
        opts = [
            '--af_gnomad', '--allele_number', '--appris', f'--assembly {vep_genome}',
            '--biotype', '--cache', f'--cache_version {vep_version}', '--canonical',
            '--ccds', '--check_ref', f'--dir {vep_dir}', '--dont_skip',
            '--flag_pick_allele_gene', f'--fasta {ref_fa}', '--force_overwrite',
            f'--fork {threads}', '--format vcf', '--gencode_basic',
            '--hgvs', '--no_escape', '--numbers', '--offline',
            '--pick_order rank,appris,biotype,tsl,ccds,canonical,length,mane',
            '--species homo_sapiens', '--symbol', '--total_length',
            '--variant_class', '--vcf', '--xref_refseq',
        ]
        opts_line = ' '.join(opts)
        out_ungz, out_gz = get_ungz_gz(output.vcf)
        # use VEP through PCGR's conda env
        shell(conda_cmd.format('pcgr') +
            'vep --input_file {input.vcf} --output_file {out_ungz} {opts_line}')
        shell('bgzip {out_ungz} && tabix {out_gz}')

# Handle SnpEff capitalising ALT - see https://github.com/pcingola/SnpEff/issues/237
# BPI and bedtools>=2.29.2 will crash if left as is.
rule fix_snpeff:
    input:
        vcf = rules.sv_snpeff_maybe.output.vcf,
    output:
        vcf = 'work/{batch}/structural/snpeff/{batch}-sv-snpeff-fix.vcf.gz',
    group: "sv_vcf"
    shell:
        # TODO: swap my monstrosity with 'sed -e pat1 -e pat2...'
        'gunzip -c {input.vcf} '
        '| sed "s/CHR/chr/" '
        '| sed "s/chrOM/CHROM/" '
        '| sed "s/V1_RANDOM/v1_random/" '
        '| sed "s/V2_RANDOM/v2_random/" '
        '| sed "s/V1_ALT/v1_alt/" '
        '| sed "s/V2_ALT/v2_alt/" '
        '| sed "s/V2_DECOY/v1_decoy/" '
        '| sed "s/V2_DECOY/v2_decoy/" '
        '| sed "s/V1/v1/" '
        '| sed "s/V2/v2/" '
        '| bgzip -c > {output.vcf} '
        '&& tabix -p vcf {output.vcf}'

rule sv_prioritize:
    input:
        vcf = rules.fix_snpeff.output.vcf
    output:
        vcf = 'work/{batch}/structural/prioritize/{batch}-sv-eff-prio.vcf.gz',
        tbi = 'work/{batch}/structural/prioritize/{batch}-sv-eff-prio.vcf.gz.tbi',
    params:
        genome = run.genome_build
    group: "sv_vcf"
    run:
        shell(f'cat {input.vcf} | prioritize_sv -g {params.genome} | bgzip -c > {output.vcf} && tabix -p vcf {output.vcf}')
        before = count_vars(input.vcf)
        after = count_vars(output.vcf)
        assert before == after, (before, after)

# If there are too many SV calls (FFPE?), also remove unprioritized SVs
# Also removing older very cluttered ANN field
rule sv_subsample_if_too_many:
    input:
        vcf = rules.sv_prioritize.output.vcf
    output:
        vcf = 'work/{batch}/structural/sv_subsample_if_too_many/{batch}-manta.vcf'
    group: "sv_vcf"
    run:
        if count_vars(input.vcf) < MAX_SVS:
            shell(f'bcftools view {input.vcf} -o {output.vcf}')
        else:
            if count_vars(input.vcf, bcftools_filter_expr='-i "SV_TOP_TIER < 4"') < MAX_SVS:
                cmd = f'bcftools filter -i "SV_TOP_TIER < 4" {input.vcf}'
            elif count_vars(input.vcf, bcftools_filter_expr='-i "SV_TOP_TIER < 3"') < MAX_SVS:
                cmd = f'bcftools filter -i "SV_TOP_TIER < 3" {input.vcf}'
            else:
                cmd = f'bcftools filter -i "SV_TOP_TIER < 2" {input.vcf}'
            cmd += f' -o {output.vcf}'
            shell(cmd)

# if BPI was disabled in bcbio
rule sv_bpi_maybe:
    input:
        vcf = rules.sv_subsample_if_too_many.output.vcf,
        tumor_bam = lambda wc: batch_by_name[wc.batch].tumors[0].bam,
        tumor_bai = lambda wc: batch_by_name[wc.batch].tumors[0].bam + '.bai',
        normal_bam = lambda wc: batch_by_name[wc.batch].normals[0].bam,
        normal_bai = lambda wc: batch_by_name[wc.batch].normals[0].bam + '.bai',
    output:
        vcf = 'work/{batch}/structural/maybe_bpi/{batch}-manta.vcf'
    group: "sv_vcf"
    log:
        'log/structural/{batch}/{batch}-bpi_stats.txt'
    params:
        xms = 1000,
        xmx = 30000,
        tmp_dir = '{batch}/structural/maybe_bpi/tmp_dir'
    resources:
        mem_mb = 30000
    run:
        if vcf_contains_field(input.vcf, 'BPI_AF', 'INFO'):  # already BPI'ed
            shell('cp {input.vcf} {output.vcf}')
        else:
            # if not is_ffpe:  # running BPI only for non-FFPE samples
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
            # else:  # fake BPI_AF from the original AF
            #     def func(rec, vcf):
            #         rec.INFO['BPI_AF'] = rec.INFO['AF']
            #         return rec
            #     iter_vcf(input.vcf, output.vcf, func)

# Keep all with read support above 10x; or allele frequency above 10%, but only if read support is above 5x
rule filter_sv_vcf:
    input:
        vcf = rules.sv_bpi_maybe.output.vcf
    output:
        vcf = 'work/{batch}/structural/filt/{batch}-manta.vcf.gz',
        tbi = 'work/{batch}/structural/filt/{batch}-manta.vcf.gz.tbi',
    group: "sv_vcf"
    run:
        t_name = batch_by_name[wildcards.batch].tumors[0].rgid
        vcf_samples = VCF(input.vcf).samples
        assert t_name in vcf_samples, f"Tumor name {t_name} not in VCF {input.vcf}, available: {vcf_samples}"
        tumor_id = vcf_samples.index(t_name)
        # tumor_id = VCF(input.vcf).samples.index(batch_by_name[wildcards.batch].tumor.name)
        print(f'Derived tumor VCF index: {tumor_id}')
        shell('''
bcftools view -f.,PASS {input.vcf} |
bcftools filter -e "SVTYPE == 'BND' & FORMAT/SR[{tumor_id}:1] - FORMAT/PR[{tumor_id}:1] > 0" |
bcftools filter -e "SV_TOP_TIER > 2 & FORMAT/SR[{tumor_id}:1]<5  & FORMAT/PR[{tumor_id}:1]<5" |
bcftools filter -e "SV_TOP_TIER > 2 & FORMAT/SR[{tumor_id}:1]<10 & FORMAT/PR[{tumor_id}:1]<10 & (BPI_AF[0] < 0.1 | BPI_AF[1] < 0.1)" |
bcftools view -s {t_name} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}

''')

rule reprioritize_rescued_svs:
    input:
        dummy = lambda wc: 'work/{batch}/purple/{batch}.purple.purity.tsv' \
            if ('purple' in stages or isfile('work/{batch}/purple/{batch}.purple.sv.vcf.gz')) \
            else [],
    output:
        vcf = 'work/{batch}/structural/sv_after_purple/{batch}-manta.vcf.gz',
        tbi = 'work/{batch}/structural/sv_after_purple/{batch}-manta.vcf.gz.tbi',
    params:
        vcf = 'work/{batch}/purple/{batch}.purple.sv.vcf.gz',
        genome = run.genome_build
    group: "sv_after_purple"
    run:
        cmd = f'gunzip -c {params.vcf}' \
            f' | prioritize_sv -g {params.genome}' \
            f' | bcftools annotate -x "INFO/ANN"' \
            f' -Oz -o {output.vcf}' \
            f' && tabix -p vcf {output.vcf}'
        shell(cmd)
        before = count_vars(params.vcf)
        after = count_vars(output.vcf)
        assert before == after, (before, after)


rule copy_sv_vcf_ffpe_mode:
    input:
        vcf = lambda wc: rules.filter_sv_vcf.output.vcf \
            if ((batch_by_name[wc.batch].somatic_caller == 'strelka2') or ('purple' not in stages)) \
            else f'work/{wc.batch}/structural/sv_after_purple/{wc.batch}-manta.vcf.gz'
    output:
        vcf = '{batch}/structural/{batch}-manta.vcf.gz',
        tbi = '{batch}/structural/{batch}-manta.vcf.gz.tbi',
    shell:
        'cp {input.vcf} {output.vcf} ; '
        'cp {input.vcf}.tbi {output.tbi}'


def parse_info_field(rec, name):
    val = rec.INFO.get(name)
    if val is None:
        return ''
    elif isinstance(val, float) or isinstance(val, int) or isinstance(val, bool) or isinstance(val, str):
        return str(val)
    else:
        return ','.join(map(str, val))


# Produce a TSV file for further analysis in Rmd
# caller  sample                chrom   start       end         svtype  lof  annotation                                                                 split_read_support  paired_support_PE  paired_support_PR  somaticscore  tier
# manta   PRJ180253_E190-T01-D  1       161513440   161595209   DUP          DUP|GENE_FUSION|FCGR2B&RP11-25K21.6|ENSG00000273112|NOT_PRIORITISED|3,...                                         67,8
# manta   PRJ180253_E190-T01-D  8       33320739    33321344    DEL          DEL|UPSTREAM_GENE_VARIANT|FUT10|ENST00000518672|NOT_PRIORITISED|3          76,7                                   33,1
# manta   PRJ180253_E190-T01-D  11      118802640   118803304   DEL          DEL|DOWNSTREAM_GENE_VARIANT|RN7SL688P|ENST00000471754|NOT_PRIORITISED|3    61,8                                   07,2
rule prep_sv_tsv:
    input:
        vcf = '{batch}/structural/{batch}-manta.vcf.gz',
    output:
        '{batch}/structural/{batch}-manta.tsv'
    group: "sv_after_purple"
    run:
        sample_name = batch_by_name[wildcards.batch].tumors[0].name
        rgid = batch_by_name[wildcards.batch].tumors[0].rgid
        vcf_samples = VCF(input.vcf).samples
        assert rgid in vcf_samples, f"Tumor sample {rgid} is not in VCF {input.vcf}, available: {vcf_samples}"
        tumor_id = VCF(input.vcf).samples.index(rgid)
        with open(output[0], 'w') as out:
            header = ["caller", "sample", "chrom", "start", "end", "svtype",
                      "split_read_support", "paired_support_PE", "paired_support_PR", "AF_BPI", "somaticscore",
                      "tier", "annotation",
                      "AF_PURPLE", "CN_PURPLE", "CN_change_PURPLE", "Ploidy_PURPLE", "PURPLE_status", "START_BPI", "END_BPI", "ID", "MATEID", "ALT"]
            out.write('\t'.join(header) + '\n')
            for rec in VCF(input.vcf):
                tier = parse_info_field(rec, 'SV_TOP_TIER')
                if not tier:
                    tier = '4'
                simple_ann = parse_info_field(rec, 'SIMPLE_ANN')

                PURPLE_status = ''
                if rec.FILTER and 'INFERRED' in rec.FILTER:
                    PURPLE_status = 'INFERRED'
                    if not simple_ann:
                        simple_ann = f'{rec.INFO["SVTYPE"]}||||From_CNV|{tier}'
                elif rec.INFO.get('RECOVERED'):
                    PURPLE_status = 'RECOVERED'

                data = ['manta',
                        sample_name,
                        rec.CHROM,
                        rec.POS,
                        rec.INFO.get('END', ''),
                        rec.INFO['SVTYPE'],
                        ','.join(map(str, rec.format('SR')[tumor_id])) if 'SR' in rec.FORMAT else '',
                        ','.join(map(str, rec.format('PE')[tumor_id])) if 'PE' in rec.FORMAT else '',
                        ','.join(map(str, rec.format('PR')[tumor_id])) if 'PR' in rec.FORMAT else '',
                        parse_info_field(rec, 'BPI_AF'),
                        parse_info_field(rec, 'SOMATICSCORE'),
                        tier,
                        simple_ann,
                        parse_info_field(rec, 'PURPLE_AF'),
                        parse_info_field(rec, 'PURPLE_CN'),
                        parse_info_field(rec, 'PURPLE_CN_CHANGE'),
                        parse_info_field(rec, 'PURPLE_PLOIDY'),
                        PURPLE_status,
                        parse_info_field(rec, 'BPI_START'),
                        parse_info_field(rec, 'BPI_END'),
                        rec.ID,
                        parse_info_field(rec, 'MATEID'),
                        rec.ALT[0]
                        ]
                out.write('\t'.join(map(str, data)) + '\n')


# At least for the most conservative manta calls, generate a file for viewing in Ribbon
rule ribbon_filter_manta:
    input:
        manta_vcf = '{batch}/structural/{batch}-manta.vcf.gz',
    output:
        'work/{batch}/structural/ribbon/manta.vcf'
    group: "sv_after_purple"
    shell:
        'bcftools view {input.manta_vcf} > {output}'

rule ribbon_filter_vcfbedtope_starts:
    input:
        bed = rules.ribbon_filter_manta.output[0],
        fai = refdata.get_ref_file(run.genome_build, key='fa') + '.fai'
    output:
        'work/{batch}/structural/ribbon/manta-starts.bed'
    params:
        vcftobedpe = vcftobedpe
    group: "sv_after_purple"
    shell:
        'cat {input.bed} | {params.vcftobedpe}'
        ' | cut -f 1-3'
        ' | bedtools slop -b 5000 -i stdin -g {input.fai}'
        ' > {output}'

rule ribbon_filter_vcfbedtope_ends:
    input:
        bed = rules.ribbon_filter_manta.output[0],
        fai = refdata.get_ref_file(run.genome_build, key='fa') + '.fai'
    output:
        'work/{batch}/structural/ribbon/manta-ends.bed'
    params:
        vcftobedpe = vcftobedpe
    group: "sv_after_purple"
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
    group: "sv_after_purple"
    shell:
        'cat {input.starts} {input.ends} | bedtools sort -i stdin | bedtools merge -i stdin > {output}'


#### Convert matna VCF to bedpe ####
rule bedpe:
    input:
        manta_vcf = '{batch}/structural/{batch}-manta.vcf.gz',
    output:
        '{batch}/structural/{batch}-manta.bedpe'
    params:
        vcftobedpe = vcftobedpe
    group: "sv_after_purple"
    shell:
        'bcftools view {input.manta_vcf}'
        ' | {params.vcftobedpe}'
        ' | cut -f 1-7'
        ' > {output}'


#############

rule structural_batch:
    input:
        lambda wc: rules.bedpe.output       if batch_by_name[wc.batch].sv_vcf else [],
        lambda wc: rules.ribbon.output      if batch_by_name[wc.batch].sv_vcf else [],
        lambda wc: rules.prep_sv_tsv.output if batch_by_name[wc.batch].sv_vcf else [],
    output:
        temp(touch('log/structural_{batch}.done'))


rule structural:
    input:
        expand(rules.structural_batch.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/structural.done'))
