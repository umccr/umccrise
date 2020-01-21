"""
Structural variants
------------------
Re-do the CNV plots. This will need lots of love (drop gene names, make the scatterplot viable again, etc.).
"""
import glob
from os.path import join, dirname, basename
from cyvcf2 import VCF
from ngs_utils.file_utils import safe_mkdir, verify_dir, get_ungz_gz
from ngs_utils.logger import critical
from vcf_stuff import count_vars, vcf_contains_field, iter_vcf


vcftobedpe = 'vcfToBedpe'


localrules: structural


rule sv_vep:
    input:
        vcf = lambda wc: batch_by_name[wc.batch].sv_vcf,
    output:
        vcf = 'work/{batch}/structural/vep/{batch}-sv-vep.vcf.gz',
        tbi = 'work/{batch}/structural/vep/{batch}-sv-vep.vcf.gz.tbi',
    params:
        genome = run.genome_build,
        pcgr_dir = hpc.get_ref_file(key='pcgr_data'),
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

        opts = [
            '--hgvs',               # Add HGVS nomenclature based on Ensembl stable identifiers to the output.
                                    # Both coding and protein sequence names are added where appropriate.
                                    # HGVS notations given on Ensembl identifiers are versioned.
            '--af_gnomad',          # Include allele frequency from Genome Aggregation Database (gnomAD)
                                    # exome populations. Note only data from the gnomAD exomes are included; to
                                    # retrieve data from the additional genomes data set, see this guide.
                                    # Must be used with --cache
            '--variant_class',      # Output the Sequence Ontology
                                    # [variant class](https://asia.ensembl.org/info/genome/variation/prediction/classification.html#classes)
            '--symbol',             # Adds the gene symbol (e.g. HGNC) (where available) to the output
            '--ccds',               # Adds the CCDS transcript identifer (where available) to the output
            '--appris',             # Adds the [APPRIS](https://asia.ensembl.org/Help/Glossary?id=521)
                                    # isoform annotation for this transcript to the output
            '--biotype',            # Adds the biotype of the transcript or regulatory feature
            '--canonical',          # Adds a flag indicating if the transcript is the canonical transcript for the gene
            '--gencode_basic',      # Limit your analysis to transcripts belonging to the GENCODE basic set.
                                    # This set has fragmented or problematic transcripts removed
            '--cache',              # Enables use of the cache. Add --refseq or --merged to use the refseq or merged cache
            '--numbers',            # Adds affected exon and intron numbering to to output. Format is Number/Total.
            '--total_length',       # Give cDNA, CDS and protein positions as Position/Length
            '--allele_number',      # dentify allele number from VCF input, where 1 = first ALT allele,
                                    # 2 = second ALT allele etc. Useful when using --minimal
            '--no_escape',          # Don't URI escape HGVS strings
            '--xref_refseq',        # Output aligned RefSeq mRNA identifier for transcript
            '--vcf',                # Writes output in VCF format. Consequences are added in the INFO field of
                                    # the VCF file, using the key "CSQ". Data fields are encoded separated by "|";
                                    # the order of fields is written in the VCF header. Output fields in the "CSQ"
                                    # INFO field can be selected by using --fields.
            '--check_ref',          # Force VEP to check the supplied reference allele against the sequence stored in
                                    # the Ensembl Core database or supplied FASTA file. Lines that do not match are
                                    # skipped
            '--dont_skip',          # Don't skip input variants that fail validation, e.g. those that fall on
                                    # unrecognised sequences. Combining --check_ref with --dont_skip will add a
                                    # CHECK_REF output field when the given reference does not match the underlying
                                    # reference sequence.
            '--flag_pick_allele_gene',  # As per --pick_allele_gene, but adds the PICK flag to the chosen block of
                                        # consequence data and retains others
            '--pick_order rank,appris,biotype,tsl,ccds,canonical,length,mane',
            '--force_overwrite',
            '--species homo_sapiens',
            f'--assembly {vep_genome}',
            '--offline',
            f'--fork {threads}',
            f'--dir {vep_dir}',
            f'--cache_version {vep_version}',
            f'--fasta {ref_fa}',
            '--format vcf'          # By default, VEP auto-detects the input file format. However it fails
                                    # to autodetect with SVs in VCF. So specifying explicitly as VCF
        ]
        opts_line = ' '.join(opts)
        out_ungz, out_gz = get_ungz_gz(output.vcf)
        shell(conda_cmd.format('pcgr') +
            'vep --input_file {input.vcf} --output_file {out_ungz} {opts_line}')
        shell('bgzip {out_ungz} && tabix {out_gz}')

rule sv_prioritize:
    input:
        vcf = 'work/{batch}/structural/vep/{batch}-sv-vep.vcf.gz',
    output:
        vcf = 'work/{batch}/structural/prioritize/{batch}-sv-vep-prio.vcf.gz',
        tbi = 'work/{batch}/structural/prioritize/{batch}-sv-vep-prio.vcf.gz.tbi',
    params:
        genome = run.genome_build
    group: "sv_vcf"
    run:
        assert vcf_contains_field(input.vcf, 'INFO/ANN'), f'Manta {input.vcf} must be annotated with SnpEff'

        cmd = f'cat {input.vcf}'
        # remove previous annotation
        filts_to_remove = [f'{f}' for f in ['INFO/SIMPLE_ANN', 'INFO/SV_HIGHEST_TIER',
                                            'FILTER/Intergenic', 'FILTER/MissingAnn', 'FILTER/REJECT']
                           if vcf_contains_field(input.vcf, f)]
        if filts_to_remove:
            cmd += f' | bcftools annotate -x "' + ','.join(f'{f}' for f in filts_to_remove) + '"'

        cmd += f' | prioritize_sv -g {params.genome} | bgzip -c > {output.vcf} && tabix -p vcf {output.vcf}'
        # keeping INFO/ANN here because later (after PURPLE we reprioritizing)
        shell(cmd)
        before = count_vars(input.vcf)
        after = count_vars(output.vcf)
        assert before == after, (before, after)

# Keep passed variants.
# If there are too many SV calls (FFPE?), also remove unprioritized SVs
# Also removing older very cluttered ANN field
rule sv_keep_pass:
    input:
        vcf = rules.sv_prioritize.output.vcf
    output:
        vcf = 'work/{batch}/structural/keep_pass/{batch}-manta.vcf'
    group: "sv_vcf"
    run:
        cmd = f'bcftools view -f.,PASS {input.vcf}'
        if count_vars(input.vcf, filter='.,PASS') > 100000:
            cmd += ' | bcftools filter -i "SV_TOP_TIER <= 3"'
        cmd += f' -o {output.vcf}'
        shell(cmd)

# rule sv_maybe_keep_prioritize:
#     input:
#         vcf = rules.sv_keep_pass.output.vcf
#     output:
#         vcf = 'work/{batch}/structural/maybe_keep_prio/{batch}-manta.vcf'
#     group: "sv_vcf"
#     run:
#         if count_vars(input.vcf, filter='.,PASS') > 1000:
#             # Still too cluttered (FFPE?) - removing all tier=4 as well
#             def func(rec):
#                 ann = rec.INFO.get('SIMPLE_ANN')
#                 if ann:
#                     anns = [a for a in ann.split(',') if '|NOT_PRIORITISED|' not in a]
#                     if anns:
#                         ann = ','.join(anns)
#                         rec.INFO['SIMPLE_ANN'] = ann
#                         return rec
#             iter_vcf(input.vcf, output.vcf, func)
#         else:
#             shell('cp {input.vcf} {output.vcf}')

# if BPI was disabled in bcbio
rule sv_maybe_bpi:
    input:
        vcf = rules.sv_keep_pass.output.vcf,
        tumor_bam = lambda wc: batch_by_name[wc.batch].tumor.bam,
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
    output:
        vcf = 'work/{batch}/structural/maybe_bpi/{batch}-manta.vcf'
    group: "sv_vcf"
    log:
        'log/structural/{batch}/{batch}-bpi_stats.txt'
    params:
        xms = 1000,
        xmx = 16000,
        tmp_dir = '{batch}/structural/maybe_bpi/tmp_dir'
    resources:
        mem_mb = 16000
    run:
        # Handle SnpEff capitalising ALT (see https://github.com/pcingola/SnpEff/issues/237).
        # BPI and bedtools>=2.29.2 will crash if left as is.
        shell('sed -i "s/CHR/chr/" {input.vcf} && sed -i "s/chrOM/CHROM/" {input.vcf}; ')
        if vcf_contains_field(input.vcf, 'BPI_AF', 'INFO'):  # already BPI'ed
            shell('cp {input.vcf} {output.vcf}')
        else:
            if not is_ffpe:  # running BPI only for non-FFPE samples
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
            else:  # fake BPI_AF from the original AF
                def func(rec):
                    rec.INFO['BPI_AF'] = rec.INFO['AF']
                    return rec
                iter_vcf(input.vcf, output.vcf, func)

# Keep all with read support above 10x; or allele frequency above 10%, but only if read support is above 5x
rule filter_sv_vcf:
    input:
        vcf = rules.sv_maybe_bpi.output.vcf
    output:
        vcf = 'work/{batch}/structural/filt/{batch}-manta.vcf'
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name
    group: "sv_vcf"
    run:
        print(f'VCF samples: {VCF(input.vcf).samples}')
        print(f'Bcbio batch tumor name: {batch_by_name[wildcards.batch].tumor.name}')
        tumor_id = VCF(input.vcf).samples.index(params.sample)
        # tumor_id = VCF(input.vcf).samples.index(batch_by_name[wildcards.batch].tumor.name)
        print(f'Derived tumor VCF index: {tumor_id}')
        shell('''
bcftools view -f.,PASS {input.vcf} |
bcftools filter -e "SVTYPE == 'BND' & FORMAT/SR[{tumor_id}:1] - FORMAT/PR[{tumor_id}:1] > 0" |
bcftools filter -e "SV_TOP_TIER > 2 & FORMAT/SR[{tumor_id}:1]<5  & FORMAT/PR[{tumor_id}:1]<5" |
bcftools filter -e "SV_TOP_TIER > 2 & FORMAT/SR[{tumor_id}:1]<10 & FORMAT/PR[{tumor_id}:1]<10 & (BPI_AF[0] < 0.1 | BPI_AF[1] < 0.1)" |
bcftools view -s {params.sample} > {output.vcf}
''')

if not is_ffpe:
    rule reprioritize_rescued_svs:
        input:
            vcf = 'work/{batch}/purple/{batch}.purple.sv.vcf.gz',
        output:
            vcf = '{batch}/structural/{batch}-manta.vcf.gz',
            tbi = '{batch}/structural/{batch}-manta.vcf.gz.tbi',
        params:
            genome = run.genome_build
        group: "sv_after_purple"
        run:
            cmd = f'gunzip -c {input.vcf}' \
                f' | prioritize_sv -g {params.genome}' \
                f' | bcftools annotate -x "INFO/ANN"' \
                f' -Oz -o {output.vcf}' \
                f' && tabix -p vcf {output.vcf}'
            shell(cmd)
            before = count_vars(input.vcf)
            after = count_vars(output.vcf)
            assert before == after, (before, after)

else:
    rule copy_sv_vcf_ffpe_mode:
        input:
            vcf = rules.filter_sv_vcf.output.vcf
        output:
            vcf = '{batch}/structural/{batch}-manta.vcf.gz',
            tbi = '{batch}/structural/{batch}-manta.vcf.gz.tbi',
        group: "sv_after_purple"
        shell:
            'bgzip -c {input.vcf} > {output.vcf} && tabix -p vcf {output.vcf}'


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
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name
    group: "sv_after_purple"
    run:
        tumor_id = VCF(input.vcf).samples.index(params.sample)
        with open(output[0], 'w') as out:
            header = ["caller", "sample", "chrom", "start", "end", "svtype",
                      "split_read_support", "paired_support_PE", "paired_support_PR", "AF_BPI", "somaticscore",
                      "tier", "annotation",
                      'AF_PURPLE', 'CN_PURPLE', 'CN_change_PURPLE', 'Ploidy_PURPLE', 'PURPLE_status', 'START_BPI', 'END_BPI', 'ID', 'MATEID']
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

                data = ['manta', params.sample, rec.CHROM, rec.POS, rec.INFO.get('END', ''),
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
                        parse_info_field(rec, 'MATEID')
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
        fai = hpc.get_ref_file(run.genome_build, key='fa') + '.fai'
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
        fai = hpc.get_ref_file(run.genome_build, key='fa') + '.fai'
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

rule structural:
    input:
        expand(rules.bedpe.output, batch=batch_by_name.keys()),
        expand(rules.ribbon.output, batch=batch_by_name.keys()),
        expand(rules.prep_sv_tsv.output, batch=batch_by_name.keys()),
    output:
        temp(touch('log/structural.done'))
