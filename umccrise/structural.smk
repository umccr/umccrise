"""
Structural variants
------------------
Re-do the CNV plots. This will need lots of love (drop gene names, make the scatterplot viable again, etc.).
"""
import itertools
import re

from cyvcf2 import VCF
from ngs_utils.utils import flatten
from ngs_utils.file_utils import safe_mkdir
from ngs_utils.reference_data import get_key_genes_txt, get_key_tsgenes_txt, \
    get_known_fusion_pairs, get_fusioncatcher_pairs, get_known_fusion_heads, get_known_fusion_tails
from vcf_stuff import count_vars, vcf_contains_field, iter_vcf
from bed_annotation import canon_transcript_per_gene

vcftobedpe = 'vcfToBedpe'


localrules: structural


def get_manta_path(b):
    return join(batch_by_name[b].tumor.dirpath, f'{batch_by_name[b].name}-sv-prioritize-manta.vcf.gz')
def get_sv_tsv_path(b):
    return join(batch_by_name[b].tumor.dirpath, f'{batch_by_name[b].name}-sv-prioritize.tsv')
if not all(isfile(get_manta_path(b)) for b in batch_by_name.keys()):
    # CWL?
    def get_manta_path(b):
        return join(run.date_dir, batch_by_name[b].tumor.name + '-manta-prioritized.vcf.gz')
    def get_sv_tsv_path(b):
        return join(run.date_dir, batch_by_name[b].tumor.name + '-prioritize.tsv')

    if not all(isfile(get_manta_path(b)) for b in batch_by_name.keys()):
        critical('Could not find manta files for all batches neither under sample folders as '
                 '<tumor>/<batch>-sv-prioritize-manta.vcf.gz (conventional bcbio), nor in the project folder as'
                 'project/<tumor>-manta-prioritized.vcf.gz (CWL bcbio).')

# TODO: ? Rerun SnpEFF as well to target canonical transcripts, so we don't miss intergenic variants touching non-canonical transripts?

rule sv_subset_to_canonical:
    input:
        vcf = lambda wc: get_manta_path(wc.batch),
    output:
        vcf = 'work/{batch}/structural/canonical_ann/{batch}-manta.vcf.gz',
        tbi = 'work/{batch}/structural/canonical_ann/{batch}-manta.vcf.gz.tbi',
    group: "sv_vcf"
    run:
        assert vcf_contains_field(input.vcf, 'INFO/ANN'), f'Manta {input.vcf} must be annotated with SnpEff'

        principal_by_gids = canon_transcript_per_gene(run.genome_build, use_gene_id=True, only_principal=True)
        all_princ_transcripts = set(principal_by_gids.values())

        def proc_line(rec, **kwargs):
            ann = rec.INFO.get('ANN')
            if ann is None:
                return rec

            if isinstance(ann, str):
                anns = ann.split(',')
            else:
                anns = ann

            new_anns = []
            for ann_line in anns:
                ann_fields = ann_line.split('|')
                assert len(ann_fields) >= 11, f'rec: {rec}, ann_line: {ann_line}'
                allele, effect, impact, genename, geneid, feature, featureid, biotype, rank, c_change, p_change = ann_fields[:11]

                if feature == 'chromosome':
                    new_anns.append(ann_line)

                elif feature == 'transcript':
                    feature_ids = set(featureid.split('&'))
                    assert all(re.fullmatch(r'ENST\d+', fid) for fid in feature_ids), ann
                    if not feature_ids - all_princ_transcripts:
                        new_anns.append(ann_line)

                else:  # 'interaction'?
                    # Splitting 'interaction' feature ids like the following:
                    # "3UVN:C_107-D_1494:ENST00000358625-ENST00000262519"
                    # -> {'3UVN', 'C_107', 'D_1494', 'ENST00000262519', 'ENST00000358625'}
                    feature_ids = set(fid for fid in flatten(fids.split('-') for fids in featureid.split(':')) if fid.startswith('ENST'))
                    assert all(re.fullmatch(r'ENST\d+', fid) for fid in feature_ids), ann
                    if not feature_ids - all_princ_transcripts:
                        new_anns.append(ann_line)

            if new_anns:
                if isinstance(ann, str):
                    new_ann = ','.join(new_anns)
                else:
                    new_ann = new_anns
                rec.INFO['ANN'] = new_ann
            else:
                del rec.INFO['ANN']

            # Also remove older SIMPLE_ANN as we are rerunning for canonical transcripts
            if rec.INFO.get('SIMPLE_ANN') is not None:
                del rec.INFO['SIMPLE_ANN']
                del rec.INFO['SV_HIGHEST_TIER']

            return rec

        iter_vcf(input.vcf, output.vcf, proc_line)
        before = count_vars(input.vcf)
        after = count_vars(output.vcf)
        assert before == after, (before, after)

rule prep_sv_prio_lists:
    input:
        fusion_pairs = get_known_fusion_pairs(),
        fusion_heads = get_known_fusion_heads(),
        fusion_tails = get_known_fusion_tails(),
    output:
        pairs = 'work/fusion_pairs.txt',
        promisc = 'work/fusion_promisc.txt',
    run:
        with open(input.fusion_pairs) as f, open(output.pairs, 'w') as out:
            for l in f:
                l = l.strip().replace('"', '')
                if l:
                    g1, g2 = l.split(',')[0:2]
                    if g1 and g2 and g1 != 'H_gene':
                        out.write(f'{g1},{g2}\n')
        with open(input.fusion_heads) as f1, open(input.fusion_tails) as f2, open(output.promisc, 'w') as out:
            for l in itertools.chain(f1, f2):
                l = l.strip().replace('"', '')
                if l:
                    gene = l.split(',')[0]
                    if gene and gene != 'gene':
                        out.write(f'{gene}\n')

rule sv_prioritize:
    input:
        vcf = rules.sv_subset_to_canonical.output.vcf,
        known_pairs = rules.prep_sv_prio_lists.output.pairs,
        known_pairs_fusioncatcher = get_fusioncatcher_pairs(),
        known_promisc = rules.prep_sv_prio_lists.output.promisc,
        cancer_genes = get_key_genes_txt(),
        tsgenes = get_key_tsgenes_txt(),
    output:
        vcf = 'work/{batch}/structural/prioritize/{batch}-manta.vcf.gz',
        tbi = 'work/{batch}/structural/prioritize/{batch}-manta.vcf.gz.tbi',
    group: "sv_vcf"
    run:
        cmd = f'cat {input.vcf}'
        # remove previous annotation
        filts_to_remove = [f'{f}' for f in ['INFO/SIMPLE_ANN', 'INFO/SV_HIGHEST_TIER',
                                            'FILTER/Intergenic', 'FILTER/MissingAnn', 'FILTER/REJECT']
                           if vcf_contains_field(input.vcf, f)]
        if filts_to_remove:
            cmd += f' | bcftools annotate -x "' + ','.join(f'{f}' for f in filts_to_remove) + '"'

        cmd += (f' | simple_sv_annotation.py - '
                f'--tier2_fusion_pairs {input.known_pairs_fusioncatcher} '
                f'--known_fusion_pairs {input.known_pairs} '
                f'--known_fusion_promiscuous {input.known_promisc} '
                f'--gene_list {input.cancer_genes} '
                f'--tumor_suppressors {input.tsgenes} '
                f'-o - | bgzip -c > {output.vcf} && tabix -p vcf {output.vcf}')
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
        cmd = f'bcftools view -f.,PASS {input.vcf} | bcftools annotate -x INFO/ANN'
        if count_vars(input.vcf, filter='.,PASS') > 1000:
            cmd += ' | bcftools filter -i "SV_HIGHEST_TIER <= 3"'
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
        xmx = 2800,
        tmp_dir = '{batch}/structural/maybe_bpi/tmp_dir'
    resources:
        mem_mb = 3000
    run:
        if not vcf_contains_field(input.vcf, 'BPI_AF', 'INFO'):
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
        else:
            shell('cp {input.vcf} {output.vcf}')

# Keep all with read support above 10x; or allele frequency above 10%, but only if read support is above 5x
rule filter_sv_vcf:
    input:
        vcf = rules.sv_maybe_bpi.output.vcf
    output:
        vcf = '{batch}/structural/{batch}-manta.vcf'
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name
    group: "sv_vcf"
    run:
        print(f'VCF samples: {VCF(input.vcf).samples}')
        print(f'Bcbio batch tumor name: {batch_by_name[wildcards.batch].tumor.name}')
        tumor_id = VCF(input.vcf).samples.index(params.sample)
        tumor_id = VCF(input.vcf).samples.index(batch_by_name[wildcards.batch].tumor.name)
        print(f'Derived tumor VCF index: {tumor_id}')
        shell('''
bcftools view -f.,PASS {input.vcf} | \
bcftools filter -e "FORMAT/SR[{tumor_id}:1]<5  & FORMAT/PR[{tumor_id}:1]<5" | \
bcftools filter -e "FORMAT/SR[{tumor_id}:1]<10 & FORMAT/PR[{tumor_id}:1]<10 & (BPI_AF[0] < 0.1 | BPI_AF[1] < 0.1)" \
> {output.vcf}
''')

# Produce a TSV file for further analysis in Rmd
# caller  sample                chrom   start       end         svtype  lof  annotation                                                                 split_read_support  paired_support_PE  paired_support_PR  somaticscore  tier
# manta   PRJ180253_E190-T01-D  1       161513440   161595209   DUP          DUP|GENE_FUSION|FCGR2B&RP11-25K21.6|ENSG00000273112|NOT_PRIORITISED|3,...                                         67,8
# manta   PRJ180253_E190-T01-D  8       33320739    33321344    DEL          DEL|UPSTREAM_GENE_VARIANT|FUT10|ENST00000518672|NOT_PRIORITISED|3          76,7                                   33,1
# manta   PRJ180253_E190-T01-D  11      118802640   118803304   DEL          DEL|DOWNSTREAM_GENE_VARIANT|RN7SL688P|ENST00000471754|NOT_PRIORITISED|3    61,8                                   07,2
rule prep_sv_tsv:
    input:
        vcf = rules.filter_sv_vcf.output.vcf
    output:
        '{batch}/structural/{batch}-manta.tsv'
    params:
        sample = lambda wc: batch_by_name[wc.batch].tumor.name
    group: "sv_vcf"
    run:
        tumor_id = VCF(input.vcf).samples.index(params.sample)
        with open(output[0], 'w') as out:
            header = ["caller", "sample", "chrom", "start", "end", "svtype",
                      "split_read_support", "paired_support_PE", "paired_support_PR", "BPI_AF", "somaticscore", "tier",
                      "annotation", "lof"]
            out.write('\t'.join(header) + '\n')
            for rec in VCF(input.vcf):
                data = ['manta', params.sample, rec.CHROM, rec.POS, rec.INFO.get('END', ''),
                        rec.INFO['SVTYPE'],
                        ','.join(map(str, rec.format('SR')[tumor_id])) if 'SR' in rec.FORMAT else '',
                        ','.join(map(str, rec.format('PE')[tumor_id])) if 'PE' in rec.FORMAT else '',
                        ','.join(map(str, rec.format('PR')[tumor_id])) if 'PR' in rec.FORMAT else '',
                        ','.join(map(str, rec.INFO['BPI_AF'])) if rec.INFO.get('BPI_AF') else '',
                        rec.INFO.get('SOMATICSCORE', ''),
                        rec.INFO.get('SV_HIGHEST_TIER', ''),
                        rec.INFO.get('SIMPLE_ANN', ''),
                        rec.INFO.get('LOF_GENES', ''),
                        ]
                out.write('\t'.join(map(str, data)) + '\n')

# At least for the most conservative manta calls, generate a file for viewing in Ribbon
rule ribbon_filter_manta:
    input:
        manta_vcf = rules.filter_sv_vcf.output.vcf
    output:
        'work/{batch}/structural/ribbon/manta.vcf'
    group: "sv_vcf"
    shell:
        'bcftools view {input.manta_vcf} > {output}'

rule ribbon_filter_vcfbedtope_starts:
    input:
        bed = rules.ribbon_filter_manta.output[0],
        fai = ref_fa + '.fai'
    output:
        'work/{batch}/structural/ribbon/manta-starts.bed'
    params:
        vcftobedpe = vcftobedpe
    group: "sv_vcf"
    shell:
        'cat {input.bed} | {params.vcftobedpe}'
        ' | cut -f 1-3'
        ' | bedtools slop -b 5000 -i stdin -g {input.fai}'
        ' > {output}'

rule ribbon_filter_vcfbedtope_ends:
    input:
        bed = rules.ribbon_filter_manta.output[0],
        fai = ref_fa + '.fai'
    output:
        'work/{batch}/structural/ribbon/manta-ends.bed'
    params:
        vcftobedpe = vcftobedpe
    group: "sv_vcf"
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
    group: "sv_vcf"
    shell:
        'cat {input.starts} {input.ends} | bedtools sort -i stdin | bedtools merge -i stdin > {output}'


#### Convert matna VCF to bedpe ####
rule bedpe:
    input:
        manta_vcf = rules.filter_sv_vcf.output.vcf
    output:
        '{batch}/structural/{batch}-manta.bedpe'
    params:
        vcftobedpe = vcftobedpe
    group: "sv_vcf"
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
