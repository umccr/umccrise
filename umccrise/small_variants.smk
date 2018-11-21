#################
#### Somatic ####
from os.path import isfile, join
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.reference_data import get_key_genes_set
from umccrise import package_path
import cyvcf2
import toml
import csv
from ngs_utils.file_utils import which


localrules: small_variants


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
#   Remove DP<25 & AF<5% in tricky regions: gc15, gc70to75, gc75to80, gc80to85, gc85, low_complexity_51to200bp, low_complexity_gt200bp, non-GIAB confident, unless coding in cancer genes


# Preparations: annotate TUMOR_X and NORMAL_X fields for PCGR, remove non-standard chromosomes and mitochondria, remove non-PASSed calls
rule somatic_vcf_prep:
    input:
        vcf = lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].name}-ensemble-annotated.vcf.gz')
    output:
        vcf = 'work/{batch}/small_variants/somatic-ensemble-prep.vcf.gz',
        tbi = 'work/{batch}/small_variants/somatic-ensemble-prep.vcf.gz.tbi'
    group: "small_variants_1round"
    shell:
        'pcgr_prep {input.vcf} |'
        ' bcftools view -f.,PASS -Oz -o {output.vcf}'
        ' && tabix -p vcf {output.vcf}'

rule somatic_vcf_pon_annotate:
    input:
        vcf = rules.somatic_vcf_prep.output.vcf,
        tbi = rules.somatic_vcf_prep.output.tbi,
    params:
        genome_build = run.genome_build,
    output:
        vcf = 'work/{batch}/small_variants/somatic-ensemble-prep-pon.vcf.gz',
        tbi = 'work/{batch}/small_variants/somatic-ensemble-prep-pon.vcf.gz.tbi',
    group: "small_variants_1round"
    shell:
        'pon_anno {input.vcf} -o {output.vcf} --pon-dir ' + pon_dir + '&& tabix -p vcf {output.vcf}'

rule somatic_vcf_pon_1round:
    input:
        rules.somatic_vcf_pon_annotate.output.vcf
    output:
        vcf = 'work/{batch}/small_variants/somatic-ensemble-pon_round1.vcf.gz',
        tbi = 'work/{batch}/small_variants/somatic-ensemble-pon_round1.vcf.gz.tbi',
    params:
        pon_hits = 3
    group: "small_variants_1round"
    shell:
        'bcftools filter {input} -e "INFO/PoN_CNT>={params.pon_hits}" -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'

def _iter_vcf(inp_path, out_path, func, func_hdr=None):
    inp_vcf = cyvcf2.VCF(inp_path)
    ungz = get_ungz_gz(out_path)[0]
    if func_hdr is not None:
        func_hdr(inp_vcf)
    out_vcf = cyvcf2.Writer(ungz, inp_vcf)
    out_vcf.write_header()
    for rec in inp_vcf:
        res_rec = func(rec)
        if res_rec is not None:
            out_vcf.write_record(res_rec)
    out_vcf.close()
    shell('bgzip {ungz} && tabix -p vcf {out_path}')

rule somatic_vcf_pon_pass_keygenes:
    input:
        vcf = rules.somatic_vcf_pon_1round.output.vcf,
    output:
        vcf = 'work/{batch}/small_variants/somatic-ensemble-pon_round1-keygenes.vcf.gz',
    group: "small_variants_1round"
    run:
        genes = get_key_genes_set()
        def func(rec):
            if rec.INFO.get('ANN') is not None and rec.INFO['ANN'].split('|')[3] in genes:
                return rec
        _iter_vcf(input.vcf, output.vcf, func)

rule somatic_vcf_pcgr_ready:
    input:
        full_vcf = rules.somatic_vcf_pon_annotate.output.vcf,
        keygenes_vcf = rules.somatic_vcf_pon_1round.output.vcf,
    output:
        vcf = 'work/{batch}/small_variants/somatic-ensemble-pon_round1-ready.vcf.gz',
    group: "small_variants_1round"
    run:
        total_vars = int(subprocess.check_output(f'bcftools view -H {input.full_vcf} | wc -l', shell=True).strip())
        vcf = input.full_vcf if total_vars <= 500_000 else input.keygenes_vcf  # to avoid PCGR choking on too many variants
        shell(f'ln -s {basename(vcf)} {output.vcf}')
        shell(f'ln -s {basename(vcf)}.tbi {output.vcf}.tbi')

rule somatic_vcf_pcgr_1round:
    input:
        vcf = rules.somatic_vcf_pcgr_ready.output.vcf,
        pcgr_data = pcgr_data,
    output:
        tiers = 'work/{batch}/small_variants/pcgr_1round/{batch}-somatic.pcgr.snvs_indels.tiers.tsv',
    params:
        output_dir = 'work/{batch}/small_variants/pcgr_1round',
        genome = run.genome_build,
        sample_name = '{batch}-somatic',
        opt='--no-docker' if not which('docker') else '',
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
    shell:
        conda_cmd.format('pcgr') + which('pcgr') +
        ' {input.vcf} -g {params.genome} -o {params.output_dir} -s {params.sample_name} '
        '{params.opt} --pcgr-data {input.pcgr_data}'

rule somatic_vcf_pcgr_anno:
    input:
        vcf = rules.somatic_vcf_pcgr_ready.output.vcf,
        tiers = rules.somatic_vcf_pcgr_1round.output.tiers,
    output:
        vcf = 'work/{batch}/small_variants/somatic-ensemble-pon_round1-ready-pcgr_anno.vcf.gz',
    group: "small_variants_2round"
    run:
        gene_by_snp = dict()
        tier_by_snp = dict()
        with open(input.tiers) as f:
            reader = csv.DictReader(f, delimiter='\t', fieldnames=f.readline().strip().split('\t'))
            for row in reader:
                k = row['GENOMIC_CHANGE']
                gene_by_snp[k] = row['SYMBOL']
                tier_by_snp[k] = row['TIER']
        def func_hdr(vcf):
            vcf.add_info_to_header({'ID': 'PCGR_TIER', 'Description': 'TIER as reported by PCGR in .snvs_indels.tiers.tsv file', 'Type': 'String', 'Number': '1'})
            vcf.add_info_to_header({'ID': 'PCGR_GENE', 'Description': 'Gene symbol as reported by PCGR in .snvs_indels.tiers.tsv file', 'Type': 'String', 'Number': '1'})
        def func(rec):
            k = f'{rec.CHROM}:g.{rec.POS}{rec.REF}>{rec.ALT[0]}'
            rec.INFO['PCGR_TIER'] = tier_by_snp.get(k, '.')
            rec.INFO['PCGR_GENE'] = gene_by_snp.get(k, '.')
            return rec
        _iter_vcf(input.vcf, output.vcf, func, func_hdr)

rule prep_giab_bed:
    input:
        get_ref_file(run.genome_build, ['truth_sets', 'giab', 'bed'])
    output:
        'work/giab_conf.bed.gz'
    group: "small_variants_2round"
    shell:
        'cat {input} | bgzip -c > {output} && tabix -p bed {output}'

rule prep_anno_toml:
    input:
        tricky_bed = os.path.join(loc.extras, 'GRCh37_tricky.bed.gz'),
        giab_conf_bed = rules.prep_giab_bed.output[0],
    output:
        'work/tricky_vcfanno.toml'
    group: "small_variants_2round"
    params:
        toml_text  = lambda wc, input, output: f'''
[[annotation]]
file="{input.tricky_bed}"
names=["TRICKY"]
columns=[4]
ops=["self"]

[[annotation]]
file="{input.giab_conf_bed}"
names=["GIAB_CONF"]
columns=[3]
ops=["flag"]
'''.replace('\n', r'\\n').replace('"', r'\"'),
    shell:
        'printf "{params.toml_text}" > {output}'

rule somatic_vcf_regions_anno:
    input:
        vcf = rules.somatic_vcf_pcgr_anno.output.vcf,
        toml = rules.prep_anno_toml.output[0]
    output:
        vcf = 'work/{batch}/small_variants/somatic-ensemble-pon_round1-ready-pcgr_anno-regions_anno.vcf.gz',
    group: "small_variants_2round"
    shell:
        'vcfanno {input.toml} {input.vcf} | bgzip -c > {output} && tabix -p vcf -f {output}'

def _add_cyvcf2_filter(rec, filt):
    filters = rec.FILTER.split(';') if rec.FILTER else []
    filters.append(filt)
    rec.FILTER = ';'.join(filters)
    return rec

rule somatic_vcf_filter:
    input:
        vcf = rules.somatic_vcf_regions_anno.output.vcf,
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-ensemble-filt.vcf.gz',
    group: "small_variants_2round"
    run:
        def func_hdr(vcf):
            vcf.add_filter_to_header({'ID': 'gnomAD_common', 'Description': 'Occurs in gnomAD with frequency above 1%'})
            vcf.add_filter_to_header({'ID': 'PoN', 'Description': 'Panel of normals hits 1 or more'})
            vcf.add_filter_to_header({'ID': 'bad_promoter', 'Description': 'Indel overlapping bad promoter tricky region'})
            vcf.add_filter_to_header({'ID': 'LowQual_TRICKY', 'Description': 'DP<25 & AF<5%, and: GC<=15% or GC>=70 or low complexity region >51bp long'})
            vcf.add_filter_to_header({'ID': 'LowQual_GIAB_LCR', 'Description': 'DP<25 & AF<5%, and does not overlap GiaB high confidence regions'})
        def func(rec):
            t = rec.INFO['PCGR_TIER']
            int_tier = int(t.split()[1]) if 'TIER' in t else 5  # "TIER 2" -> 2
            # Keeping all variants with tier 1, 2, 3:
            # Tier 1 - variants of strong clinical significance
            # Tier 2 - variants of potential clinical significance
            # Tier 3 - variants of unknown clinical significance
            if 1 >= int_tier >= 3:
                return rec
            # Applying LC, PoN, depth and AF filters to tier 4 and non-coding:
            # Tier 4 - other coding variants
            # Noncoding variants
            # Remove gnomad_AF >0.02
            # Remove PoN_CNT>0         # {0 if issnp else 1}'
            # Remove indels in "bad_promoter" tricky regions
            # Remove DP<25 & AF<5% in tricky regions:
            #        gc15, gc70to75, gc75to80, gc80to85, gc85, low_complexity_51to200bp, low_complexity_gt200bp,
            #        non-GIAB confident,
            #    unless coding in cancer genes
            else:
                for rgn in ['GLOBAL', 'SAS', 'AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']:
                    if rec.INFO.get(f'{rgn}_AF_GNOMAD', 0) >= 0.01:
                        _add_cyvcf2_filter(rec, 'gnomAD_common')
                # second round of panel of normals
                pon = rec.INFO.get('PoN_CNT')
                if pon is not None and pon >= 1:
                    _add_cyvcf2_filter(rec, 'PoN')
                # removing indels in bad promoter regions
                tricky_set = set(rec.INFO.get('TRICKY', '').split(','))
                if not rec.is_snp and 'bad_promoter' in tricky_set:
                    _add_cyvcf2_filter(rec, 'bad_promoter')
                # removing low AF and low DP variants in low complexity regions
                if rec.INFO['TUMOR_DP'] < 25 and rec.INFO['TUMOR_AF'] < 0.05:
                    if tricky_set & {'gc15', 'gc70to75', 'gc75to80', 'gc80to85', 'gc85', 'low_complexity_51to200bp',
                                     'low_complexity_gt200bp'}:
                        _add_cyvcf2_filter(rec, 'LowQual_TRICKY')
                    if not rec.INFO.get('GIAB_CONF'):
                        _add_cyvcf2_filter(rec, 'LowQual_GIAB_LCR')
            return rec
        _iter_vcf(input.vcf, output.vcf, func, func_hdr)

rule somatic_vcf_filter_pass:
    input:
        vcf = rules.somatic_vcf_filter.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-somatic-ensemble-filt.pass.vcf.gz',
    group: "small_variants_2round"
    shell:
        'bcftools view -f.,PASS {input.vcf} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}'

##################
#### Germline ####
# Annotate any events found in ~200 cancer predisposition gene set.
rule germline_vcf_subset:  # {batch}
    input:
        vcf = lambda wc: join(run.date_dir, f'{batch_by_name[wc.batch].normal.name}{GERMLINE_SUFFIX}-ensemble-annotated.vcf.gz'),
    output:
        vcf = 'work/{batch}/small_variants/raw_normal-ensemble-predispose_genes.vcf.gz',
    params:
        ungz = lambda wc, output: get_ungz_gz(output[0])[0]
    group: "small_variants_germline"
    run:
        pcgr_toml_fpath = join(package_path(), 'pcgr', 'cpsr.toml')
        genes = [g for g in toml.load(pcgr_toml_fpath)['cancer_predisposition_genes']]
        def func(rec):
            if rec.INFO.get('ANN') is not None and rec.INFO['ANN'].split('|')[3] in genes:
                return rec
        _iter_vcf(input.vcf, output.vcf, func)

# Preparations: annotate TUMOR_X and NORMAL_X fields, remove non-standard chromosomes and mitochondria, remove non-PASSed calls.
# Suites for PCGR, but for all other processing steps too
rule germline_vcf_prep:
    input:
        vcf = rules.germline_vcf_subset.output.vcf
    output:
        vcf = '{batch}/small_variants/{batch}-normal-ensemble-predispose_genes.vcf.gz',
    group: "small_variants_germline"
    shell:
        'pcgr_prep {input.vcf} |'
        ' bcftools view -f.,PASS -Oz -o {output.vcf}'
        ' && tabix -p vcf {output.vcf}'


#############

rule small_variants:
    input:
        expand(rules.somatic_vcf_filter.output.vcf, batch=batch_by_name.keys()),
        expand(rules.somatic_vcf_filter_pass.output.vcf, batch=batch_by_name.keys()),
        expand(rules.germline_vcf_prep.output, batch=batch_by_name.keys())
    output:
        temp(touch('log/small_variants.done'))