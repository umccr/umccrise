from os.path import isfile, join
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.vcf_utils import get_sample_ids, get_sample_names
from umccrise import package_path
import toml
from vcf_stuff import iter_vcf
import cyvcf2
import yaml
from vcf_stuff.filtering import add_cyvcf2_filter


localrules: sage


rule run_sage:
    input:
        tumor_bam  = lambda wc: batch_by_name[wc.batch].tumor.bam,
        normal_bam = lambda wc: batch_by_name[wc.batch].normal.bam,
        coding_bed = get_ref_file(run.genome_build, 'coding_regions'),
        ref_fa = ref_fa,
        hotspots_vcf = get_ref_file(run.genome_build, key='hotspots'),
    output:
        sage_vcf = 'work/{batch}/sage/call/{batch}.vcf.gz',
        sage_tbi = 'work/{batch}/sage/call/{batch}.vcf.gz.tbi',
    params:
        jar = join(package_path(), 'jars', 'sage-1.0-jar-with-dependencies.jar'),
        rundir = 'work/{batch}/purple',
        outdir = 'work/{batch}/purple',
        normal_sname = lambda wc: batch_by_name[wc.batch].normal.name,
        tumor_sname  = lambda wc: batch_by_name[wc.batch].tumor.name,
        xms = 2000,
        xmx = 19000,
    resources:
        mem_mb = 20000
    group: "sage"
    shell:
        'java -Xms{params.xms}m -Xmx{params.xmx}m -cp {params.jar} com.hartwig.hmftools.sage.SageHotspotApplication '
        '-tumor {params.tumor_sname} -tumor_bam {input.tumor_bam} '
        '-reference {params.normal_sname} -reference_bam {input.normal_bam} '
        '-known_hotspots <(bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\n" {input.hotspots_vcf}) '
        '-coding_regions {input.coding_bed} '
        '-ref_genome {input.ref_fa} '
        '-out {output.sage_vcf} '
        '&& tabix -f -p vcf {output.sage_vcf}'


rule sage_rename_anno:
    input:
        sage_vcf = rules.run_sage.output.sage_vcf,
    output:
        sage_vcf = 'work/{batch}/sage/rename_anno/{batch}.vcf.gz',
        sage_tbi = 'work/{batch}/sage/rename_anno/{batch}.vcf.gz.tbi',
    group: "sage"
    shell: """
bcftools view {input.sage_vcf} | \
sed 's/HOTSPOT/SAGE_HOTSPOT/g' | \
sed 's/Hotspot Type: known, inframe/SAGE Hotspot Type: known, inframe'/g | \
bgzip -c > {output.sage_vcf} && tabix -p vcf {output.sage_vcf}
"""


rule sage_reorder_samples:
    input:
        sage_vcf = rules.sage_rename_anno.output.sage_vcf,
        sage_tbi = rules.sage_rename_anno.output.sage_tbi,
        vcf = rules.somatic_vcf_pass_sort.output.vcf,
    output:
        sage_vcf = 'work/{batch}/sage/reoder_samples/{batch}.vcf.gz',
        sage_tbi = 'work/{batch}/sage/reoder_samples/{batch}.vcf.gz.tbi',
    group: "sage"
    run:
        tumor_index, normal_index = get_sample_ids(input.vcf)
        assert sorted([tumor_index, normal_index]) == [0, 1]
        sample_in_order = [None, None]
        tumor_name, normal_name = get_sample_names(input.vcf)
        sample_in_order[tumor_index] = tumor_name
        sample_in_order[normal_index] = normal_name
        shell(f'bcftools view -s {",".join(sample_in_order)} {input.sage_vcf} -Oz -o {output.sage_vcf} '
              f'&& tabix -p vcf {output.sage_vcf}')


rule sage_pass:
    input:
        sage_vcf = rules.sage_reorder_samples.output.sage_vcf,
        sage_tbi = rules.sage_reorder_samples.output.sage_tbi,
    output:
        sage_vcf = 'work/{batch}/sage/pass/{batch}.vcf.gz',
        sage_tbi = 'work/{batch}/sage/pass/{batch}.vcf.gz.tbi',
    group: "sage"
    shell:
        'bcftools view -f.,PASS {input.sage_vcf} -Oz -o {output.sage_vcf} '
        '&& tabix -p vcf {output.sage_vcf}'


rule sage_pass_novel:
    input:
        sage_vcf = rules.sage_pass.output.sage_vcf,
        sage_tbi = rules.sage_pass.output.sage_tbi,
        vcf = rules.somatic_vcf_pass_sort.output.vcf,
    output:
        sage_vcf = 'work/{batch}/sage/pass_novel/{batch}.vcf.gz',
        sage_tbi = 'work/{batch}/sage/pass_novel/{batch}.vcf.gz.tbi',
    group: "sage"
    shell:
        'bcftools isec {input.sage_vcf} {input.vcf} -C -w1 -Oz -o {output.sage_vcf} '
        '&& tabix -p vcf {output.sage_vcf}'


rule add_novel_sage_calls:
    input:
        vcf = rules.somatic_vcf_pass_sort.output.vcf,
        tbi = rules.somatic_vcf_pass_sort.output.vcf,
        sage_vcf = rules.sage_pass_novel.output.sage_vcf,
        sage_tbi = rules.sage_pass_novel.output.sage_tbi,
    output:
        vcf = 'work/{batch}/small_variants/sage_add/{batch}-somatic-' + run.somatic_caller + '-sage.vcf.gz',
        tbi = 'work/{batch}/small_variants/sage_add/{batch}-somatic-' + run.somatic_caller + '-sage.vcf.gz.tbi',
    group: "sage"
    run:
         shell('bcftools concat -a {input.sage_vcf} {input.vcf} -Oz -o {output.vcf} '
               '&& tabix -p vcf {output.vcf}')
         assert len(cyvcf2.VCF(output.vcf).samples) == 2


rule sort_saged:
    input:
        vcf = rules.add_novel_sage_calls.output.vcf,
    output:
        vcf = 'work/{batch}/small_variants/sorted/{batch}-somatic-' + run.somatic_caller + '.vcf.gz',
        tbi = 'work/{batch}/small_variants/sorted/{batch}-somatic-' + run.somatic_caller + '.vcf.gz.tbi',
    group: "sage"
    shell:
        '(bcftools view -h {input.vcf} ; bcftools view -H -f.,PASS {input.vcf} | sort -k1,1V -k2,2n) | '
        'bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'


rule annotate_from_sage:
    input:
        vcf = rules.sort_saged.output.vcf,
        tbi = rules.sort_saged.output.tbi,
        sage_vcf = rules.sage_reorder_samples.output.sage_vcf,
        sage_tbi = rules.sage_reorder_samples.output.sage_tbi,
    output:
        vcf = 'work/{batch}/small_variants/sage_anno/{batch}-somatic-' + run.somatic_caller + '-sage.vcf.gz',
        tbi = 'work/{batch}/small_variants/sage_anno/{batch}-somatic-' + run.somatic_caller + '-sage.vcf.gz.tbi',
    group: "sage"
    run:
        sage_calls = dict()
        for rec in cyvcf2.VCF(input.sage_vcf):
            key = (rec.CHROM, rec.POS, rec.REF, rec.ALT[0])
            sage_calls[key] = rec

        def proc_hdr(vcf):
            vcf.add_filter_to_header({'ID': 'SAGE_lowconf', 'Description': 'SAGE assigned low confidence to this call'})

        def proc_rec(rec, tumor_index, normal_index):
            key = (rec.CHROM, rec.POS, rec.REF, rec.ALT[0])
            sage_call = sage_calls.get(key)
            if sage_call is not None:
                # pcgr_prep should handle the existing SAGE fields, so we just need to
                # figure out how to populate FORMAT fields from cyvcf2
                if not sage_call.FILTER or sage_call.FILTER == 'PASS':
                    assert sage_call.INFO.get('SAGE_HOTSPOT') is not None, sage_call
                    rec.INFO['SAGE_HOTSPOT'] = sage_call.INFO['SAGE_HOTSPOT']
                    rec.FILTER = 'PASS'
                else:
                    add_cyvcf2_filter(rec, 'SAGE_lowconf')
                rec.set_format('DP', sage_call.format('DP'))
                rec.set_format('AD', sage_call.format('AD'))
            return rec

        tumor_index, normal_index = get_sample_ids(input.vcf)
        iter_vcf(input.vcf, output.vcf, proc_rec, proc_hdr=proc_hdr, tumor_index=tumor_index, normal_index=normal_index)

#         with open(input.sage_vcf + '_toml', 'w') as f:
#             f.write(f"""[[annotation]]
# file="{input.sage_vcf}"
# fields = ["FILTER", "AF", "SAGE_HOTSPOT"]
# names = ["SAGE_FILTER", "AF", "SAGE_HOTSPOT"]
# ops = ["self", "self", "self"]
# """)
#         shell(f'vcfanno {input.sage_vcf}_toml {input.vcf} | bgzip -c > {output.vcf} && tabix -p vcf {output.vcf}')


#############

rule sage:
    input:
        expand(rules.annotate_from_sage.output.vcf, batch=batch_by_name.keys()),
    output:
        temp(touch('log/sage.done'))

