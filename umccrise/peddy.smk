import csv
from os.path import isfile, join, dirname


PEDDY_OUT_EXTENSIONS = [
    ".background_pca.json", ".het_check.csv", ".pca_check.png",
    ".ped_check.png", ".ped_check.rel-difference.csv",
    ".ped_check.csv", ".peddy.ped", ".sex_check.csv", ".ped_check.png",
    ".html"]

rule create_ped_file:
    output:
        ped = 'work/{batch}/peddy.ped'
    group: 'peddy'
    params:
        sample_name = lambda wc: getattr(batch_by_name[wc.batch], wc.phenotype + 's')[0].name,
    run:
        header = ['#Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype']
                  # P003         CCR180088_NH18T002P003    -9             -9      0            2
        with open(output.ped, 'w') as f:
            writer = csv.writer(f, dialect='excel-tab')
            writer.writerow(header)
            family_id = wildcards.batch
            writer.writerow([family_id, params.sample_name, -9, -9,
                             {'male': 1, 'female': 2, 'unknown': 0}.get('unknown'),
                             {'tumor': 2, 'normal': 1}.get(wildcards.phenotype, 0)])

rule run_peddy:
    input:
        vcf = 'work/{batch}/small_variants/germline/{batch}-germline-PASS.vcf.gz',
        ped = 'work/{batch}/peddy.ped',
    output:
        dir = directory('work/{batch}/peddy')
    group: 'peddy'
    params:
        prefix = lambda wc, input, output: \
            join(output.dir, getattr(batch_by_name[wc.batch], wc.phenotype + 's')[0].name),
        genome = run.genome_build,
        phenotype = lambda wc: wc.phenotype,
    threads:
        threads_per_sample
    resources:
        mem_mb=8000
    shell:
        'mkdir -p {output.dir} ; peddy -p {threads} --sites {params.genome} --plot '
        '--prefix {params.prefix} {input.vcf} {input.ped}'

#
# rule preddy_multiqc_files:
#     input:
#         dirs = expand(rules.run_peddy.output.dir, phenotype=['tumor', 'normal'])
#     output:


rule peddy:
    input:
        expand(rules.run_peddy.output[0], batch=batch_by_name.keys(), phenotype=['normal'])
    output:
        temp(touch('log/peddy.done'))

