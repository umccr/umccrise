"""
PCGR
-------------
Prepare somatic, germline variant files, and configuration TOMLs for PCGR; tarball and upload to the AWS instance
"""
localrules: pcgr_somatic_vcf, pcgr_germline_vcf, pcgr_cns, pcgr_symlink_somatic, pcgr_symlink_germline, pcgr_prep, pcgr


import subprocess
from ngs_utils.file_utils import which


######################
### Preparing inputs

# PCGR struggles with anything but the basic chromosome setup.
# It also ignores any variant not marked `PASS` so might as well remove others to save on transfer times.
rule pcgr_somatic_vcf:
    input:
        vcf = rules.somatic_vcf_pon_pass.output.vcf,
        keygenes_vcf = rules.somatic_vcf_pon_pass_keygenes.output.vcf,
    output:
        vcf = '{batch}/pcgr/input/{batch}-somatic.vcf.gz',
        tbi = '{batch}/pcgr/input/{batch}-somatic.vcf.gz.tbi'
    run:
        total_vars = int(subprocess.check_output(f'bcftools view -H {input.vcf} | wc -l', shell=True).strip())
        vcf = input.vcf if total_vars <= 500_000 else input.keygenes_vcf  # to avoid PCGR choking on too many variants
        shell(f'cp {vcf} {output.vcf}')
        shell(f'cp {vcf}.tbi {output.vcf}.tbi')

rule pcgr_germline_vcf:
    input:
        vcf = rules.germline_vcf_prep.output.vcf,
        tbi = rules.germline_vcf_prep.output.tbi
    output:
        vcf = '{batch}/pcgr/input/{batch}-normal.vcf.gz',
        tbi = '{batch}/pcgr/input/{batch}-normal.vcf.gz.tbi'
    shell:
        'cp {input.vcf} {output.vcf} && cp {input.tbi} {output.tbi}'

# PCGR also wants a slightly different format for the CNS data:
rule pcgr_cns:
    input:
        '{batch}/purple/{batch}.purple.cnv',
        # lambda wc: join(batch_by_name[wc.batch].tumor.dirpath, f'{batch_by_name[wc.batch].name}-cnvkit-call.cns')
    output:
        '{batch}/pcgr/input/{batch}-somatic-cna.tsv'
    shell:
        'echo -e "Chromosome\\tStart\\tEnd\\tSegment_Mean" > {output} && ' \
        'cat {input} | grep -v ^# | grep -v ^GL | cut -f1,2,3,4 >> {output}'

if pcgr_data:
    ######################
    ### Running PCGR
    rule run_pcgr_local_somatic:
        input:
            vcf = rules.pcgr_somatic_vcf.output.vcf,
            tbi = rules.pcgr_somatic_vcf.output.tbi,
            cns = rules.pcgr_cns.output[0],
            pcgr_data = pcgr_data
        output:
            '{batch}/pcgr/{batch}-somatic.pcgr.html'
        params:
            output_dir = '{batch}/pcgr',
            genome_build = run.genome_build,
            sample_name = '{batch}-somatic',
            opt='--no-docker' if not which('docker') else ''
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 20000
            # TODO: memory based on the mutation number. E.g. over 455k tumor mutations need over 10G
        shell:
            conda_cmd.format('pcgr') +
            'pcgr {input.vcf} {input.cns} -g {params.genome_build} -o {params.output_dir} -s {params.sample_name} '
            '{params.opt} --pcgr-data {input.pcgr_data}'

    rule run_pcgr_local_germline:
        input:
            vcf = rules.pcgr_germline_vcf.output.vcf,
            tbi = rules.pcgr_germline_vcf.output.tbi,
            pcgr_data = pcgr_data
        output:
            '{batch}/pcgr/{batch}-normal.cpsr.html'
        params:
            output_dir = '{batch}/pcgr',
            genome_build = run.genome_build,
            sample_name = '{batch}-normal',
            opt='--no-docker' if not which('docker') else ''
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 20000
        shell:
            conda_cmd.format('pcgr') +
            'pcgr {input.vcf} -g {params.genome_build} -o {params.output_dir} -s {params.sample_name} --germline '
            '{params.opt} --pcgr-data {input.pcgr_data} ' \
            # '|| echo "Failed to run CPSR" >> {output}'

    rule pcgr_symlink_somatic:
        input:
            rules.run_pcgr_local_somatic.output[0]
        output:
            '{batch}/{batch}-somatic.pcgr.html'
        params:
            source = 'pcgr/{batch}-somatic.pcgr.html'
        shell:
            'ln -s {params.source} {output}'

    rule pcgr_symlink_germline:
        input:
            rules.run_pcgr_local_germline.output[0]
        output:
            '{batch}/{batch}-normal.cpsr.html'
        params:
            source = 'pcgr/{batch}-normal.cpsr.html'
        shell:
            'ln -s {params.source} {output}'

    ######################
    ###  Target rule

    rule pcgr:
        input:
            expand(rules.pcgr_symlink_somatic.output, batch=batch_by_name.keys()),
            expand(rules.pcgr_symlink_germline.output, batch=batch_by_name.keys())
        output:
            temp(touch('log/pcgr.done'))

else:

    ######################
    ###  Target rule

    rule pcgr_prep:
        input:
            expand(rules.pcgr_somatic_vcf.output, batch=batch_by_name.keys()),
            expand(rules.pcgr_cns.output, batch=batch_by_name.keys()),
            expand(rules.pcgr_germline_vcf.output, batch=batch_by_name.keys())
        output:
            temp(touch('log/pcgr_prep.done'))





#################################################
############# Rules for AWS PCGR ################

# rule prep_tomls:
#     input:
#         somatic = join(package_path(), 'pcgr', 'pcgr_configuration_somatic.toml'),
#         germline = join(package_path(), 'pcgr', 'pcgr_configuration_normal.toml')
#     output:
#         somatic = '{batch}/pcgr/input/{batch}' + uuid_suffix + '-somatic.toml',
#         germline = '{batch}/pcgr/input/{batch}' + uuid_suffix + '-normal.toml'
#     shell:
#         'cp {input.somatic} {output.somatic} && cp {input.germline} {output.germline}'
#
# ######################
# ###  Making tarballs
#
# rule somatic_tar_gz:
#     input:
#         vcf = rules.pcgr_somatic_vcf.output.vcf,
#         tbi = rules.pcgr_somatic_vcf.output.tbi,
#         cns = rules.pcgr_cns.output[0],
#         toml = rules.prep_tomls.output.somatic
#     output:
#         '{batch}/pcgr/input/{batch}' + uuid_suffix + '-somatic.tar.gz'
#     params:
#         basedir = '{batch}/pcgr/input/',
#         vcf = lambda wc, input: basename(input.vcf),
#         tbi = lambda wc, input: basename(input.tbi),
#         cns = lambda wc, input: basename(input.cns),
#         toml = lambda wc, input: basename(input.toml)
#     shell:
#         'tar -czf {output} -C {params.basedir} {params.vcf} {params.tbi} {params.cns} {params.toml}'
#
# rule germline_tar_gz:
#     input:
#         vcf = rules.pcgr_germline_vcf.output.vcf,
#         tbi = rules.pcgr_germline_vcf.output.tbi,
#         toml = rules.prep_tomls.output.germline
#     output:
#         '{batch}/pcgr/input/{batch}' + uuid_suffix + '-normal.tar.gz'
#     params:
#         basedir = '{batch}/pcgr/input/',
#         vcf = lambda wc, input: basename(input.vcf),
#         tbi = lambda wc, input: basename(input.tbi),
#         toml = lambda wc, input: basename(input.toml)
#     shell:
#         'tar -czf {output} -C {params.basedir} {params.vcf} {params.tbi} {params.toml}'
#
# ######################
# ###  Uploading tarballs
#
# upload_cmd = upload_proxy + 'aws s3 ls s3://pcgr/{params.fname}' \
#     ' && touch {output} && echo "Tarball already uploaded, or cannot access the bucket"' \
#     ' || ' + upload_proxy + 'aws s3 cp {input} s3://pcgr && touch {output}'
#
# rule upload_somatic_to_pcgr:
#     priority: 50
#     input:
#         rules.somatic_tar_gz.output[0]
#     output:
#         '{batch}/pcgr/input/upload-somatic.done'
#     params:
#         fname = lambda w, output, input: basename(input[0])
#     shell:
#         upload_cmd
#
# rule upload_germline_to_pcgr:
#     priority: 50
#     input:
#         rules.germline_tar_gz.output[0]
#     output:
#         '{batch}/pcgr/input/upload-normal.done'
#     params:
#         fname = lambda w, output, input: basename(input[0])
#     shell:
#         upload_cmd

######################
###  Target rules

# rule pcgr_prep:
#     input:
#         expand(rules.somatic_tar_gz.output, batch=batch_by_name.keys())
#     output:
#         temp(touch('pcgr_prep.done'))

# ###################
# ### Downloading results
#
# rule download_pcgr:
#     input:
#         expand(rules.igv_upload.output, phenotype=['tumor', 'normal'], batch=batch_by_name.keys())\
#             if not config.get('pcgr_download')\
#             else []  # unless download_pcgr is set explicitly, run only after done with IGV - to make sure instance is finished
#     output:
#          # Even if failed, create some output anyway - we cannot be sure that PCGR has
#          # finished, so shouldn't just fail the run because of that.
#          temp('{batch}/pcgr/{batch}' + uuid_suffix + '-{phenotype}.pcgr.download_status')
#     # output:
#     #     (('{batch}/pcgr/{batch}' + uuid_suffix + '-{phenotype}.pcgr.html')
#     #      if config.get('download_pcgr')
#     #      else ('{batch}/pcgr/{batch}' + uid_suffix + '-{phenotype}.pcgr.download_status'))  # if failed, create some output anyway - we cannot be sure that PCGR has finished, so shouldn't just crash
#     params:
#         targz_folder = 'work/{batch}/pcgr',
#         targz_fname = '{batch}' + uuid_suffix + '-{phenotype}-output.tar.gz',
#         untar_output_dirname = '{batch}' + uuid_suffix + '-{phenotype}-output',
#         final_output_folder = '{batch}/pcgr'
#     shell:
#         upload_proxy + 'aws s3 ls s3://pcgr/{params.targz_fname}'
#         ' && mkdir -p {params.targz_folder} '
#         ' && ' + upload_proxy + 'aws s3 cp s3://pcgr/{params.targz_fname} {params.targz_folder}'
#         ' && tar -xzf {params.targz_folder}/{params.targz_fname} -C {params.targz_folder}'
#         ' && mv {params.targz_folder}/{params.untar_output_dirname}/* {params.final_output_folder}'
#         ' && echo "ok" > {output}'
#         ' || echo "Cannot find the result tarball at s3://pcgr/{params.targz_fname}. PCGR is not finished?" >&2'
#         ' && echo "failed" > {output}'
#
# ### Target rule for downloading
# rule pcgr_download:
#     input:
#         expand(rules.download_pcgr.output, phenotype=['somatic', 'normal'], batch=batch_by_name.keys())
#     output:
#         # 'if [ -f {input[0]} -a -f {input[1]} ] ; then rm {input} ; fi'
#         temp(touch('pcgr_download.done'))
