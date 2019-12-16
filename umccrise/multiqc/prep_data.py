import os
import shutil
import subprocess
import tarfile
from inspect import getsourcefile
from os.path import join, dirname, abspath, pardir, isfile, isdir, islink, basename, relpath, getmtime, exists
from time import localtime, strftime
import yaml
import re

from ngs_utils import logger
from ngs_utils.bcbio import BcbioProject
from ngs_utils.call_process import run, run_simple
from ngs_utils.file_utils import verify_file, safe_mkdir, can_reuse, file_transaction
from ngs_utils.logger import info, warn, err, critical, timestamp, debug


def make_report_metadata(bcbio_proj, tumor_sample, normal_sample, base_dirpath, analysis_dir=None,
                         prog_versions_fpath=None, data_versions_fpath=None, new_dir_for_versions=None):
    conf = dict()
    conf['umccr'] = dict()
    additional_files = []

    conf['umccr']['tumor_name'] = tumor_sample
    conf['umccr']['normal_name'] = normal_sample

    conf['umccr']['is_rnaseq'] = bcbio_proj.is_rnaseq

    # General links
    conf['title'] = bcbio_proj.project_name
    conf['umccr']['run_section'] = get_run_info(bcbio_proj, base_dirpath, analysis_dir=analysis_dir,
                                                prog_versions_fpath=prog_versions_fpath,
                                                data_versions_fpath=data_versions_fpath,
                                                new_dir_for_versions=new_dir_for_versions)

    # if bcbio_proj.is_rnaseq:
    #     conf['umccr']['expression_links'] = _rna_general_links(bcbio_proj, base_dirpath)
    # else:
    #     mutations_links = _dna_general_links(bcbio_proj, base_dirpath)
    #     conf['umccr']['mutations_links'] = mutations_links
    #     if call_vis_html_fpath:
    #         call_vis_link = _make_link(call_vis_html_fpath, base_dirpath, 'call visualisation', blank=True)
    #         mutations_links.append(call_vis_link)
    #     if combined_ngs_rep_html_fpath:
    #         link = _make_link(combined_ngs_rep_html_fpath, base_dirpath, 'known mutations')
    #         mutations_links.append(link)

    return conf, additional_files


def _extract_qc_file(fp, new_mq_data_dir, final_dir, f_by_fp=None):
    """ Extracts QC file `fp` either by copying from `final_dir` (native bcbio),
        or from tar.gz file `tar_path` (CWL bcbio). Writes into a new file at new_mq_data_dir
    """
    if fp.startswith('report/metrics/'):
        fp = fp.replace('report/metrics/', 'project/multiqc/')  # for CWL _bcbio.txt files

    dst_fp = join(new_mq_data_dir, fp)

    fp_in_final = join(final_dir, fp)
    if isfile(fp_in_final):
        safe_mkdir(dirname(dst_fp))
        shutil.copy2(fp_in_final, dst_fp)
        return dst_fp

    elif f_by_fp and fp in f_by_fp:
        safe_mkdir(dirname(dst_fp))
        with open(dst_fp, 'wb') as out:
            out.write(f_by_fp[fp].read())
        return dst_fp


def multiqc_prep_data(bcbio_mq_filelist, bcbio_final_dir, new_mq_data_dir,
                      generated_conf, out_filelist_file, out_conf_yaml,
                      additional_files, exclude_files=None, include_files=None,
                      bcbio_mq_yaml=None, new_bcbio_mq_yaml=None):

    if generated_conf:
        with file_transaction(None, out_conf_yaml) as tx:
            with open(tx, 'w') as f:
                yaml.dump(generated_conf, f, default_flow_style=False)

    if bcbio_mq_yaml and new_bcbio_mq_yaml:
        shutil.copy(bcbio_mq_yaml, new_bcbio_mq_yaml)

    ########################################
    ###### Processing bcbio file list ######

    verify_file(bcbio_mq_filelist, is_critical=True)
    qc_files_not_found = []

    # Cromwell?
    cwl_targz = join(dirname(bcbio_mq_filelist), 'multiqc-inputs.tar.gz')
    tar_f_by_fp = dict()
    if isfile(cwl_targz):
        logger.info(f'Found CWL MultiQC output {cwl_targz}, extracting required QC files from the archive')
        if cwl_targz:
            tar = tarfile.open(cwl_targz)
            for member in tar.getmembers():
                rel_fp = member.name
                if 'call-multiqc_summary/execution/qc/multiqc/' in rel_fp:
                    rel_fp = rel_fp.split('call-multiqc_summary/execution/qc/multiqc/')[1]
                tar_f_by_fp[rel_fp] = tar.extractfile(member)

    with file_transaction(None, out_filelist_file) as tx:
        try:
            with open(bcbio_mq_filelist) as inp, open(tx, 'w') as out:
                for fp in [l.strip() for l in inp if l.strip()]:
                    if fp == 'trimmed' or fp.endswith('/trimmed'):
                        continue  # back-compatibility with bcbio
                    if exclude_files:
                        if isinstance(exclude_files, str):
                            exclude_files = [exclude_files]
                        if any(re.search(ptn, fp) for ptn in exclude_files):
                            continue
                    if include_files:
                        if isinstance(include_files, str):
                            include_files = [include_files]
                        if not any(re.search(ptn, fp) for ptn in include_files):
                            continue

                    new_fp = _extract_qc_file(fp, new_mq_data_dir, bcbio_final_dir, tar_f_by_fp)
                    if not new_fp:
                        qc_files_not_found.append(fp)
                        continue
                    else:
                        out.write(new_fp + '\n')

                if qc_files_not_found:
                    warn('-')
                    warn(f'Some QC files from list {bcbio_mq_filelist} were not found:' +
                        ''.join('\n  ' + fpath for fpath in qc_files_not_found))
                for fp in additional_files:
                    if fp:
                        out.write(fp + '\n')
        except OSError as e:
            err(str(e))


def make_multiqc_report(conf_yamls, filelist, multiqc_fpath, extra_params=''):
    postproc_multiqc_dir = safe_mkdir(dirname(multiqc_fpath))

    cmdl = (f'multiqc '
            f'-f{" -v" if logger.is_debug else ""} '
            f'-o {postproc_multiqc_dir} '
            f'-l {filelist}')

    if conf_yamls:
        for conf_yaml in conf_yamls:
            cmdl += ' -c ' + conf_yaml

    if extra_params:
        cmdl += ' ' + extra_params

    if isfile(multiqc_fpath):
        try:
            os.remove(multiqc_fpath)
        except OSError as e:
            err(e)
    run(cmdl, env_vars={'https_proxy': None})
    return multiqc_fpath


def _make_url(fpath, base_dirpath):
    return relpath(fpath, base_dirpath) if verify_file(fpath) else None


def _make_link(fpath, base_dirpath, text=None, blank=False):
    url = _make_url(fpath, base_dirpath)
    if url:
        return '<a href="' + url + '"' + (' target="_blank"' if blank else '') + '>' + (text or basename(fpath)) + '</a>'
    else:
        return '<span>' + (text or basename(fpath)) + '</span>'


def get_run_info(bcbio_proj, base_dirpath, analysis_dir=None,
                 prog_versions_fpath=None, data_versions_fpath=None,
                 new_dir_for_versions=None):
    info('Getting run and codebase information...')
    run_info_dict = dict()
    cur_fpath = abspath(getsourcefile(lambda: 0))
    reporting_suite_dirpath = dirname(dirname(dirname(cur_fpath)))

    run_info_dict["run_date"] = strftime('%d %b %Y, %H:%M (GMT%z)', localtime())

    try:
        from umccrise import _version
    except ImportError:
        err('Cannot import __version__ from umccrise._version and get version')
        umccrsie_version = None
    else:
        umccrsie_version = _version.__version__

    last_modified_datestamp = ''
    try:
        py_fpaths = set()
        for rootpath, dirnames, fnames in os.walk(reporting_suite_dirpath):
            for fn in fnames:
                if fn.endswith('.py'):
                    fpath = abspath(join(rootpath, fn))
                    if isfile(fpath):
                        py_fpaths.add(fpath)
        last_modification_time = max(getmtime(fpath) for fpath in py_fpaths)
    except Exception:
        warn('Cannot get codebase files datestamps')
    else:
        last_modified_datestamp = strftime('%d %b %Y, %H:%M (GMT%z)', localtime(last_modification_time))

    if last_modified_datestamp or umccrsie_version:
        version_text = ''
        if umccrsie_version:
            version_text += umccrsie_version
        if umccrsie_version and last_modified_datestamp:
            version_text += ', '
        if last_modified_datestamp:
            version_text += 'last modified ' + last_modified_datestamp
        run_info_dict['umccrise_version'] = version_text

    # prog versions
    # if prog_versions_fpath and new_dir_for_versions and verify_file(prog_versions_fpath, silent=True):
    info('Adding umccrise versoin into prog_versions file ' + prog_versions_fpath)
    with open(prog_versions_fpath) as f:
        program_versions = dict()
        for l in f:
            l = l.strip()
            if l:
                try:
                    program_versions[l.split(',')[0]] = l.split(',')[1]
                except:
                    pass
    new_prog_versions_fpath = join(new_dir_for_versions, basename(prog_versions_fpath))
    if umccrsie_version:
        program_versions['umccrise'] = umccrsie_version
    with open(new_prog_versions_fpath, 'w') as f:
        for p, v in sorted(program_versions.items(), key=lambda kv: kv[0]):
            f.write(p + ',' + v + '\n')
    info('Saved prog_versions file into ' + new_prog_versions_fpath)
    assert exists(new_prog_versions_fpath)

    programs_url = relpath(new_prog_versions_fpath, base_dirpath)
    run_info_dict['program_versions'] = f'<a href="{programs_url}">program versions</a>'

    # data version
    # if data_versions_fpath and new_dir_for_versions and verify_file(data_versions_fpath, silent=True):
    new_data_versions_fpath = join(new_dir_for_versions, basename(data_versions_fpath).replace(".csv", ".txt"))
    run_simple(f'cp {data_versions_fpath} {new_data_versions_fpath}')
    info('Saved data_versions file into ' + new_data_versions_fpath)
    assert exists(new_data_versions_fpath)
    data_versions_url = relpath(new_data_versions_fpath, base_dirpath)
    run_info_dict['data_versions'] = f'<a href="{data_versions_url}">data versions</a>'

    run_info_dict['analysis_dir'] = analysis_dir or bcbio_proj.final_dir
    return run_info_dict


