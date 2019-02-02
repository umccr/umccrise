import os
import shutil
import subprocess
from inspect import getsourcefile
from os.path import join, dirname, abspath, pardir, isfile, isdir, islink, basename, relpath, getmtime, exists
from time import localtime, strftime
import yaml
import re

from ngs_utils import logger
from ngs_utils.bcbio import BcbioProject
from ngs_utils.call_process import run
from ngs_utils.file_utils import verify_file, safe_mkdir, can_reuse, file_transaction
from ngs_utils.logger import info, warn, err, critical, timestamp, debug


def make_report_metadata(bcbio_proj, tumor_sample, normal_sample, base_dirpath, analysis_dir=None,
                         program_versions_fpath=None, data_versions_fpath=None):
    conf = dict()
    conf['umccr'] = dict()
    additional_files = []

    conf['umccr']['tumor_name'] = tumor_sample
    conf['umccr']['normal_name'] = normal_sample

    conf['umccr']['is_rnaseq'] = bcbio_proj.is_rnaseq

    # General links
    conf['title'] = bcbio_proj.project_name
    conf['umccr']['run_section'] = get_run_info(bcbio_proj, base_dirpath, analysis_dir=analysis_dir,
                                                program_versions_fpath=program_versions_fpath,
                                                data_versions_fpath=data_versions_fpath)

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


def multiqc_prep_data(bcbio_mq_filelist, bcbio_mq_yaml, bcbio_final_dir,
                      new_mq_data_dir, generated_conf, out_filelist_file, out_conf_yaml,
                      new_bcbio_mq_yaml, additional_files, exclude_files=None):
    if generated_conf:
        with file_transaction(None, out_conf_yaml) as tx:
            with open(tx, 'w') as f:
                yaml.dump(generated_conf, f, default_flow_style=False)

    shutil.copy(bcbio_mq_yaml, new_bcbio_mq_yaml)

    verify_file(bcbio_mq_filelist, is_critical=True)
    qc_files_not_found = []

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
                    old_fpath = join(bcbio_final_dir, fp)
                    new_fpath = join(new_mq_data_dir, fp)
                    if not verify_file(old_fpath, silent=True):
                        fp = fp.replace('inputs/', '').replace('report/metrics/', 'project/multiqc/')  # CWL?
                        old_fpath = join(bcbio_final_dir, fp)
                        if not verify_file(old_fpath, silent=True):
                            qc_files_not_found.append(old_fpath)
                            continue
                    safe_mkdir(dirname(new_fpath))
                    shutil.copy2(old_fpath, new_fpath)
                    out.write(new_fpath + '\n')

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
                 program_versions_fpath=None, data_versions_fpath=None):
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

    if verify_file(program_versions_fpath, silent=True):
        with open(program_versions_fpath) as f:
            program_versions = dict(l.strip().split(',')[:2] for l in f.readlines())
        if umccrsie_version:
            program_versions['umccrise'] = umccrsie_version
        try:
            with open(program_versions_fpath, 'w') as f:
                for p, v in sorted(program_versions.items(), key=lambda kv: kv[0]):
                    f.write(p + ',' + v + '\n')
        except OSError as e:
            err(e)
        programs_url = relpath(program_versions_fpath, base_dirpath)
        run_info_dict['program_versions'] = '<a href="{programs_url}">program versions</a>'.format(**locals())

    if verify_file(data_versions_fpath, silent=True):
        datas_url = relpath(data_versions_fpath, base_dirpath)
        run_info_dict['data_versions'] = '<a href="{datas_url}">data versions</a>'.format(**locals())

    run_info_dict['analysis_dir'] = analysis_dir or bcbio_proj.final_dir
    return run_info_dict


