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
from ngs_utils.Sample import BaseProject, BaseBatch
from ngs_utils.call_process import run_simple
from ngs_utils.file_utils import verify_file, safe_mkdir, can_reuse, file_transaction
from ngs_utils.logger import info, warn, err, critical, timestamp, debug


def make_report_metadata(proj: BaseProject, batch: BaseBatch=None, base_dirpath=None,
                         prog_versions_fpath=None, data_versions_fpath=None,
                         new_dir_for_versions=None):
    conf = dict()
    conf['umccr'] = dict()
    additional_files = []

    if batch:
        conf['umccr']['tumor_name']         = batch.tumors[0].name
        conf['umccr']['normal_name']        = batch.normals[0].name
        conf['umccr']['tumor_rgid']         = batch.tumors[0].rgid
        conf['umccr']['normal_rgid']        = batch.normals[0].rgid
        conf['umccr']['batch']              = batch.name

    conf['umccr']['is_rnaseq'] = proj.is_rnaseq

    # General links
    conf['title'] = proj.project_name
    conf['umccr']['run_section'] = get_run_info(
        proj,
        base_dirpath=base_dirpath,
        prog_versions_fpath=prog_versions_fpath,
        data_versions_fpath=data_versions_fpath,
        new_dir_for_versions=new_dir_for_versions
    )

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


def multiqc_prep_data(generated_conf, out_filelist_file, out_conf_yaml, qc_files):
    if generated_conf:
        with file_transaction(None, out_conf_yaml) as tx:
            with open(tx, 'w') as f:
                yaml.dump(generated_conf, f, default_flow_style=False, width=float("inf"))

    with file_transaction(None, out_filelist_file) as tx:
        try:
            with open(tx, 'w') as out:
                for fp in qc_files:
                    if fp:
                        out.write(fp + '\n')
        except OSError as e:
            err(str(e))


def make_multiqc_report(conf_yamls, filelist, multiqc_fpath, extra_params=''):
    postproc_multiqc_dir = safe_mkdir(dirname(multiqc_fpath))

    cmdl = (f'unset https_proxy; multiqc '
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
    run_simple(cmdl)
    return multiqc_fpath


def _make_url(fpath, base_dirpath):
    return relpath(fpath, base_dirpath) if verify_file(fpath) else None


def _make_link(fpath, base_dirpath, text=None, blank=False):
    url = _make_url(fpath, base_dirpath)
    if url:
        return '<a href="' + url + '"' + (' target="_blank"' if blank else '') + '>' + (text or basename(fpath)) + '</a>'
    else:
        return '<span>' + (text or basename(fpath)) + '</span>'


def get_run_info(proj: BaseProject, base_dirpath=None,
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
    if base_dirpath:
        if prog_versions_fpath:
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
            safe_mkdir(new_dir_for_versions)
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
        if data_versions_fpath:
            safe_mkdir(new_dir_for_versions)
            new_data_versions_fpath = join(new_dir_for_versions, basename(data_versions_fpath).replace(".csv", ".txt"))
            run_simple(f'cp {data_versions_fpath} {new_data_versions_fpath}')
            info('Saved data_versions file into ' + new_data_versions_fpath)
            assert exists(new_data_versions_fpath)
            data_versions_url = relpath(new_data_versions_fpath, base_dirpath)
            run_info_dict['data_versions'] = f'<a href="{data_versions_url}">data versions</a>'

    return run_info_dict


