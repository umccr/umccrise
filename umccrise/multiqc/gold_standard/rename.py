import os
from os.path import join
from ngs_utils.file_utils import safe_mkdir
import shutil
import sys

inp_dir = sys.argv[1]
out_dir = sys.argv[2]
out_dir_hg38 = sys.argv[3]

sample_mapping = {
    "E194__PRJ180506_E194-T01-D": "Alice_T",
    "E194__PRJ180507_E194-B01-D": "Alice_B",
    "E201__PRJ180492_E201-T01-D": "Bob_T",
    "E201__PRJ180493_E201-B01-D": "Bob_B",
    "E199__PRJ180494_E199-T01-D": "Chen_T",
    "E199__PRJ180495_E199-B01-D": "Chen_B",
    "E190__PRJ180253_E190-T01-D": "Dakota_T",
    "E190__PRJ180254_E190-B01-D": "Dakota_B",
    "E202__PRJ180499_E202-T01-D": "Elon_T",
    "E202__PRJ180500_E202-B01-D": "Elon_B",
    "PRJ180506_E194-T01-D": "Alice_T",
    "PRJ180507_E194-B01-D": "Alice_B",
    "PRJ180492_E201-T01-D": "Bob_T",
    "PRJ180493_E201-B01-D": "Bob_B",
    "PRJ180494_E199-T01-D": "Chen_T",
    "PRJ180495_E199-B01-D": "Chen_B",
    "PRJ180253_E190-T01-D": "Dakota_T",
    "PRJ180254_E190-B01-D": "Dakota_B",
    "PRJ180499_E202-T01-D": "Elon_T",
    "PRJ180500_E202-B01-D": "Elon_B",
}

def replace_sample_names_in_file_content(cont: str):
    for k, v in sample_mapping.items():
        cont = cont.replace(k, v)
    return cont

def replace_chrom_names_in_file_content(cont: str):
    new_cont = ''
    for s in cont.splitlines():
        # fixing chromosome names in -idxstats.txt or .mosdepth.region.dist.txt
        if any(s.startswith(str(chrom_name)) for chrom_name in list(range(1, 23)) + ['X', 'Y']):
            s = 'chr' + s
        elif s.startswith('MT'):
            s = 'chrM' + s[2:]
        elif s.startswith('GL'):
            s = None
        if s is not None:
            new_cont += s + '\n'
    return new_cont


safe_mkdir(out_dir)
safe_mkdir(out_dir_hg38)

for root, dirs, files in os.walk(inp_dir):
    rn_root = replace_sample_names_in_file_content(root.replace(inp_dir, out_dir))
    rn_root_hg38 = replace_sample_names_in_file_content(root.replace(inp_dir, out_dir_hg38))
    for dname in dirs:
        dpath = os.path.join(root, dname)
        rn_dpath = join(rn_root, replace_sample_names_in_file_content(dname))
        rn_dpath_hg38 = join(rn_root_hg38, replace_sample_names_in_file_content(dname))
        safe_mkdir(rn_dpath)
        safe_mkdir(rn_dpath_hg38)

    for fname in files:
        fpath = os.path.join(root, fname)
        rn_fpath = join(rn_root, replace_sample_names_in_file_content(fname))
        rn_fpath_hg38 = join(rn_root_hg38, replace_sample_names_in_file_content(fname))

        if not fpath.endswith('.zip'):
            cont = open(fpath).read()
            rn_cont = replace_sample_names_in_file_content(cont)
            with open(rn_fpath, 'w') as out:
                out.write(rn_cont)

            if rn_fpath.endswith('-idxstats.txt') or\
                    rn_fpath.endswith('.mosdepth.region.dist.txt') or\
                    rn_fpath.endswith('-indexcov.roc'):
                rn_cont = replace_chrom_names_in_file_content(rn_cont)
            with open(rn_fpath_hg38, 'w') as out:
                out.write(rn_cont)

        else:
            shutil.copy(fpath, rn_fpath)
            shutil.copy(fpath, rn_fpath_hg38)

