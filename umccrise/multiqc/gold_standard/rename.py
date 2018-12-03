import os
from os.path import join
from ngs_utils.file_utils import safe_mkdir
import shutil
import sys

final_dir = sys.argv[1]
rn_final_dir = sys.argv[2]

d = {
    "PRJ180506_E194-T01-D": "Alice_T",
    "PRJ180507_E194-B01-D": "Alice_B",
    "PRJ180492_E201-T01-D": "Bob_T",
    "PRJ180493_E201-B01-D": "Bob_B",
    "PRJ180494_E199-T01-D": "Chen_T",
    "PRJ180495_E199-B01-D": "Chen_B",
    "PRJ180253_E190-T01-D": "Dakota_T",
    "PRJ180254_E190-B01-D": "Dakota_B",
    "PRJ180499_E202-T01-D": "Eugene_T",
    "PRJ180500_E202-B01-D": "Eugene_B",
}

def rn(s):
    for k, v in d.items():
        s = s.replace(k, v)
    return s

safe_mkdir(rn_final_dir)

for root, dirs, files in os.walk(final_dir):
    rn_root = rn(root.replace(final_dir, rn_final_dir))
    for dname in dirs:
        dpath = os.path.join(root, dname)
        rn_dpath = join(rn_root, rn(dname))
        safe_mkdir(rn_dpath)

    for fname in files:
        fpath = os.path.join(root, fname)
        rn_fpath = join(rn_root, rn(fname))
        if not rn_fpath.endswith('.zip'):
            cont = open(fpath).read()
            rn_cont = rn(cont)
            with open(rn_fpath, 'w') as out:
                out.write(rn_cont)
        else:
            shutil.copy(fpath, rn_fpath)

