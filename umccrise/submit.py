#!/usr/bin/env python3
import os
import sys
from python_utils.hpc import get_loc
from snakemake.utils import read_job_properties


submit_cmd = sys.argv[1]
job_cmd = sys.argv[2]

job_properties = read_job_properties(job_cmd)

job_cmd = job_cmd.replace('{threads}', str(job_properties.get('threads', 1)))
job_cmd = job_cmd.replace('{resources.mem_mb}', str(job_properties.get('resources', {}).get('mem_mb', 2000)))
job_cmd = job_cmd.replace('{job_name}', job_properties.get('rule', job_properties.get('groupid', 'umccrise')))

cmd = f'{submit_cmd} {job_cmd}'
sys.stderr.write(cmd + '\n')
os.system(cmd)
