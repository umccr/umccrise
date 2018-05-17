import os
import sys
from python_utils.hpc import get_loc
from snakemake.utils import read_job_properties

loc = get_loc()
submit_cmd = loc.submit_job_cmd

job_cmd = sys.argv[1]

job_properties = read_job_properties(job_cmd)

submit_cmd = submit_cmd.replace('{threads}', str(job_properties.get('threads', 1)))
submit_cmd = submit_cmd.replace('{resources.mem_mb}', str(job_properties.get('resources', {}).get('mem_mb', 2000)))
job_name = job_properties.get('rule') or job_properties.get('groupid') or 'umccrise'
submit_cmd = submit_cmd.replace('{job_name}', job_name)

cmd = f'{submit_cmd} {job_cmd}'
sys.stderr.write(cmd + '\n')
os.system(cmd)
