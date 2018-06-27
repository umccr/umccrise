import os
import sys
from python_utils.hpc import get_loc
from snakemake.utils import read_job_properties

loc = get_loc()
submit_cmd = loc.submit_job_cmd

timestamp = sys.argv[1]
logs_dir = sys.argv[2]
job_cmd = sys.argv[3]

job_properties = read_job_properties(job_cmd)

submit_cmd = submit_cmd\
    .replace('{threads}', str(job_properties.get('threads', 1)))\
    .replace('{resources.mem_mb}', str(job_properties.get('resources', {}).get('mem_mb', 2000)))

job_name = job_properties.get('rule') or job_properties.get('groupid') or 'umccrise'
job_name += '.' + '__'.join(job_properties['wildcards'].values())

submit_cmd = submit_cmd\
    .replace('{job_name}', job_name)\
    .replace('{log_file}', os.path.join(logs_dir, f'{timestamp}_{job_name}.cluster.log'))

cmd = f'{submit_cmd} {job_cmd}'
sys.stderr.write(cmd + '\n')
os.system(cmd)
