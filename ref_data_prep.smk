from os.path import join, splitext, isfile
from pybedtools import BedTool
from ngs_utils.key_genes_utils import get_genes_from_file
from ngs_utils.logger import info
from ngs_utils import gtf
from python_utils.hpc import get_ref_file
from umccrise import package_path
from umccrise.reference_data import get_key_genes_fpath_300


GENOMES = ['GRCh37', 'hg38']
ref_data_dir = 'ref_data'


rule all:
    input:
        expand(join(ref_data_dir, 'key_genes.{genome}.bed'), genome=GENOMES)


rule extract_transcripts:
    input:
        az_gene_file = get_key_genes_fpath_300(),
        umccr_gene_file = 'umccr_key_genes.txt'
    output:
        temp(join(ref_data_dir, 'key_genes.{genome}.transcripts.unmerged.bed'))
    run:
        info(f'Getting coordinates for key/target genes for {wildcards.genome}')
        gene_names = get_genes_from_file(input.az_gene_file) | get_genes_from_file(input.umccr_gene_file)

        gtf_fpath = get_ref_file(wildcards.genome, key='gtf')
        info('Reading the GTF database')
        db = gtf.get_gtf_db(gtf_fpath)

        def _get(_rec, _key):
            val = _rec.attributes.get(_key)
            if val is None:
                return None
            assert len(val) == 1, (_key, str(val))
            return val[0]

        lines = []
        info(f'Extracting transcripts for {len(gene_names)} genes')
        for rec in db.all_features(order_by=('seqid', 'start', 'end')):
            gname = _get(rec, 'gene_name')
            if gname not in gene_names: continue

            if rec.featuretype != 'transcript': continue
            if rec.end - rec.start < 0: continue

            fs = [rec.chrom,
                  str(rec.start - 1),
                  str(rec.end),
                  gname]
            lines.append(fs)
        BedTool(lines).saveas(output[0])


rule merge_transcripts_bed:
    input:
        rules.extract_transcripts.output[0]
    output:
        join(ref_data_dir, 'key_genes.{genome}.bed')
    run:
        info('Sorting and merging results')
        BedTool(input[0]).sort().merge(c=4, o='distinct').saveas(output[0])
