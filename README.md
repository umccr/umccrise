UMCCRisation of Bcbio results. Filter, plot, put together, report
-----------------------------------------------------------------

## Installation
Install conda
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ./miniconda
export PATH=$(pwd)/miniconda/bin:$PATH
```

Install umccrise
```
conda env create -p $(pwd)/miniconda/envs/umccrise --file environment.yml
source activate $(pwd)/miniconda/envs/umccrise
pip install -e .
```

Create a loader script
```
echo "export PATH=$(pwd)/miniconda/bin:\$PATH" > load_umccrise.sh
echo "source activate $(pwd)/miniconda/envs/umccrise" >> load_umccrise.sh
```

## Loading

*Raijin:*
```
source /g/data3/gx8/extras/umccrise/load_umccrise.sh
```

*Spartan:*
```
source /data/cephfs/punim0010/extras/umccrise/load_umccrise.sh
```

## Patient analysis
```
cd /path/to/bcbio/project/final
umccrise . -j 30  # run using 30 CPUs
```
The output will be created in `umccrised` folder.

To just run a particular part of the workflow, use:
```
umccrise . <part_name>
```
Where `<part_name>` is one of `pcgr`, `coverage`, `structural`, `igv`, `sig`, `symlink_multiqc`, `copy_logs`.

E.g.:
```
umccrise . structural
```


## Somatic filtering

This repository provides description, code and results for the approaches to somatic variant filtering in UMCCR.

We summarize the filtering ideas into a [spreadsheet](https://docs.google.com/spreadsheets/d/1Xbz4nW76mofKb9ym3C725W035qkA7JuUu8_FvYSCOT0/edit#gid=0), where for problem and idea we provide a corresponding solution (if available) used in [CaVEMan](https://github.com/cancerit/cgpCaVEManWrapper), [bcbio](http://bcb.io/2015/03/05/cancerval/), and [AZ VarDict pipeline](vardict/vardict_filtering.md). The basic factors and ideas are from Matt Eldridge's [slides](https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2017/Day3/somatic_snv_filtering.html).

### Panel of normals

Removing variants detected as germline in a set of unrelated normal tissue samples helps to reduce the FP rate when it was caused by unbalanced coverage in matching regions normals.

On the image, evaluation of the ICGC MB T/N variant calling with 300x tumor coverage, and 50 normal coverage. The number in `vardict_n1` means how many heats in PoN we allow before we filter out the variant.

![Evaluation of the ICGC MB T/N variant calling with 300x tumor coverage, and 50 normal coverage](benchmark_50_300.png)

Usage:

1. Create a `config.yaml` file. For instance, with the following config, the pipeline will pick all `.vcf.gz` files in `/path/to/normals` as a panel of normals, then it will apply filtering to 2 samples 3 times, once for each threshold 1, 2, and 3. You can then evaluate which threshold suits you most. Output files will use the names `syn3-vardict_n1.vcf.gz`, `syn3-vardict_n2.vcf.gz`, etc.
```
samples:
    ensemble:    MB_100vs50-ensemble-annotated.vcf.gz
    mutect2:     MB_100vs50-mutect2-annotated.vcf.gz
    vardict:     MB_100vs50-vardict-annotated.vcf.gz
    strelka2:    MB_100vs50-strelka2-annotated.vcf.gz

normals_dir: /path/to/normals

hits_thresholds: [1, 2, 3]
```

2. Then run the pipeline with:
```
snakemake -s /somatic_filtering/panel_of_normals/Snakefile --configfile config.yaml
```
The filtered VCF files will be written to `work/pon/filter/`.

### Validation

To generate the table as above for the PoN, use the following pipeline. 

1. Create `config.yaml` like the following:
```
samples:
    ensemble:    MB_100vs50-ensemble-annotated.vcf.gz
    ensemble_n1: pon/work/filter/ensemble-ann-n1.vcf.gz
    ensemble_n2: pon/work/filter/ensemble-ann-n2.vcf.gz
    ensemble_n3: pon/work/filter/ensemble-ann-n3.vcf.gz
    mutect2:     MB_100vs50-mutect2-annotated.vcf.gz
    mutect2_n1:  pon/work/filter/mutect2-ann-n1.vcf.gz
    mutect2_n2:  pon/work/filter/mutect2-ann-n2.vcf.gz
    mutect2_n3:  pon/work/filter/mutect2-ann-n3.vcf.gz
    vardict:     MB_100vs50-vardict-annotated.vcf.gz
    vardict_n1:  pon/work/filter/vardict-ann-n1.vcf.gz
    vardict_n2:  pon/work/filter/vardict-ann-n2.vcf.gz
    vardict_n3:  pon/work/filter/vardict-ann-n3.vcf.gz
    strelka2:    MB_100vs50-strelka2-annotated.vcf.gz
    strelka2_n1: pon/work/filter/strelka-ann-n1.vcf.gz
    strelka2_n2: pon/work/filter/strelka-ann-n2.vcf.gz
    strelka2_n3: pon/work/filter/strelka-ann-n3.vcf.gz

#regions: /Users/vsaveliev/Analysis/snv_validation/dream_syn3_grch37/NGv3.bed

truth_variants:  /data/cephfs/punim0010/data/External/Reference/ICGC_MB/MB-benchmark.vcf.gz
#truth_regions:  truth_regions.bed
reference_fasta: GRCh37/seq/GRCh37.fa
```

`regions` and `truth_regions` are the optional fields for the validation target BED files. The BED files, if any provided, will be merged together and used to subset the variants in both truth and query variant sets.

`reference_fasta` can be an absolute path, or a relative path in Spartan or Raijin servers.

2. Run the pipeline with:
```
snakemake -s /validation/Snakefile --configfile config.yaml
```
The result will be written to a tab-separated file report.tsv, which can be inserted into Excel for highliting.

### Normalization

Validation pipeline uses the following variant normalization approach:

1. Split multi-allelic variants into single sample records.

For instance, split one record 
```
#CHROM  POS     ID      REF     ALT
1       10       .      A       T,C
```
Into 2 separate records
```
#CHROM  POS     ID      REF     ALT
1       10       .      A       T
1       10       .      A       C
```
For that, we are using [vt tools](https://github.com/atks/vt):
```
vt decompose -s vcf_file
```

2. Decompose biallelic block substitutions. 

For instance, split the following one records:

```
#CHROM  POS     ID      REF     ALT
1       20      .       AG      CT       
```
into 2 separate ones:
```
#CHROM  POS     ID      REF     ALT
1       20       .      A       G
1       20       .      G       T
```

We are using for that vcflib's [`vcfallelicprimitives`](https://github.com/vcflib/vcflib#vcfallelicprimitives):

```
vcfallelicprimitives -t DECOMPOSED --keep-geno vcf_file
```

3. Left-align and normalize indels, check if REF alleles match the reference.

For instance, given that the reference chromosome 1 starts with `GCTCCG`, split the following records
```
#CHROM  POS     ID      REF     ALT
1       2       .       CTCC    CCC,C,CCCC
```
into the following 3:
```
#CHROM  POS     ID      REF     ALT
1       1       .       GCTC    G
1       2       .       CT      C
1       3       .       T       C
```

### VarDict VCF filtering

Commands to filter VarDict VCF files.

Moves main sample AF and DP fields into INFO, for PCGR post-processing:
```
proc_vardict_vcf fix_fields vardict.vcf.gz > out.vcf
```

Applies special AF threshold filtering to homopolimers based on `MSI` `INFO` fields generated by VarDict. Writes `MSI_FILTER` into the `FILTER` field.
```
proc_vardict_vcf filter_low_af_msi vardict.vcf.gz > out.vcf
```








