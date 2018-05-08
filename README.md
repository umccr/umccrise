UMCCRization of bcbio results: filter, normalise, generate plots and reports
----------------------------------------------------------------------------

[![Build Status](https://travis-ci.org/umccr/umccrise.svg?branch=master)](https://travis-ci.org/umccr/umccrise)

Umccrise post-processess an output from [bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen) somatic variant calling pipeline for cancer samples:

- Filters small somatic calls with [panel of normals](https://github.com/umccr/vcf_stuff#panel-of-normals)
- Filters small germline calls with key genes
- Runs [PCGR](https://github.com/sigven/pcgr) for somatic and germline variants
- Generates an Rmd report with mutational signatures and strand bias analysis
- QCs coverage for 300 key cancer genes
- Filters CNV and plots a diagram
- Filters SV and generates files to view in Ribbon
- Generates mini-bams to view in IGV
- Copies MultiQC reports and summaries from bcbio

Contents:

- [Installation](#installation)
- [Updating](#updating)
- [Testing](#testing)
- [Loading](#loading)
- [Usage](#usage)
    + [Run selected steps](#run-selected-steps)
    + [Run on selected samples](#run-on-selected-samples)
    + [Use HPC cluster](#use-hpc-cluster)
- [Output explanation](#output-explanation)
- [Version history](#version-history)


## Installation

Clone the repository

```
git clone https://github.com/umccr/umccrise
```

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
cat <<EOT > load_umccrise.sh
SCRIPTPATH=\$(dirname \$(readlink -e $(pwd)))
export PATH=\$SCRIPTPATH/miniconda/bin:\$PATH
source activate \$SCRIPTPATH/miniconda/envs/umccrise
EOT
```

Install PCGR

```
# Clone the fork that is decoupled from Docker
git clone https://github.com/vladsaveliev/pcgr

# Download the data bundle
gdrive download 1tOyPmzgXkSZjPJQOojFQxUP8JlQuQqLq  # hg19
gdrive download 1cKq-rgSNCYPCUJ38pCi_xy6_PJH-FZWD  # hg38
gunzip -c pcgr.databundle.grch37.20180422.tgz | tar xvf -
gunzip -c pcgr.databundle.grch37.20180422.tgz | tar xvf -
bash pcgr/install_no_docker/install.sh
```

## Updating

```
source load_umccrise.sh
git pull                                                             # if the code base changed
conda env update -f environment.yml                                  # if dependencies changed
./setup.py develop && source deactivate && source load_umccrise.sh   # if added/renamed packages or scripts
```

## Testing

```
source load_umccrise.sh
git clone https://github.com/umccr/umccrise_test_data
nosetests -s umccrise_test_data/test_umccrise.py
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

## Usage

Runs the patient analysis pipeline on bcbio-nextgen `final` folder.

```
umccrise /path/to/bcbio/project/final -j 30  # run using 30 CPUs
```

The output will be created in `umccrised` folder. To override, use `-o`:

```
umccrise /path/to/bcbio/project/final -o umccrised_results
```

#### Run selected steps

Umccrise workflow consists of the following steps: `pcgr`, `coverage`, `structural`, `small_variants`, `rmd`, `copy_multiqc`, `copy_logs`, `igv`.

To run just a particular step (or steps), use:

```
umccrise /path/to/bcbio/project/final <step_name>
```

Where `<step_name>` is from the list above. E.g.:

```
umccrise /path/to/bcbio/project/final pcgr
```

Note that the `igv` step (preparing minibams and uploading them to `s3://umccr-igv`) takes ~5 hours for a WGS sample compared to ~20 minutes for all other steps combined. For that reason, it is always executed in the end of the pipeline, so you can expect that when it is being executed, all other output is ready.

#### Run on selected samples

By default, Umccrise will process all batches in the run in parallel. You can submit only certain samples/batchs using `--sample` or `--batch` arguments, e.g.:

```
umccrise /path/to/bcbio/project/final --batch cup-batch
umccrise /path/to/bcbio/project/final --sample cup-tumor_1,cup-tumor_2
```

Or you might want to exclude certain samples/batches with `--exclude`:

```
umccrise /path/to/bcbio/project/final --exclude cup-tumor_1,cup-batch_2
```

#### Use HPC cluster

Set `--cluster-auto` option to submit jobs on HPC cluster. Supports Spartan for now.

```
umccrise /path/to/bcbio/project/final -j 30 --cluster-auto
```

Alternatively, you can specify a custom submission template with `--cluster-cmd`, e.g.:

```
umccrise /path/to/bcbio/project/final -j 30 --cluster-cmd "sbatch -p vccc -n {threads} -t 24:00:00 --mem {resources.mem_mb} -J umccrise"
```

Make sure to use `-j` outside of that template: this options tells snakemake how many cores is allowed to use at single moment.


## Output explanation

```
umccrised/
    {batch}/                               # - Folder with a batch {batch} (tumor/normal pair) analysis
        {batch}-{sample}-rmd_report.html   # - Rmd report with mutational signatures, AF frequencies,
                                           #   structural variants, and strand bias plots
        coverage/
            {batch}-{sample}-indexcov/index.html  # - Plots by `goleft indexcov`
            {batch}-{sample}-normal.callable.bed  # - Coverage for exons of 300 AstraZeneca key cancer genes,
            {batch}-{sample}-tumor.callable.bed   #   calculated by `goleft depth`. The "callable" coverage
            {batch}-{sample}-normal.depth.bed     #   thresholds: 10x for normal, 30x for tumor.
            {batch}-{sample}-tumor.depth.bed      #       
        igv/ 
            {batch}-{sample}-normal-mini.bam      # - Minibams, containing locations of AZ 300 cancer genes,
            {batch}-{sample}-tumor-mini.bam       #   and areas around passed variants (SNV, indels, SV, CNV)
        pcgr/
            {batch}-{sample}-somatic.pcgr_acmg.html
            {batch}-{sample}-normal.pcgr_acmg.html
        snv/
            {batch}-{sample}-somatic-ensemble-pon_softfiltered.vcf.gz  # - Somatic small variants (SNV and indels), soft- and
            {batch}-{sample}-somatic-ensemble-pon_hardfiltered.vcf.gz  #   hard-filtered against the panel of normals
            {batch}-{sample}-normal-ensemble-cancer_genes.vcf.gz       # - Germline small variants, subset to 106 cancer genes
        structural/
            {batch}-{sample}-sv-prioritize-manta-pass.bedpe       # - Prioritized Manta SV calls in differrent
            {batch}-{sample}-sv-prioritize-manta-pass.ribbon.bed  #   formats (e.g. to view in Ribbon)
            {batch}-{sample}-sv-prioritize-manta-pass.tsv         #
            {batch}-{sample}-sv-prioritize-manta-pass.vcf         #
            {batch}-{sample}-cnvkit-diagram.pdf                   # - Diagram of CNV called by CNVkit
    log/
        config/                            # - Copy of config folder from the original bcbio-nextgen run
            {project-name}-template.yaml   #
            {project-name}.csv             #
            {project-name}.yaml            #
        {project-name}-data_versions.csv   # - Copy of the bcbio-nextgen file with reference data versions
        {project-name}-programs.txt        # - Copy of the bcbio-nextgen file with software versions
    {project-name}-multiqc_report.html     # - Project-level MultiQC summary report: coverage stats and more
```


## Version history

0.5.0:
- PCGR is deployed directly on Spartan, so no AWS dependency.
  - Add pcgr wrapper: `pcgr variants.vcf.gz cnv.tsv -o results [-g hg38]`
- Correctly providing memory resources on HPC to avoid oom-kill
- On Spartan, support --cluster-auto to automatically substitute proper cluster parameters

0.4.0: 
- Propagate snakemake's cluster options to the wrapper
- Propagate snakemake's `--rerun-incomplete` and to the wrapper
- PCGR: stop using the UID, skip uploading if already uploaded
- Fix `HTTP_PROXY` setting for running on Spartan worker nodes
- New panel of normals: add the large A5 cohort, normalize normals, check only POS for indels
- Correctly set threads to mibibams samtools runs

0.3.0:
- Refactor output folder structure (see docs)
- `pcgr_download` taget to automatically pull results
- Automatically attempt downloading PCGR results after minibams
- Minibams always run in the ned
- Use `--unlock` to automate the snakemake unlock pattern to continue unterrupted runs
- Move `vcf_stuff` to a separate repo
- Remove GL* chromosomes from the cnvkit diagram
- PCGR tomls: changed mutsignatures_normalization to genome
- Rmd: add SV-prioritize table






