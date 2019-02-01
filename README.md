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
bash miniconda.sh -b -p ./miniconda && rm miniconda.sh
export PATH=$(pwd)/miniconda/bin:$PATH
conda update conda
```

Install environments

```
ENV_NAME=umccrise
conda env create -p $(pwd)/miniconda/envs/${ENV_NAME} --file umccrise/envs/umccrise.yml
conda env create -p $(pwd)/miniconda/envs/${ENV_NAME}_purple --file umccrise/envs/purple.yml
# conda env create -p $(pwd)/miniconda/envs/${ENV_NAME}_pcgr --file umccrise/envs/pcgr_macos.yml    # macos
conda env create -p $(pwd)/miniconda/envs/${ENV_NAME}_pcgr --file umccrise/envs/pcgr_linux.yml    # linux
export PATH=$(pwd)/miniconda/envs/${ENV_NAME}/bin:$PATH
pip install -e umccrise
```

Dirty fix for Raijin [DEPRECATED]

```
ENV_NAME=umccrise
cd ./miniconda/envs/${ENV_NAME}_purple/lib
if [ ! -e libwebp.so.7 ] ; then ln -s libwebp.so.6 libwebp.so.7; fi
cd -
```

To automate sourcing in the future, you can create a loader script

```
ENV_NAME=umccrise
cat <<EOT > load_umccrise.sh
unset PYTHONPATH
unset PERL5LIB
export PATH=$(pwd)/miniconda/envs/${ENV_NAME}/bin:$(pwd)/miniconda/bin:\$PATH
export CONDA_PREFIX=$(pwd)/miniconda/envs/${ENV_NAME}
EOT
```

## Updating

```
source load_umccrise.sh
git pull                                                                  # if the code base changed
conda env update -f envs/umccrise.yml                                     # if dependencies changed
python setup.py develop && source deactivate && source load_umccrise.sh   # if added/renamed packages or scripts
```

## Reference data bundle

Umccrise automatically finds reference data on Spartan and NCI environments, as well as the reference data bundle 
mounted to the Docker image under `/genomes`.

To sync the reference data from Spartan to NCI, use:

```
cd /data/cephfs/punim0010/extras/umccrise
rsync -rv --size-only genomes/ rjn:/g/data3/gx8/extras/umccrise/genomes
```

### Sources for the reference data

#### GNOMAD

Version 2.0.1 (26G, hosted by ensmble, same as in bcbio):

```
wget http://ftp.ensemblorg.ebi.ac.uk/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz
bcftools filter -i 'FILTER="PASS" & (AN_POPMAX>=500 & AF_POPMAX>=0.01 | AN_POPMAX>=100 & AF_POPMAX>=0.01)' gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz -Ob | \
    bcftools annotate -x ID,^INFO/AN_POPMAX,^INFO/AF_POPMAX,FORMAT -Oz -o gnomad_genome.common_pass_clean.vcf.gz
tabix -p vcf gnomad_genome.common_pass_clean.vcf.gz
```

Version 2.1 (latest, 500G, hosted by broad):

```
wget -c https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgzwget -c https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz
bcftools filter -i 'FILTER="PASS" & !segdup & !lcr & !decoy & (AN_popmax>=500 & AF_popmax>=0.01 | AN_popmax>=100 & AF_popmax>=0.01)' gnomad.genomes.r2.1.sites.vcf.bgz -Ob | \
    bcftools annotate -x ID,^INFO/AN_popmax,^INFO/AF_popmax,FORMAT -Oz -o gnomad_genome.r2.1.common_pass_clean.vcf.gz
tabix -p vcf gnomad_genome.r2.1.common_pass_clean.vcf.gz
```

Normalise

```
ref = genomes/GRCh37/GRCh37.fa
export TMPDIR=`pwd`
vt decompose -s gnomad_genome.r2.1.common_pass_clean.vcf.gz | vt normalize -r $ref -n - | gsort -m 3000 /dev/stdin $ref.fai | bgzip -c > gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz
# vcf_norm gnomad_genome.r2.1.common_pass_clean.vcf.gz -o gnomad_genome.r2.1.common_pass_clean.norm_mine.vcf.gz

tabix -f -p vcf gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz
```

#### PCGR

The PCGR data bundle gets refreshed every release, so please select the appropriate one from [PCGR's README](https://github.com/sigven/pcgr#step-2-download-pcgr-and-data-bundle)!

```bash
# Download the data bundles
pip install gdown  # or use  `$(pwd)/miniconda/envs/${ENV_NAME}_pcgr/bin/gdown`
gdown https://drive.google.com/uc?id=<GDOCS_ID_SEE_PCGR_DATABUNDLE_README> -O - | tar xvfz - # hg19
gdown https://drive.google.com/uc?id=<GDOCS_ID_SEE_PCGR_DATABUNDLE_README> -O - | tar xvfz - # hg38

# (Optional) if you are running on AWS, upload the PCGR data bundles to S3 like this:
gdown https://drive.google.com/uc?id=<GDOCS_ID_SEE_PCGR_DATABUNDLE_README> -O - | aws s3 cp - s3://umccr-umccrise-refdata-dev/Hsapiens/GRCh37/PCGR/pcgr.databundle.grch37.YYYMMDD.tgz
gdown https://drive.google.com/uc?id=<GDOCS_ID_SEE_PCGR_DATABUNDLE_README> -O - | aws s3 cp - s3://umccr-umccrise-refdata-dev/Hsapiens/hg38/PCGR/pcgr.databundle.grch38.YYYMMDD.tgz
```

#### SegDup

```
cd problem_regions
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz .
gunzip -c genomicSuperDups.txt.gz | cut -f2,3,4 >> segdup.bed_tmp
gunzip -c genomicSuperDups.txt.gz | cut -f8,9,10 >> segdup.bed_tmp
grep -v gl segdup.bed_tmp | sed 's/chr//' | bedtools sort -i - | bedtools merge -i - > segdup.bed
bgzip -f segdup.bed && tabix -f -p bed segdup.bed.gz
rm segdup.bed_tmp genomicSuperDups.txt.gz
```


## Testing

Tests are stored in a separate repository https://github.com/umccr/umccrise_test_data

```
source load_umccrise.sh
git clone https://github.com/umccr/umccrise_test_data
nosetests -s umccrise_test_data/test.py -a normal
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

Umccrise workflow consists of the following steps: `pcgr`, `coverage`, `structural`, `small_variants`, `rmd`, `multiqc`, `copy_logs`, `igv`.

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

Set `--cluster-auto` (`-c`) option to submit jobs on HPC cluster. Supports Spartan for now.

```
umccrise /path/to/bcbio/project/final -j 30 -c
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

## Dockerized installation

Download conda

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p ./miniconda && rm miniconda.sh
. miniconda/etc/profile.d/conda.sh
```

Create minimal environment

```
conda env create --file deploy/env_docker_wrapper.yml
conda activate umccrise
pip install -e .
```

Pull docker

```
docker pull umccr/umccrise:latest
```

Testing

```
git clone https://github.com/umccr/umccrise_test_data
nosestests -s umccrise_test_data/test.py -a docker
```

Usage

```
umccrise --docker \
    umccrise_test_data/data/bcbio_test_project \
    -o umccrise_test_data/results/dockerized \
    -j 2 \
    --bcbio-genomes umccrise_test_data/data/genomes \
    --pon umccrise_test_data/data/panel_of_normals
```






