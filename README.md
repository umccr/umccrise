UMCCRization of bcbio results: filter, normalise, generate plots and reports
----------------------------------------------------------------------------

[![Build Status](https://travis-ci.org/umccr/umccrise.svg?branch=master)](https://travis-ci.org/umccr/umccrise)

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

To update
```
source load_umccrise.sh
git pull                                                             # if the code base changed
conda env update -f environment.yml                                  # if dependencies changed
./setup.py develop && source deactivate && source load_umccrise.sh   # if added/renamed packages or scripts
```

To test
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

To just run a particular part of the workflow, use:
```
umccrise /path/to/bcbio/project/final <part_name>
```
Where `<part_name>` is one of `pcgr`, `pcgr_download`, `coverage`, `structural`, `small_variants`, `igv`, `sig`, `copy_multiqc`, `copy_logs`.

E.g.:
```
umccrise /path/to/bcbio/project/final pcgr
```

Umccrise submits a request to PCGR AWS instance. To download the results back, use the `pcgr_download` target, and specify the unique ID from the original umccrise run when the PCGR request was submitted:
```
umccrise /path/to/bcbio/project/final pcgr_download --uid f725ab
```

## Output

```
umccrised/
    {batch}/                               # - Folder with a batch {batch} (tumor/normal pair) analysis
        {batch}-rmd_report.html            # - Rmd report with mutational signatures, AF frequencies,
                                           #   structural variants, and strand bias plots
        coverage/
            {batch}-indexcov/index.html    # - Plots by `goleft indexcov`
            {batch}-normal.callable.bed    # - Coverage for exons of 300 AstraZeneca key cancer genes,
            {batch}-tumor.callable.bed     #   calculated by `goleft depth`. The "callable" coverage
            {batch}-normal.depth.bed       #   thresholds: 10x for normal, 30x for tumor.
            {batch}-tumor.depth.bed        #       
        igv/ 
            {batch}-normal_mini.bam        # - Minibams, containing locations of AZ 300 cancer genes,
            {batch}-tumor_mini.bam         #   and areas around passed variants (SNV, indels, SV, CNV)
        pcgr/
            {batch}-5a6808-normal.tar.gz   # - Tarballs ready to be processed by PCGR. `5a6808` is
            {batch}-5a6808-somatic.tar.gz  #   a unique ID of a umccrise run.
        snv/
            {batch}-somatic-ensemble-pon_softfiltered.vcf.gz  # - Somatic small variants (SNV and indels), soft- and
            {batch}-somatic-ensemble-pon_hardfiltered.vcf.gz  #   hard-filtered against the panel of normals
            {batch}-normal-ensemble-cancer_genes.vcf.gz       # - Germline small variants, subset to 106 cancer genes
        structural/
            {batch}-sv-prioritize-manta-pass.bedpe       # - Prioritized Manta SV calls in differrent
            {batch}-sv-prioritize-manta-pass.ribbon.bed  #   formats (e.g. to view in Ribbon)
            {batch}-sv-prioritize-manta-pass.tsv         #
            {batch}-sv-prioritize-manta-pass.vcf         #
            {batch}-cnvkit-diagram.pdf                   # - Diagram of CNV called by CNVkit
    log/
        config/                            # - Copy of config folder from the original bcbio-nextgen run
            {project-name}-template.yaml   #
            {project-name}.csv             #
            {project-name}.yaml            #
        {project-name}-data_versions.csv   # - Copy of the bcbio-nextgen file with reference data versions
        {project-name}-programs.txt        # - Copy of the bcbio-nextgen file with software versions
    {project-name}-multiqc_report.html     # - Project-level MultiQC summary report: coverage stats and more
```
