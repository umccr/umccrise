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
# soruce ~/reset_path.sh
unset PYTHONPATH
unset CONDA_PREFIX
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
conda activate $(pwd)/miniconda/envs/${ENV_NAME}  # export PATH=$(pwd)/miniconda/envs/${ENV_NAME}/bin:$PATH
pip install -e umccrise
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

## Testing

```
git clone https://github.com/umccr/umccrise_test_data
TEST_OPTS="-c -j2" nosetests -s umccrise_test_data/test.py
```

## Updating

```
source load_umccrise.sh
git pull                                                                  # if the code base changed
conda env update -f envs/umccrise.yml                                     # if dependencies changed
python setup.py develop && source deactivate && source load_umccrise.sh   # if added/renamed packages or scripts
```

## Development

Changes pulled in `umccrise` repository clone folder will affect immidiately due to use of the `-e` option in `pip install -e`. To do the same for other related packages, you can clone them as well (or move already cloned repos from `./umccrise/envs/src`, and run `pip install -e` on them as well:

```
source load_umccrise.sh
mv ./umccrise/envs/src/* .
rm pip-delete-this-directory.txt
pip install -e ngs-utils
pip install -e hpc-utils
pip install -e vcf-stuff
pip install -e multiqc
pip install -e multiqc-bcbio
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

Version 2.1 (latest, 500G, hosted by broad):

```
wget -c https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz |\
    #can optionally remove LCR becaues we annoatate vs LCR futher anyway, but the file will be already small enough:
       #bcftools filter -i 'FILTER="PASS" & segdup=0 & lcr=0 & decoy=0 & (AN_popmax>=500 & AF_popmax>=0.01 | AN_popmax>=100 & AF_popmax>=0.01)' gnomad.genomes.r2.1.sites.vcf.bgz -Ob |\
    bcftools annotate -x ID,^INFO/AN_popmax,^INFO/AF_popmax,FORMAT -Oz -o gnomad_genome.r2.1.common_pass_clean.vcf.gz
tabix -p vcf gnomad_genome.r2.1.common_pass_clean.vcf.gz
```

Normalise (see https://github.com/chapmanb/cloudbiolinux/pull/279, however after all just 5 indels will be changed, so not a bit deal)

```
ref=GRCh37.fa
norm_vcf gnomad_genome.r2.1.common_pass_clean.vcf.gz -o gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz --ref-fasta $ref
```

Counts: 

```
bcftools view gnomad_genome.r2.1.common_pass_clean.vcf.gz -H | wc    # with lcr&segdup&decoy
24671774
bcftools view gnomad_genome.common_pass_clean.vcf.gz -H | wc
21273673
```

Convert to hg38

```
# spartan
CrossMap.py vcf /data/cephfs/punim0010/extras/hg19ToHg38.over.chain.gz ../GRCh37/gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz hg38.fa gnomad_genome.
r2.1.common_pass_clean.norm.vcf.unsorted
bcftools view gnomad_genome.r2.1.common_pass_clean.norm.vcf.unsorted | bcftools sort -Oz -o gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz
tabix -p vcf gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz
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


#### Problem regions

Copy GRCh37 from bcbio

```
cp -r /g/data/gx8/local/development/bcbio/genomes/Hsapiens/GRCh37/coverage/problem_regions problem_regions
```

Generate SegDup

```
cd problem_regions
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz .
gunzip -c genomicSuperDups.txt.gz | cut -f2,3,4 >> segdup.bed_tmp
gunzip -c genomicSuperDups.txt.gz | cut -f8,9,10 >> segdup.bed_tmp
grep -v gl segdup.bed_tmp | sed 's/chr//' | bedtools sort -i - | bedtools merge -i - > segdup.bed
bgzip -f segdup.bed && tabix -f -p bed segdup.bed.gz
rm segdup.bed_tmp genomicSuperDups.txt.gz
```


Lift over to hg38:

```
convert () {
    f=$(basename $1)
    echo "Processing $f"
    zless $1 | py -x "('chr' + x) if not x.startswith('MT') else 'chrM'" | grep -v chrG > $f.hg19
    CrossMap.py bed /g/data3/gx8/extras/hg19ToHg38.over.chain.gz $f.hg19 $f.unsorted
    bedtools sort -i $f.unsorted | bgzip -c > $f
    tabix -p vcf $f
}

convert ../../GRCh37/problem_regions/segdup.bed.gz

mkdir GA4GH ENCODE repeats

cd GA4GH
for fp in $(ls ../../../GRCh37/problem_regions/GA4GH/*.bed.gz) ; do convert $fp ; done

cd ../ENCODE
convert ../../../GRCh37/problem_regions/ENCODE/wgEncodeDacMapabilityConsensusExcludable.bed.gz

cd ../repeats
convert ../../../GRCh37/problem_regions/repeats/LCR.bed.gz
convert ../../../GRCh37/problem_regions/repeats/polyx.bed.gz

cat ../../../GRCh37/problem_regions/repeats/sv_repeat_telomere_centromere.bed | py -x "('chr' + x) if not x.startswith('MT') else 'chrM'" | grep -v chrG > sv_repeat_telomere_centromere.bed_hg19
CrossMap.py bed /g/data3/gx8/extras/hg19ToHg38.over.chain.gz sv_repeat_telomere_centromere.bed_hg19 sv_repeat_telomere_centromere.bed_unsorted
bedtools sort -i sv_repeat_telomere_centromere.bed_unsorted | > sv_repeat_telomere_centromere.bed
```


#### Coding regions (SAGE)

```
cd GRCh37/hmf
generate_bed.py -g GRCh37\
   --principal --key-genes --features CDS | sort -k1,1V -k2,2n | grep -v ^MT | grep -v ^GL\
   | bedtools merge -c 4 -o collapse -i - > coding_regions.bed

cd hg38/hmf   
generate_bed.py -g hg38\
   --principal --key-genes --features CDS | sort -k1,1V -k2,2n | grep -v ^chrM\
   | bedtools merge -c 4 -o collapse -i -\
   > coding_regions.bed
```

#### Ensembl annotation

Using pyensembl package

```
# use ENSEMBL_VERSION=75 for GRCh37, ENSEMBL_VERSION=95 for hg38
export PYENSEMBL_CACHE_DIR=$ENSEMBL_DIR
if [ ! -d $PYENSEMBL_CACHE_DIR/pyensembl ] ; then
    # In 2 steps: first on loging node to make it download the files:
    pyensembl install --release $ENSEMBL_VERSION --species human
    # when it starts `Reading GTF from`, go into a worker node and run again.
fi
```


#### Hotspots

Combining Hartwig's and PCGR hotspots. Stats:
* Hartwigs: 10211 changes in 3650 locations, 
* PCGR: 10627 changes in 2494 locations.
* Overlap: 2960 changes in 968 locations.

The overlap is small, so we better merge sources into a single VCF.

First, download HMF TSV file and convert to VCF:

```
wget https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMF-Pipeline-Resources&files=KnownHotspots.tsv.gz -O KnownHotspots.tsv.gz

echo "##fileformat=VCFv4.2" > hmf.vcf
echo '##INFO=<ID=HMF,Number=.,Type=Flag,Description="Hotspot is from HMF">' >> hmf.vcf
echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> hmf.vcf
gunzip -c KnownHotspots.tsv.gz | py -x "print('\t'.join([x.split()[0], x.split()[1], '.', x.split()[2], x.split()[3], '.', '.', 'HMF']))" >> hmf.vcf
bgzip hmf.vcf
tabix -p vcf hmf.vcf.gz
```

Prepare PCGR hotspots:

```
SRC=/Users/vsaveliev/bio/genomes/pcgr/data/grch37/cancer_hotspots/cancer_hotspots.vcf.gz
bcftools view -h $SRC | grep ^## > pcgr.vcf
echo '##INFO=<ID=PCGR,Number=.,Type=Flag,Description="Hotspot is from PCGR (cancerhotspots.org_v2)">' >> pcgr.vcf
bcftools view -h $SRC | grep ^#CRHOM >> pcgr.vcf
bcftools view -H $SRC | bioawk -t '{ print $1,$2,$3,$4,$5,$6,$7,$8";PCGR" }' >> pcgr.vcf

bgzip pcgr.vcf
tabix -p vcf pcgr.vcf.gz
```

Merge:

```
bcftools merge -m none hmf.vcf.gz pcgr.vcf.gz -Oz -o merged.vcf.gz
tabix -p vcf merged.vcf.gz
```

Convert to hg38

```
cd ../hg38/hotspots
CrossMap.py vcf /g/data3/gx8/extras/hg19ToHg38.over.chain.gz ../../GRCh37/hotspots/merged.vcf.gz ../hg38.fa merged.vcf.gz
```

#### Other HMF files

Download NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz and out_150_hg19.mappability.bed.gz 
from https://resources.hartwigmedicalfoundation.nl

To hg38:

```
convert ../../GRCh37/hmf/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz
convert ../../GRCh37/hmf/out_150_hg19.mappability.bed.gz
```

#### Fusions

We use [HMF fusions](https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd?path=%2FHMF-Pipeline-Resources) for SV prioritization. See `NGS_Utils/ngs_utils/refernece_data/__init__.py` for details.


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






