**UMCCR WGS tumor/normal reporting**

umccrise is a [Snakemake](https://github.com/snakemake/snakemake) workflow that
post-processes results from the
[Illumina DRAGEN](https://sapac.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html)
WGS tumor/normal pipeline and generates HTML reports helpful for researchers and
curators at UMCCR.

- [Summary](#summary)
- [Detailed Workflow](#detailed-workflow)
- [History](#history)
- [Example Reports](#example-reports)
- [Usage](#usage)
- [Installation](#installation)
- [Reference data](#reference-data)
  - [Versioning](#versioning)
  - [Syncing with Spartan](#syncing-with-spartan)
  - [Syncing with AWS S3](#syncing-with-aws-s3)
- [Testing](#testing)
- [AWS](#aws)
- [HPC (NCI Gadi)](#hpc-nci-gadi)
- [Advanced usage](#advanced-usage)
  - [Inputs with named arguments](#inputs-with-named-arguments)
  - [Controlling the number of CPUs](#controlling-the-number-of-cpus)
  - [Running selected stages](#running-selected-stages)
  - [Custom input](#custom-input)
  - [Running on selected samples](#running-on-selected-samples)
- [Updating](#updating)
- [Development](#development)
- [Docker](#docker)
- [Building reference data](#building-reference-data)
    - [PURPLE](#purple)
    - [GNOMAD](#gnomad)
    - [PCGR](#pcgr)
    - [Problem regions](#problem-regions)
    - [Coding regions (SAGE)](#coding-regions-sage)
    - [Ensembl annotation](#ensembl-annotation)
    - [Hotspots](#hotspots)
    - [Other HMF files](#other-hmf-files)
    - [Fusions](#fusions)
    - [SnpEff](#snpeff)
  - [DVC](#dvc)
- [GRIDSS and LINX](#gridss-and-linx)

## Summary

In summary, umccrise can:

- Filter artefacts and germline leakage from somatic variant calls
- Run [PCGR](https://github.com/sigven/pcgr) to annotate, prioritize and report
  somatic variants
- Run [CPSR](https://github.com/sigven/cpsr) to annotate, prioritize and report
  germline variants
- Filter, annotate, prioritize and report structural variants (SVs) from
  [Manta](https://github.com/Illumina/manta)
- Run [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) to
  call copy number variants (CNVs), recover SVs, and infer tumor purity & ploidy
- Generate a [MultiQC](https://github.com/ewels/MultiQC) report that summarizes
  quality control statistics in context of background "gold standard" samples
- Generate a cancer report with mutational signatures, inferred HRD status,
  circos plots, prioritized copy number and structural variant calls
- Run [CACAO](https://github.com/sigven/cacao) to calculate coverage in common
  hotspots, as well as
  [goleft indexcov](https://github.com/brentp/goleft/tree/master/indexcov) to
  estimate coverage problems
- Run [Conpair](https://github.com/nygenome/Conpair) to estimate tumor/normal
  concordance and sample contamination
- Run [oviraptor](https://github.com/umccr/oviraptor) to detect viral
  integration sites and affected genes

## Detailed Workflow

See [workflow.md](workflow.md) for a detailed description of the workflow.

## History

See [HISTORY.md](HISTORY.md) for the version history.

## Example Reports

Below are example reports for a HCC1395/HCC1395BL cell line tumor/normal pair
[sequenced and validated by the SEQC-II consortium](https://sites.google.com/view/seqc2/home/data-analysis/high-confidence-somatic-snv-and-indel-v1-1).

[1. MultiQC (quality control metrics and plots) ![MultiQC](docs/multiqc.png)](https://umccr.github.io/umccrise/SEQC-II-50pc__SEQC-II_Tumor_50pc-multiqc_report.html)

<br>

[2. Cancer report (mutational signatures, circos plots, CNV, SV, oncoviruses) ![Cancer report](docs/cancer_report.png)](https://umccr.github.io/umccrise/SEQC-II-50pc__SEQC-II_Tumor_50pc_cancer_report.html)

<br>

[3. CPSR (germline variants) ![CPSR](docs/cpsr.png)](https://umccr.github.io/umccrise/SEQC-II-50pc__SEQC-II_Tumor_50pc-normal.cpsr.html)

<br>

[4. PCGR (somatic variants) ![PCGR](docs/pcgr.png)](https://umccr.github.io/umccrise/SEQC-II-50pc__SEQC-II_Tumor_50pc-somatic.pcgr.html)

<br>

[5. CACAO (coverage reports) ![CACAO](docs/cacao_tumor.png)](https://umccr.github.io/umccrise/SEQC-II-50pc__SEQC-II_Tumor_50pc-tumor.cacao.html)

<br>

## Usage

Given input data from DRAGEN `somatic` and `germline` output folders, or a
custom set of BAM or VCF files, umccrise can be run with:

```
umccrise <input-data ...> -o umccrised
```

For more options, see [Advanced usage](#advanced-usage).

## Installation

Create a `umccrise` directory and install the umccrise GitHub repo along with
the required [conda](https://docs.conda.io/projects/conda/en/latest/index.html)
environments with the following:

```shell
mkdir umccrise
cd umccrise
git clone https://github.com/umccr/umccrise umccrise.git
bash umccrise.git/install.sh
```

The above will generate a `load_umccrise.sh` script that can be sourced to load
the umccrise conda environment on demand:

```shell
source load_umccrise.sh
```

## Reference data

umccrise needs a 64G bundle of reference data to run. From within the UMCCR AWS setup, sign
in to AWS, and run `umccrise_refdata_pull`

```shell
aws sso login --profile sso-dev-admin
umccrise_refdata_pull
export UMCCRISE_GENOMES=${PWD}/refdata/genomes
```

Alternatively, you can specify a custom path with `--genomes <path>`. The path
can be a tarball and will be automatically extracted.

The path can also be a location on S3 or GDS, prefixed with `s3://` or `gds://`.
E.g.:

```shell
umccrise /input --genomes s3://umccr-refdata-dev/genomes
```

Versioned locations would also be checked. For the case above, umccrise will
check the following locations in the order specified:

- `s3://umccr-refdata-dev/genomes_102`
- `s3://umccr-refdata-dev/genomes_10`, and
- `s3://umccr-refdata-dev/genomes`, assuming that the
  [reference_data](https://github.com/umccr/reference_data) package version is
  `1.0.2`.

umccrise will sync the reference data locally into a `~/umccrise_genomes`
directory. You can symlink any other path to that path if you want a different
location. If the data is already downloaded, umccrise will only attempt to
update the changed files or upload new ones. To avoid attempts to check S3/GDS
again at all, specify the downloaded location directly:
`--genomes ~/umccrise_genomes`

Another option to specify the reference data is through an environment variable
`$UMCCRISE_GENOMES`

If you have access to UMCCR's AWS account, you can sync the reference data from
`s3://umccr-refdata-dev`. If you have access to UMCCR's NCI Gadi account, you
can sync the data from `/g/data3/gx8/extras/umccrise/genomes`. Otherwise, you
can build the bundle from scratch following the
[details below](#building-reference-data).

### Versioning

The reference data is versioned as a python package at
<https://github.com/umccr/reference_data>

### Syncing with AWS S3

```shell
ref_data_version=1.0.0
aws s3 sync hg38 s3://umccr-refdata-dev/genomes_${ref_data_version//./}/hg38
aws s3 sync hg38-manifest.txt s3://umccr-refdata-dev/genomes_${ref_data_version//./}/hg38-manifest.txt
```

## Testing

Load the umccrise environment, clone the repo with toy test data, and run
nosetests:

```shell
source load_umccrise.sh
git clone https://github.com/umccr/umccrise_test_data
TEST_OPTS="-c -j2" nosetests -s umccrise_test_data/test.py
```

## AWS

umccrise on AWS is run via AWS Batch in a defined compute environment. This is
set up and maintained via the
[umccrise Terraform Stack](https://github.com/umccr/infrastructure/tree/master/terraform/stacks/umccrise).
This stack also defines the version of umccrise that is used within AWS and how
umccrise jobs are triggered.

## Advanced usage

### Inputs with named arguments

Inputs can be provided to umccrise as a positional argument (see
[Usage](#usage)) or alternatively as named arguments (see examples below). This
is useful when dealing with DRAGEN input, which have two paired input
directories (somatic and germline). The patient and sample identifiers can also
be explicitly set for DRAGEN data - in some instances this is required as these
identifiers cannot be automatically inferred.

```bash
# DRAGEN input with named arguments
umccrise --dragen_somatic_dir PATH --dragen_germline_dir PATH -o umccrised/

# Explicitly setting subject identifier for provided DRAGEN input
umccrise --dragen_somatic_dir PATH --dragen_germline_dir PATH --dragen_subject_id IDENTIFIER -o umccrised/
```

### Controlling the number of CPUs

To set the number of allowed CPUs to use, set the `-j` option:

```
umccrise <input-folder> -j30
```

### Running selected stages

The umccrise workflow includes multiple processing stages, that can optionally
be run in isolation. The following stages are run by default:

- `conpair`
- `structural`
- `somatic`, `germline` (part of `small_variants`)
- `pcgr`
- `cpsr`
- `purple`
- `mosdepth`, `goleft`, `cacao` (part of `coverage`)
- `oncoviruses`
- `cancer_report`
- `multiqc`

The following stages are optionally available and can be enabled with `-T`:

- `microbiome`
- `immuno`

Example:

```shell
# Run only multiqc and PCGR:
umccrise /bcbio/final/ -T multiqc -T pcgr
```

To exclude stages, use `-E`:

```shell
# Runs all default stages excluding `conpair` report for contamination and T/N concordance
umccrise /bcbio/final/ -E conpair
```

### Custom input

umccrise supports bcbio-nextgen and DRAGEN projects as input. However, you can
also feed custom files as multiple positional arguments. VCF and BAM files are
supported. The sample name will be extracted from VCF and BAM headers. For now,
the VCF file is assumed to contain T/N somatic small variant calls, and the BAM
file is assumed to be from the tumor.

```shell
umccrise sample1.bam sample2.bam sample1.vcf.gz sample3.vcf.gz -o umccrised -j10
```

You can also provide a TSV file as input. If any input file has an extention
`.tsv` (e.g. `umccrise input.tsv`) the file is assumed as a TSV file with a
header, and any of the following columns in arbitrary order:

- `sample`
- `wgs` (WGS tumor BAM, required)
- `normal` (WGS normal BAM, required)
- `exome` (optional tumor BAM)
- `exome_normal` (optional normal BAM)
- `rna` (optional WTS BAM, however required for neoantigens)
- `rna_bcbio` (optional path to RNAseq bcbio workflow, required for neoantigens)
- `rna_sample` (sample name in the RNAseq bcbio workflow above, required for
  neoantigens)
- `somatic_vcf` (tumor/normal somatic VCF calls, optional. If not provided, SAGE
  will be run)
- `germline_vcf` (germline variant calls, optional)
- `sv_vcf` (SV calls, optional)

### Running on selected samples

By default, umccrise will process all batches in the run in parallel. You can
submit only certain samples/batches using `-s`/`--sample` arguments, e.g.:

```shell
# of all samples in a project, takes only sample1 and sample3, plus all corresponding normal/tumor matches:
umccrise /input/project/final -s sample1 -s sample3
```

Or you might want to exclude certain samples/batches with `-e`/`--exclude`:

```shell
# takes all samples in a project, excluding sample1 and sample2 and corresponding normal/tumor matches:
umccrise /input/project -e sample1 -e sample2
```

## Updating

```shell
source load_umccrise.sh            # load the umccrise environment
cd umccrise ; git pull ; cd ..     # if the umccrise codebase changed
```

If dependencies changed:

```shell
conda activate miniconda/envs/umccrise
conda env update -f umccrise/envs/umccrise.yml -p miniconda/envs/umccrise
conda env update -f umccrise/envs/pcgr_linux.yml -p miniconda/envs/umccrise_pcgr
# conda env update -f umccrise/envs/pcgr_macos.yml -p miniconda/envs/umccrise_pcgr  # for macos
conda env update -f umccrise/envs/hmf.yml -p miniconda/envs/umccrise_hmf
```

## Development

Changes pulled in `umccrise` repository clone folder will affect immidiately due
to use of the `-e` option in `pip install -e`. To do the same for other related
packages, you can clone them as well (or move already cloned repos from
`./umccrise/envs/src`, and run `pip install -e` on them as well:

```
source load_umccrise.sh
git clone https://github.com/vladsaveliev/NGS_Utils ngs_utils   ; pip install -e ngs_utils
git clone https://github.com/umccr/reference_data               ; pip install -e reference_data
git clone https://github.com/umccr/vcf_stuff                    ; pip install -e vcf_stuff
```

## Docker

You can pull the ready-to-run docker image from DockerHub:

```shell
docker pull umccr/umccrise:latest
```

An example command to run umccrise on docker could be (although
<abbr title="Your Mileage May Vary">YMMV</abbr>):

```shell
docker run -t --cpus 4 \
    -v=$PWD/umccrise_test_data/results/bcbio_test_project_docker:/output_dir \
    -v=$PWD/umccrise_test_data/data/bcbio_test_project:/bcbio_project \
    -v=/codebuild/output/refdata/genomes:/work/genomes \
    umccr/umccrise /bcbio_project -o /output_dir --genomes /work/genomes
```

This example assumes that:

1. You are running this umccrise container against the
   [umccrise_test_data](https://github.com/umccr/umccrise_test_data)
2. You have figured out the genome data files and directory hierarchy for
   `/work/genomes`. See the
   [building reference data section](#building-reference-data) below.

## Building reference data

To build the bundle from scratch, follow instructions for each kind of data
below.

#### PURPLE

- Download hg19 and hg38 versions of the likely heterozygous sites for AMBER
  from the HMF website at <https://resources.hartwigmedicalfoundation.nl/> `->`
  HMFTools-Resources `->` Amber3 (link valid as of May 2020)

```shell
mv GermlineHetPon.hg19.vcf.gz genomes/GRCh37/hmf
mv GermlineHetPon.hg38.vcf.gz genomes/hg38/hmf
```

- Download hg19 and hg38 versions of GC profile for COBALT from the HMF website
  at <https://resources.hartwigmedicalfoundation.nl/> `->` HMFTools-Resources
  `->` Cobalt (link valid as of May 2020)

```shell
mv GC_profile.hg19.1000bp.cnp.gz genomes/GRCh37/hmf
mv GC_profile.hg38.1000bp.cnp.gz genomes/hg38/hmf
```

#### GNOMAD

Version 2.1 (latest, 500G, hosted by Broad Institute):

```shell
wget -c https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz | \
    # can optionally remove LCR becaues we annoatate vs LCR futher anyway, but the file will be already small enough:
    #bcftools filter -i 'FILTER="PASS" & segdup=0 & lcr=0 & decoy=0 & (AN_popmax>=500 & AF_popmax>=0.01 | \
    #AN_popmax>=100 & AF_popmax>=0.01)' gnomad.genomes.r2.1.sites.vcf.bgz -Ob | \
    bcftools annotate -x ID,^INFO/AN_popmax,^INFO/AF_popmax,FORMAT -Oz -o gnomad_genome.r2.1.common_pass_clean.vcf.gz
tabix -p vcf gnomad_genome.r2.1.common_pass_clean.vcf.gz
```

Normalise (see https://github.com/chapmanb/cloudbiolinux/pull/279, however after
all just 5 indels will be changed, so not a big deal):

```shell
ref=GRCh37.fa
norm_vcf gnomad_genome.r2.1.common_pass_clean.vcf.gz -o gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz --ref-fasta $ref
```

Counts:

```
$ bcftools view gnomad_genome.r2.1.common_pass_clean.vcf.gz -H | wc    # with lcr&segdup&decoy
24671774

$ bcftools view gnomad_genome.common_pass_clean.vcf.gz -H | wc
21273673
```

Convert to hg38

```shell
# spartan
CrossMap.py vcf /data/cephfs/punim0010/extras/hg19ToHg38.over.chain.gz ../GRCh37/gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz hg38.fa gnomad_genome.r2.1.common_pass_clean.norm.vcf.unsorted
bcftools view gnomad_genome.r2.1.common_pass_clean.norm.vcf.unsorted | bcftools sort -Oz -o gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz
tabix -p vcf gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz
```

#### PCGR

The PCGR data bundle gets refreshed every release, so please select the
appropriate one from
[PCGR's README](https://github.com/sigven/pcgr)!

```shell
# Download the data bundles
pip install gdown  # or use  `$(pwd)/miniconda/envs/${ENV_NAME}_pcgr/bin/gdown`
gdown https://drive.google.com/uc?id=<GDOCS_ID_SEE_PCGR_DATABUNDLE_README> -O - | tar xvfz - # hg19
gdown https://drive.google.com/uc?id=<GDOCS_ID_SEE_PCGR_DATABUNDLE_README> -O - | tar xvfz - # hg38

# (Optional) if you are running on AWS, upload the PCGR data bundles to S3 like this:
gdown https://drive.google.com/uc?id=<GDOCS_ID_SEE_PCGR_DATABUNDLE_README> -O - | aws s3 cp - s3://umccr-umccrise-refdata-dev/Hsapiens/GRCh37/PCGR/pcgr.databundle.grch37.YYYMMDD.tgz
gdown https://drive.google.com/uc?id=<GDOCS_ID_SEE_PCGR_DATABUNDLE_README> -O - | aws s3 cp - s3://umccr-umccrise-refdata-dev/Hsapiens/hg38/PCGR/pcgr.databundle.grch38.YYYMMDD.tgz
```

#### Problem regions

Copy GRCh37 from bcbio-nextgen:

```shell
cp -r /g/data/gx8/local/development/bcbio/genomes/Hsapiens/GRCh37/coverage/problem_regions problem_regions
```

Generate SegDup:

```shell
cd problem_regions
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz .
gunzip -c genomicSuperDups.txt.gz | cut -f2,3,4 >> segdup.bed_tmp
gunzip -c genomicSuperDups.txt.gz | cut -f8,9,10 >> segdup.bed_tmp
grep -v gl segdup.bed_tmp | sed 's/chr//' | bedtools sort -i - | bedtools merge -i - > segdup.bed
bgzip -f segdup.bed && tabix -f -p bed segdup.bed.gz
rm segdup.bed_tmp genomicSuperDups.txt.gz
```

Generate ENCODE:

```shell
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
gunzip -c ENCFF356LFX.bed.gz | bgzip -c > hg38/problem_regions/ENCODE/encode4_unified_blacklist.bed.gz
tabix -p bed hg38/problem_regions/ENCODE/encode4_unified_blacklist.bed.gz
rm ENCFF356LFX.bed.gz
```

Generate the regions for variant calling: hg38-noalt chromosomes excluding the
ENCODE blacklist

```shell
bedtools subtract -a hg38_noalt.bed -b problem_regions/ENCODE/encode4_unified_blacklist.bed.gz > hg38_noalt_noBlacklist.bed
```

The blacklist removes ~2.3% of the noalt genome size:

```
$ bedsize hg38_noalt_noBlacklist.bed
3016716091

$ bedsize hg38_noalt.bed
3088286376
```

The blacklisted regions contain no hotspots:

```shell
bedtools intersect \
    -a problem_regions/ENCODE/encode4_unified_blacklist.bed.gz\
    -b hotspots/merged.vcf.gz
```

Lift over to hg38:

```shell
convert () {
    f=$(basename $1)
    echo "Processing $f"
    zless $1 | py -x "('chr' + x) if not x.startswith('MT') else 'chrM'" | grep -v chrG > $f.hg19
    CrossMap.py bed /g/data3/gx8/extras/hg19ToHg38.over.chain.gz $f.hg19 $f.unsorted
    bedtools sort -i $f.unsorted | bgzip -c > $f
    tabix -p vcf $f
}

convert ../../GRCh37/problem_regions/segdup.bed.gz

mkdir GA4GH repeats

cd GA4GH
for fp in $(ls ../../../GRCh37/problem_regions/GA4GH/*.bed.gz) ; do convert $fp ; done

cd ../repeats
convert ../../../GRCh37/problem_regions/repeats/LCR.bed.gz
convert ../../../GRCh37/problem_regions/repeats/polyx.bed.gz

cat ../../../GRCh37/problem_regions/repeats/sv_repeat_telomere_centromere.bed | py -x "('chr' + x) if not x.startswith('MT') else 'chrM'" | grep -v chrG > sv_repeat_telomere_centromere.bed_hg19
CrossMap.py bed /g/data3/gx8/extras/hg19ToHg38.over.chain.gz sv_repeat_telomere_centromere.bed_hg19 sv_repeat_telomere_centromere.bed_unsorted
bedtools sort -i sv_repeat_telomere_centromere.bed_unsorted | > sv_repeat_telomere_centromere.bed
```

#### Coding regions (SAGE)

```shell
cd GRCh37/hmf
generate_bed.py -g GRCh37 \
   --principal --key-genes --features CDS | sort -k1,1V -k2,2n | grep -v ^MT | grep -v ^GL \
   | bedtools merge -c 4 -o collapse -i - > coding_regions.bed

cd hg38/hmf
generate_bed.py -g hg38 \
   --principal --key-genes --features CDS | sort -k1,1V -k2,2n | grep -v ^chrM \
   | bedtools merge -c 4 -o collapse -i - \
   > coding_regions.bed
```

#### Ensembl annotation

Using pyensembl package:

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

- Hartwig's: 10211 changes in 3650 locations,
- PCGR: 10627 changes in 2494 locations.
- Overlap: 2960 changes in 968 locations.

The overlap is small, so we better merge sources into a single VCF.

First, download HMF TSV file and convert to VCF:

```shell
wget https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMF-Pipeline-Resources&files=KnownHotspots.tsv.gz -O KnownHotspots.tsv.gz

echo "##fileformat=VCFv4.2" > hmf.vcf
echo '##INFO=<ID=HMF,Number=.,Type=Flag,Description="Hotspot is from HMF">' >> hmf.vcf
echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> hmf.vcf
gunzip -c KnownHotspots.tsv.gz | py -x "print('\t'.join([x.split()[0], x.split()[1], '.', x.split()[2], x.split()[3], '.', '.', 'HMF']))" >> hmf.vcf
bgzip hmf.vcf
tabix -p vcf hmf.vcf.gz
```

Prepare PCGR hotspots:

```shell
SRC=/Users/vsaveliev/bio/genomes/pcgr/data/grch37/cancer_hotspots/cancer_hotspots.vcf.gz
bcftools view -h $SRC | grep ^## > pcgr.vcf
echo '##INFO=<ID=PCGR,Number=.,Type=Flag,Description="Hotspot is from PCGR (cancerhotspots.org_v2)">' >> pcgr.vcf
bcftools view -h $SRC | grep ^#CRHOM >> pcgr.vcf
bcftools view -H $SRC | bioawk -t '{ print $1,$2,$3,$4,$5,$6,$7,$8";PCGR" }' >> pcgr.vcf

bgzip pcgr.vcf
tabix -p vcf pcgr.vcf.gz
```

Merge:

```shell
bcftools merge -m none hmf.vcf.gz pcgr.vcf.gz -Oz -o merged.vcf.gz
tabix -p vcf merged.vcf.gz
# Adding into the workflows repo:
gunzip -c merged.vcf.gz | grep -v ^## > /Users/vsaveliev/git/umccr/workflows/genes/hotspots/hotspots.tsv
```

Convert to hg38

```shell
cd ../hg38/hotspots
INP=../../GRCh37/hotspots/merged.vcf.gz

gunzip -c $INP \
| py -x "x.replace('##contig=<ID=', '##contig=<ID=chr') if x.startswith('#') else 'chr' + x" \
| py -x "x.replace('chrMT', 'chrM')" \
| grep -v chrG \
| gzip -c > merged_hg19.vcf.gz

CrossMap.py vcf /g/data3/gx8/extras/hg19ToHg38.over.chain.gz merged_hg19.vcf.gz ../hg38.fa merged_unsorted.vcf
bcftools sort merged_unsorted.vcf -Oz -o merged.vcf.gz
tabix -p vcf merged.vcf.gz
```

#### Other HMF files

Download
`NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz`
and `out_150_hg19.mappability.bed.gz` from
<https://resources.hartwigmedicalfoundation.nl/> `->` HMFTools-Resources `->`
Sage (link valid as of May 2020)

To hg38:

```shell
convert ../../GRCh37/hmf/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz
convert ../../GRCh37/hmf/out_150_hg19.mappability.bed.gz
```

#### Fusions

We use
[HMF fusions](https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd?path=%2FHMF-Pipeline-Resources)
for SV prioritization. See `NGS_Utils/ngs_utils/refernece_data/__init__.py` for
details.

#### SnpEff

```shell
cd GRCh37
mkdir snpeff
cd snpeff
wget https://sourceforge.net/projects/snpeff/files/databases/v4_3/snpEff_v4_3_GRCh37.75.zip
unzip *.zip
rm *.zip

cd hg38
mkdir snpeff
cd snpeff
wget https://sourceforge.net/projects/snpeff/files/databases/v4_3/snpEff_v4_3_GRCh38.92.zip
unzip *.zip
rm *.zip
```

### DVC

We use [DVC](https://dvc.org/doc) to track reference data in the
[reference_data](https://github.com/umccr/reference_data) repo

```
git clone git@github.com:umccr/reference_data.git reference_data.git
cd reference_data.git

# The first time we did:
rsync -trv /g/data/gx8/extras/umccrise_genomes/hg38 genomes/
dvc init
dvc add genomes/hg38
# Alternatively, could used `run`:
#dvc run -n copy_hg38 -o genomes/hg38 rsync -trv /g/data/gx8/extras/umccrise_genomes/hg38 genomes/
dvc remote add -d storage s3://umccr-refdata-dev/dvc-storage
dvc push

# To pull the ref data, do:
dvc pull

# To update, do:
rsync -trv /g/data/gx8/extras/umccrise_genomes/hg38 genomes/
dvc add genomes/hg38
git add genomes/.gitignore genomes/hg38.dvc
dvc push
```
