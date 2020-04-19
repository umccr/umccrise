Umccrise
--------

UMCCR cancer reporting

[![Build Status](https://travis-ci.org/umccr/umccrise.svg?branch=master)](https://travis-ci.org/umccr/umccrise)

Umccrise is developed to post-processess outputs from cancer variant calling analysis pipelines from 2 platforms: ([bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen) and Illumina Dragen, and generate reports helpful for researchers and curators at UMCCR.

- Filter artefacts and germline leakage from somatic variant calls;
- Run [PCGR](https://github.com/sigven/pcgr) to annotate, prioritize and report somatic variants;
- Run [CPSR](https://github.com/sigven/cpsr) to annotate, prioritize and report germline variants;
- Filter, annotate, prioritize and report structural variants (SVs);
- Run [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) to call copy number calls (CNV), purity, ploidy, and recover SVs;
- Generate a MutliQC report that summarizes quality control statistics in context of background "gold standard" samples;
- Generate a cancer report with mutational signatures, circos plots, prioritized copy number and structural variant calls;
- Run [CACAO](https://github.com/sigven/cacao) to calculate coverage in common hotspots, as well as goleft to estimate coverage problems;
- Run [Conpair](https://github.com/nygenome/Conpair) to tumor/normal concordance and sample contamination.

See the [workflow.md](workflow.md) for a detailed description of the workflow.

See the [HISTORY.md](HISTORY.md) for the version history.

Contents:
- [Umccrise](#umccrise)
- [Installation](#installation)
- [Reference data](#reference-data)
- [Testing](#testing)
- [Usage](#usage)
- [AWS](#aws)
- [HPC (NCI Gadi)](#hpc--nci-gadi-)
    + [Run selected steps](#run-selected-steps)
    + [Run on selected samples](#run-on-selected-samples)
- [Updating](#updating)
- [Development](#development)
- [Docker](#docker)
- [Building reference data](#building-reference-data)
    + [PURPLE](#purple)
    + [GNOMAD](#gnomad)
    + [PCGR](#pcgr)
    + [Problem regions](#problem-regions)
    + [Coding regions (SAGE)](#coding-regions--sage-)
    + [Ensembl annotation](#ensembl-annotation)
    + [Hotspots](#hotspots)
    + [Other HMF files](#other-hmf-files)
    + [Fusions](#fusions)
- [GRIDSS and LINX](#gridss-and-linx)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

## Installation

Run the following to create a directory "umccrise" and install into it

```
mkdir umccrise
cd umccrise
git clone https://github.com/umccr/umccrise
source umccrise/install.sh
```

It will generate `load_umccrise.sh` script that can be sourced to load the umccrise environment:

```
source load_umccrise.sh
```

## Reference data

Also you will need the reference data. Umccrise automatically finds reference data on Spartan and NCI environments, as well as the reference data bundle mounted to the Docker image under `/genomes`. You can specify a custom path with `--genomes <path>` (the path can be a tarball, which is useful for runs in a pipeline). You can copy the genome bundle from NCI (`/g/data3/gx8/extras/umccrise_2019_Mar/genomes`) anywhere else. To build the bundle from scratch, see [below](#building-reference-data).

For the reference, to sync the reference data from Spartan to NCI, you can use:

```
cd /data/cephfs/punim0010/extras/umccrise
rsync -rv --size-only genomes/ rjn:/g/data3/gx8/extras/umccrise/genomes
```

## Testing

Load the umccrise environment, clone the repo with toy test data, and run nosetests: 

```
source load_umccrise.sh
git clone https://github.com/umccr/umccrise_test_data
TEST_OPTS="-c -j2" nosetests -s umccrise_test_data/test.py
```

## Usage

```
umccrise <input-folder> -o umccrised -j10
```

Where <input-folder> can be either bcbio-nextgen analysis results' `final` folder, or Dragen analysis results folder.

`-j10` specifies the number of cores the pipeline will attempt to use.


## AWS

Umccrise on AWS is run via AWS Batch in a defined compute environment. This is set up and maintained via the [umccrise Terraform Stack][https://github.com/umccr/infrastructure/tree/master/terraform/stacks/umccrise]. This stack also defines the version of umccrise that is used within AWS and how umccrise jobs are triggered.

## HPC (NCI Gadi)

The production version can be loaded by sourcing the following:

```
source /g/data3/gx8/extras/umccrise/load_umccrise.sh
```

There are also installations for each somewhat major release available in `/g/data3/gx8/extras`, for instance by sourcing the following you can load version 0.17, released on Jan 2020:

```
/g/data/gx8/extras/umccrise_017_2020_Jan/load_umccrise.sh
```

To parallelize the tool using the cluster scheduler, you can use `--cluster-auto` (`-c`) option:

```
umccrise <input-folder> -j30 -c
```

Alternatively, you can specify a custom submission template with `--cluster-cmd`, e.g.:

```
umccrise <input-folder> -j30 --cluster-cmd "sbatch -p vccc -n {threads} -t 24:00:00 --mem {resources.mem_mb} -J umccrise"
```

Make sure to use `-j` outside of that template: this options tells snakemake how many cores is allowed to use by the entire pipeline at a single moment.


#### Run selected steps

Umccrise workflow consists of the following steps: `pcgr`, `coverage`, `structural`, `small_variants`, `rmd`, `multiqc`, `purple`, `igv`.

To run just a particular step (or steps), use:

```
umccrise <input-folder> <step_name>
```

Where `<step_name>` is from the list above. E.g.:

```
umccrise <input-folder> pcgr
```

Note that the `igv` step (preparing minibams and uploading them to `s3://umccr-igv`) takes ~5 hours for a WGS sample compared to ~20 minutes for all other steps combined. For that reason, it is always executed in the end of the pipeline, so you can expect that when it is being executed, all other output is ready.

#### Run on selected samples

By default, Umccrise will process all batches in the run in parallel. You can submit only certain samples/batchs using `--sample` or `--batch` arguments, e.g.:

```
umccrise <input-folder> --batch cup-batch
umccrise <input-folder> --sample cup-tumor_1,cup-tumor_2
```

Or you might want to exclude certain samples/batches with `--exclude`:

```
umccrise <input-folder> --exclude cup-tumor_1,cup-batch_2
```

## Updating

```
source load_umccrise.sh            # load the umccrise environment
cd umccrise ; git pull ; cd ..     # if the code base changed
```

If dependencies changed:

```
conda activate miniconda/envs/umccrise
conda env update -f umccrise/envs/umccrise.yml -p miniconda/envs/umccrise
conda env update -f umccrise/envs/pcgr_linux.yml -p miniconda/envs/umccrise_pcgr
conda env update -f umccrise/envs/python2.yml -p miniconda/envs/umccrise_python2
# conda env update -f umccrise/envs/pcgr_macos.yml -p miniconda/envs/umccrise_pcgr  # for macos
conda env update -f umccrise/envs/hmf.yml -p miniconda/envs/umccrise_hmf
```

## Development

Changes pulled in `umccrise` repository clone folder will affect immidiately due to use of the `-e` option in `pip install -e`. To do the same for other related packages, you can clone them as well (or move already cloned repos from `./umccrise/envs/src`, and run `pip install -e` on them as well:

```
source load_umccrise.sh
git clone https://github.com/vladsaveliev/NGS_Utils ngs_utils   ; pip install -e ngs_utils
git clone https://github.com/umccr/hpc_utils                    ; pip install -e hpc_utils
git clone https://github.com/vladsaveliev/vcf_stuff             ; pip install -e vcf_stuff
```


## Docker

You can pull the ready-to-run docker image from DockerHub:

```
docker pull umccr/umccrise:latest
```

An example command to run umccrise on docker could be (although <abbr title="Your Mileage May Vary">YMMV</abbr>):

```shell
$ docker run -t --cpus 4 -v=$PWD/umccrise_test_data/results/bcbio_test_project_docker:/output_dir
	         -v=$PWD/umccrise_test_data/data/bcbio_test_project:/bcbio_project
             -v=/codebuild/output/refdata/genomes:/work/genomes
             umccr/umccrise /bcbio_project -o /output_dir --genomes /work/genomes
```

This example assumes that:

	1. You are running this umccrise container against the [umccrise_test_data](https://github.com/umccr/umccrise_test_data)
	2. You have figured out the genome data files and directory hierarchy for `/work/genomes`. See [building reference data section](#building-reference-data) below.


## Building reference data

To build the bundle from scratch, follow instructions for each kind of data below.


#### PURPLE

Download hg19 and hg38 versions of the likely heterozygous sites for AMBER from the HMF website:

https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMFTools-Resources%2FAmber3&files=GermlineHetPon.hg19.vcf.gz
https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMFTools-Resources%2FAmber3&files=GermlineHetPon.hg38.vcf.gz

mv GermlineHetPon.hg19.vcf.gz genomes/GRCh37/hmf
mv GermlineHetPon.hg38.vcf.gz genomes/hg38/hmf

Download hg19 and hg38 versions of GC profile for COBALT from the HMF website:

https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMFTools-Resources%2FCobalt&files=GC_profile.hg19.1000bp.cnp.gz
https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd/download?path=%2FHMFTools-Resources%2FCobalt&files=GC_profile.hg38.1000bp.cnp.gz

mv GC_profile.hg19.1000bp.cnp.gz genomes/GRCh37/hmf
mv GC_profile.hg38.1000bp.cnp.gz genomes/hg38/hmf


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
CrossMap.py vcf /data/cephfs/punim0010/extras/hg19ToHg38.over.chain.gz ../GRCh37/gnomad_genome.r2.1.common_pass_clean.norm.vcf.gz hg38.fa gnomad_genome.r2.1.common_pass_clean.norm.vcf.unsorted
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

Updating PCGR packages:

```
source /g/data/gx8/extras/umccrise_017_2020_Jan_dev/load_umccrise.sh
export PATH=/g/data/gx8/extras/umccrise_017_2020_Jan_dev/miniconda/envs/umccrise_pcgr/bin:$PATH
# conda install conda-build anaconda

cd /g/data/gx8/extras/umccrise_017_2020_Jan_dev/pcgr/install_no_docker
export VERSION=0.8.4.8
conda build conda_pkg/pcgr
conda build conda_pkg/pcgr_dockerized
conda convert --platform osx-64 \
    /g/data/gx8/extras/umccrise_017_2020_Jan_dev/miniconda/envs/umccrise/conda-bld/linux-64/pcgr_dockerized-$VERSION-*.tar.bz2 \
    --output-dir /g/data/gx8/extras/umccrise_017_2020_Jan_dev/miniconda/envs/umccrise/conda-bld/
anaconda upload --force -u pcgr /g/data/gx8/extras/umccrise_017_2020_Jan_dev/miniconda/envs/umccrise/conda-bld/*/*-$VERSION-*.tar.bz2

cd /g/data/gx8/extras/umccrise_017_2020_Jan_dev/cpsr
export VERSION=0.5.2.4
conda build conda_pkg/cpsr
conda build conda_pkg/cpsr_dockerized
conda convert --platform osx-64 \
    /g/data/gx8/extras/umccrise_017_2020_Jan_dev/miniconda/envs/umccrise/conda-bld/linux-64/cpsr_dockerized-$VERSION-*.tar.bz2 \
    --output-dir /g/data/gx8/extras/umccrise_017_2020_Jan_dev/miniconda/envs/umccrise/conda-bld/
anaconda upload --force -u pcgr /g/data/gx8/extras/umccrise_017_2020_Jan_dev/miniconda/envs/umccrise/conda-bld/*/*-$VERSION-*.tar.bz2
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

Generate ENCODE

```
git clone https://github.com/Boyle-Lab/Blacklist
gunzip -c Blacklist/lists/hg19-blacklist.v2.bed.gz | py -x "x[3:]" > Blacklist/lists/GRCh37-blacklist.v2.bed

bgzip Blacklist/lists/GRCh37-blacklist.v2.bed -c > GRCh37/problem_regions/ENCODE/blacklist.v2.bed.gz
gunzip -c Blacklist/lists/hg38-blacklist.v2.bed.gz | bgzip -c > hg38/problem_regions/ENCODE/blacklist.v2.bed.gz
tabix -p bed GRCh37/problem_regions/ENCODE/blacklist.v2.bed.gz
tabix -p bed hg38/problem_regions/ENCODE/blacklist.v2.bed.gz
rm -rf Blacklist
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
# Adding into the workflows repo:
gunzip -c merged.vcf.gz | grep -v ^## > /Users/vsaveliev/git/umccr/workflows/genes/hotspots/hotspots.tsv
```

Convert to hg38

```
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

Download NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz and out_150_hg19.mappability.bed.gz 
from https://resources.hartwigmedicalfoundation.nl

To hg38:

```
convert ../../GRCh37/hmf/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz
convert ../../GRCh37/hmf/out_150_hg19.mappability.bed.gz
```

#### Fusions

We use [HMF fusions](https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd?path=%2FHMF-Pipeline-Resources) for SV prioritization. See `NGS_Utils/ngs_utils/refernece_data/__init__.py` for details.

#### SnpEff

```
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

## GRIDSS and LINX

[GRIDSS+PURPLE+LINX](https://github.com/hartwigmedical/gridss-purple-linx) is a chain of tools developed by Hartwig Medical Foundation:

- [GRIDSS](https://github.com/PapenfussLab/gridss) a structural variant caller,
- [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) a copy number, ploidy and purity caller,
- [LINX](https://github.com/hartwigmedical/hmftools/tree/master/sv-linx) a tools for structural event classification and visualisation.

We have a Docker-free wrapper for this chain of tools, which can be run on NCI HPC cluster. To use the wrapper, first load Umccrise with:

```
source /g/data/gx8/extras/umccrise_017_2020_Jan/load_umccrise.sh
 ```

Then run the `gpl` command:

```
gpl -N <normal bam> \
    -T <tumor bam> \
    -S <small somatic variants vcf> \
    -s <tumor sample name in the vcf above> \
    -o <output dir> \
    -t<threads> 
```

For instance, to process this sample with known HPV viral integration:

```
gpl -N /g/data3/gx8/projects/Saveliev_Viral/All_WGS_SBJ00174/SBJ00174_MDX190157_L1900839-ready.bam \
    -T /g/data3/gx8/projects/Saveliev_Viral/All_WGS_SBJ00174/SBJ00174_MDX190158_L1900840-ready.bam \
    -S /g/data3/gx8/projects/Saveliev_Viral/All_WGS_SBJ00174/SBJ00174__SBJ00174_MDX190158_L1900840-somatic-ensemble-PASS.vcf.gz \
    -s SFRC01189 \
    -o /g/data3/gx8/projects/Saveliev_Viral/All_WGS_SBJ00174/GRIDSS \
    -t10 
```

The wrapper loads the [conda environment](blob/master/envs/hmf.yml) that has all the dependencies for the HMF tools (the environment is installed on Gadi at `/g/data/gx8/extras/umccrise_017_2020_Jan_dev/miniconda/envs/umccrise_hmf`), runs the `gridss-purple-linx.sh` script from [our fork of gridss-purple-linx](https://github.com/vladsaveliev/gridss-purple-linx), pointing it to the data bundle `/g/data3/gx8/extras/umccrise_017_2020_Jan_dev/genomes/GRCh37/hmf/gridss`.

Note that we are using forks of GRIDSS and gridss-purple-linx: 

- GRIDSS: https://github.com/vladsaveliev/gridss
- gridss-purple-linx: https://github.com/vladsaveliev/gridss-purple-linx

Where we have some amendments of `gridss.sh` and `gridss-purple-linx.sh` scripts to make it work with an HPC installation. Eventually we want to syncronise them with the original repos.

There are a few drawbacks that prevents us from integrating the pipeline in production. The pipeline is currently unstable (GRIDSS might occasianlly crash on certain samples, plus it is quite slow and resource demanding, so your cluster node might expire before it completes). Also it [supports only GRCh37](https://github.com/hartwigmedical/gridss-purple-linx/issues/3) - even though GRIDSS itself and PURPLE are genome-agnostic, annotations needed for LINX are put up only for GRCh37 at the moment.

The idea is to eventually complement or even subsitute the [structural variants](https://github.com/umccr/umccrise/blob/master/workflow.md#structural-variants) workflow in Umccrise, including prioritization and the section in the cancer report. In addition, we hope to use it as the viral integration tool instead of current [semi-manual workflow](https://github.com/umccr/oncoviruses). The downsides of that are:

- No hg38 support for LINX
- GRIDSS is very slow
- GRIDSS is not very stable
- HMF tools are not adapted to FFPE data
- LINX doesn't annotate events that don't affect particular gene sequence, but rather downstream or upstream
 




