Umccrise
--------

[![Build Status](https://travis-ci.org/umccr/umccrise.svg?branch=master)](https://travis-ci.org/umccr/umccrise)

Umccrise is designed to post-processess outputs from [bcbio-nextgen](https://github.com/chapmanb/bcbio-nextgen) cancer variant calling analysis pipeline.

- Filter artefacts and germline leakage from somatic calls
- Run [PCGR](https://github.com/sigven/pcgr) to annotate, prioritize and report somatic variants
- Run [CPSR](https://github.com/sigven/cpsr) to annotate, prioritize and report germline variants
- Filter, annotate, prioritize and report SV
- Run [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) to call CNV, purity, ploidy, and recover SV
- Generate a MutliQC report comparing QC to "reference" samples
- Generate a cancer report with mutational signatures, strand bias analysis, PURPLE results, and prioritized SVs
- Run [CACAO](https://github.com/sigven/cacao) to calculate coverage in common hotspots, as well as goleft to estimate coverage problems
- Run [Conpair](https://github.com/nygenome/Conpair) to tumor/normal concordance and sample contamination

See the [workflow.md](workflow.md) for the complete description of the workflow.

See the [HISTORY.md](HISTORY.md) for the version history.


Contents:

- [Umccrise](#umccrise)
- [Installation](#installation)
- [Reference data](#reference-data)
- [Testing](#testing)
- [UMCCRISE on AWS](#umccrise-on-aws)
- [UMCCR HPC](#umccr-hpc)
- [Usage](#usage)
    - [Run selected steps](#run-selected-steps)
    - [Run on selected samples](#run-on-selected-samples)
    - [Use HPC cluster](#use-hpc-cluster)
- [Updating](#updating)
- [Development](#development)
- [Docker](#docker)
- [Building reference data](#building-reference-data)
    - [GNOMAD](#gnomad)
    - [PCGR](#pcgr)
    - [Problem regions](#problem-regions)
    - [Coding regions (SAGE)](#coding-regions-sage)
    - [Ensembl annotation](#ensembl-annotation)
    - [Hotspots](#hotspots)
    - [Other HMF files](#other-hmf-files)
    - [Fusions](#fusions)


## Installation

Run the following to create a directory "umccrise" and install into it

```
mkdir umccrise
cd umccrise
source <(curl -s https://github.com/umccr/umccrise/blob/master/install.sh)
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

## UMCCRISE on AWS

umccrise on AWS is run via AWS Batch in a defined compute environment. This is set up and maintained via the [umccrise Terraform Stack][umccrise_tf_stack]. This stack also defines the version of umccrise that is used within AWS and how umccrise jobs are triggered.

## UMCCR HPC

Load in Raijin

```
source /g/data3/gx8/extras/umccrise/load_umccrise.sh
```

Load on Spartan

```
source /data/cephfs/punim0010/extras/umccrise/load_umccrise.sh
```

## Usage

Runs the patient analysis pipeline on bcbio-nextgen `final` folder.

```
umccrise /path/to/bcbio/project/final -j 30   # run using 30 CPUs
```

The output will be created in `umccrised` folder. To override, use `-o`:

```
umccrise /path/to/bcbio/project/final -o umccrised_results
```

#### Run selected steps

Umccrise workflow consists of the following steps: `pcgr`, `coverage`, `structural`, `small_variants`, `rmd`, `multiqc`, `purple`, `igv`.

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
# conda env update -f umccrise/envs/pcgr_macos.yml -p miniconda/envs/umccrise_pcgr  # for macos
conda env update -f umccrise/envs/purple.yml -p miniconda/envs/umccrise_purple
```

## Development

Changes pulled in `umccrise` repository clone folder will affect immidiately due to use of the `-e` option in `pip install -e`. To do the same for other related packages, you can clone them as well (or move already cloned repos from `./umccrise/envs/src`, and run `pip install -e` on them as well:

```
source load_umccrise.sh
git clone https://github.com/vladsaveliev/NGS_Utils ngs_utils   ; pip intall -e ngs_utils
git clone https://github.com/umccr/hpc_utils                    ; pip intall -e hpc_utils
git clone https://github.com/vladsaveliev/vcf_stuff             ; pip intall -e vcf_stuff
```


## Docker

You can pull the ready-to-run docker image from DockerHub:

```
docker pull umccr/umccrise:latest
```

## Building reference data

To build the bundle from scratch, follow instructions for each kind of data below.


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



[umccrise_tf_stack]: https://github.com/umccr/infrastructure/tree/master/terraform/stacks/umccrise

