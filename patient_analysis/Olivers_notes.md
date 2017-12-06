# Patient analysis notes

Set up result/sample names:

```
BCINSTALL="/data/projects/punim0010/local/share/bcbio/"
BCRESULT="/data/cephfs/punim0010/data/Results/Avner/WPT-013/final"
BCFINAL="2017-10-19_WPT-013"
BCPOST="/data/cephfs/punim0010/projects/Hofmann_Explore/testrun"
EXTRAS="/data/projects/punim0010/local/share/extras/"

# This information should come from the *.csv
BCTUMOR="WPT-013-organoid"
BCNORMAL="WPT-013-normal"
BCBATCH="batch1"
```


Create a subdirectory in the `final` section of a given run.

```
cd ${BCPOST}
mkdir analysis
cd analysis
```

## QC

Start with the MultiQC results and pull them into the current folder.

```
mkdir qc
cd qc
ln ${BCRESULT}/${BCFINAL}/multiqc/multiqc_report.html .
cd ..
```

## Cancer gene coverage

Looking at coverage for a limited set of (cancer) genes to assess overall reliability. Minimum coverage for normal is 10, 30 for cancer.

```
mkdir coverage
cd coverage
 
goleft depth --reference ${BCINSTALL}/genomes/Hsapiens/GRCh37/seq/GRCh37.fa --processes 32 --bed ${BCINSTALL}/genomes/Hsapiens/GRCh37/coverage/prioritize/cancer/az300.bed.gz --stats --mincov 10 --prefix ${BCNORMAL} ${BCRESULT}/${BCNORMAL}/${BCNORMAL}-ready.bam 

goleft depth --reference ${BCINSTALL}/genomes/Hsapiens/GRCh37/seq/GRCh37.fa --processes 32 --bed ${BCINSTALL}/genomes/Hsapiens/GRCh37/coverage/prioritize/cancer/az300.bed.gz --stats --mincov 30 --prefix ${BCTUMOR} ${BCRESULT}/${BCTUMOR}/${BCTUMOR}-ready.bam 
```

Also bringing in global coverage plots for review (tumor only, quick check for CNVs):

```
find ${BCRESULT}/ -name *ready.bam* -exec ln {} . \;
goleft indexcov --directory indexcov *.bam
rm *.bam*
cd ..
```

## Somatic calls

Run ensemble and raw calls through PCGR and combine with the cnv calls:

```
mkdir somatic
cd somatic

ln ${BCRESULT}/${BCTUMOR}/*.cns .
ln ${BCRESULT}/${BCFINAL}/${BCBATCH}*vcf* .
```

PCGR struggles with anything but the basic chromosome setup. It also ignores any variant not marked `PASS` so might as well remove others to save on transfer times:

```
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y ${BCBATCH}-ensemble-annotated.vcf.gz > ${BCBATCH}-ensemble-annotated-subset.vcf
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y -f .,PASS ${BCBATCH}-mutect2-annotated.vcf.gz > ${BCBATCH}-mutect2-annotated-subset.vcf
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y -f .,PASS ${BCBATCH}-vardict-annotated.vcf.gz > ${BCBATCH}-vardict-annotated-subset.vcf
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y -f .,PASS ${BCBATCH}-varscan-annotated.vcf.gz > ${BCBATCH}-varscan-annotated-subset.vcf

mkdir PCGR_Input
mv *subset* PCGR_Input
```

PCGR also wants a slightly different format for the CNS data:

```
echo -e "Chromosome\tStart\tEnd\tSegment_Mean" > PCGR_Input/${BCBATCH}-cnvkit-pcgr.tsv
cat ${BCBATCH}-cnvkit.cns | grep -v ^chromosome | cut -f 1,2,3,5 >> PCGR_Input/${BCBATCH}-cnvkit-pcgr.tsv
```

Cleanup: get rid of the original files:

```
rm *.cns
```

Ensemble calls to be submitted to CGI. Sharing functionality does not seem to work at this point. Ensemble calls only include variants that `PASS` so no additional filtering required.

Finally, for the local analysis with MutationalPatterns generate UCSC-versions of the ensemble calls:

```
mkdir rstudio

cat PCGR_Input/${BCBATCH}-ensemble-annotated-subset.vcf | awk '{ if($0 !~ /^#/) print "chr"$0; else if(match($0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print $0 }'  | grep -v chrG > rstudio/${BCBATCH}-ensemble-annotated-subset-ucsf.vcf

cd ..
```

All done; file needs to be submitted to CGI and PCGR manually for now. 


## Allelic frequencies

AF is not yet integrated into PCGR or CGI. We can extract those from VarDict or FreeBayes for plotting purposes, but still need to look up AF for genes of interest manually. Not entirely ideal but for a rough summary plot this is going to be sufficient. Can revisit once we've unified AF information across callers:

```
# NOTE: The sample name is hardcoded here. For the life of me I
# cannot get the environmental $BCTUMOR variable injected into
# the vawk command
mkdir af
cd af

## VarDict
zcat ${BCRESULT}/${BCFINAL}/${BCBATCH}-vardict-annotated.vcf.gz | bcftools view -f .,PASS | vawk '{ print S$WPT-013-organoid$AF }' > af_vardict_tumor.txt

bcftools view -f .,PASS ${BCRESULT}/${BCFINAL}/${BCBATCH}-vardict-annotated.vcf.gz | bedtools intersect -a stdin -b ${BCINSTALL}/genomes/Hsapiens/GRCh37/coverage/prioritize/cancer/az300.bed.gz -header | vawk '{ print $1,$2,$3,$4,$5,S$WPT-013-organoid$AF, INFO$ANN}' > af_vardict_az300_batch1.txt 

cd ..
```

Generating the plot manually at a later stage.


## Germline summary

Take a single germline call set (or preferably the ensemble best practice calls, if generated) and annotate any events found in Sean's 105/106 cancer predisposition gene set.

Note: need to store that cancer gene list somewhere. Happy to add to bcbio directly but then need source information.

```
mkdir germline
cd germline

cp ${EXTRAS}/cancer_genes_HUGO.txt .
cp ${EXTRAS}/cancer_genes_ENSG.txt .

zgrep ^# ${BCRESULT}/${BCFINAL}/${BCNORMAL}-ensemble-annotated.vcf.gz > ${BCNORMAL}-ensemble-annotated-cancer-only.vcf && zgrep -f cancer_genes_ENSG.txt ${BCRESULT}/${BCFINAL}/${BCNORMAL}-ensemble-annotated.vcf.gz >> ${BCNORMAL}-ensemble-annotated-cancer-only.vcf

# Pack up the subset 
bgzip ${BCNORMAL}-ensemble-annotated-cancer-only.vcf 
tabix ${BCNORMAL}-ensemble-annotated-cancer-only.vcf.gz

# Prepare it for submission to PCGR
mkdir PCGR_Input
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y -f .,PASS ${BCNORMAL}-ensemble-annotated-cancer-only.vcf.gz > PCGR_Input/${BCNORMAL}-ensemble-annotated-cancer-only-subset.vcf

cd ..
```

Can review the `.txt` file in Excel or submit the VCF to frameworks such as [gene.iobio](http://gene.iobio.io/). 


## Structural changes

Re-do the CNV plots. This will need lots of love (drop gene names, make the scatterplot viable again, etc.).

```
mkdir structural
cd structural
cp ${BCRESULT}/${BCTUMOR}/${BCBATCH}-cnvkit.cns .

# Remove gene labels
head -n 1 ${BCBATCH}-cnvkit.cns > ${BCBATCH}-cnvkit_noLabels.cns
grep -v '^chromosome' ${BCBATCH}-cnvkit.cns  | awk -v OFS='\t' '{ print $1, $2, $3, "", $5, $6, $7, $8 }' >> ${BCBATCH}-cnvkit_noLabels.cns

cnvkit.py diagram -s ${BCBATCH}-cnvkit_noLabels.cns
```

Bring in the prioritized SV calls from Manta. This will change with the inclusion of BPI and should also include a basic plot at some stage:

```
cp ${BCRESULT}/${BCTUMOR}/*prior* .
cat ${BCBATCH}-sv-prioritize.tsv | grep manta > ${BCBATCH}-sv-prioritize-manta.tsv
```

At least for the (most conservative) manta calls generate a file for viewing in Ribbon:

```
awk -v OFS='\t' {'print $1,$2'} ${BCINSTALL}/genomes/Hsapiens/hg19/seq/hg19.fa.fai > hg19.genome

bcftools view -f .,PASS ${BCBATCH}-sv-prioritize-manta.vcf.gz | ${EXTRAS}/svtools-master/vcfToBedpe | cut -f 1-3 | bedtools slop -b 5000 -i stdin -g hg19.genome > manta_regions.bed

bcftools view -f .,PASS ${BCBATCH}-sv-prioritize-manta.vcf.gz | ${EXTRAS}/svtools-master/vcfToBedpe | cut -f 4-6 | grep -v 'CHROM' | bedtools slop -i stdin -g hg19.genome -b 5000 >> manta_regions.bed

bedtools sort -i manta_regions.bed > manta_regions_sorted.bed
bedtools merge -i manta_regions_sorted.bed > ${BCBATCH}_manta_regions_merged.bed

# And one in BEDPE format
bcftools view -f .,PASS ${BCBATCH}-sv-prioritize-manta.vcf.gz | ${EXTRAS}/svtools-master/vcfToBedpe | cut -f 1-7 > ${BCBATCH}-sv-prioritize-manta.bedpe

cd ..
```


## IGV

Create BAM and VCF files suitable for moving around easily. Right now this only uses the AZ300 gene list. It also needs to include Sean's cancer predisposition list and create proper Mini-BAMs and VCFs that include regions with +/- 1kb around all somatic SNVs, CNVs and SVs.

```
mkdir igv
cd igv

# Always add a list of known cancer genes even where 
# no mutations are found
cp ${BCINSTALL}/genomes/Hsapiens/GRCh37/coverage/prioritize/cancer/az300.bed.gz .
gunzip az300.bed.gz

# Add Ensemble calls
cut -f 1-3 az300.bed > ${BCTUMOR}.bed
bcftools view -H ../somatic/${BCBATCH}-ensemble-annotated.vcf.gz | awk -v OFS="\t" '{print $1, $2-100, $2+100}' >> ${BCTUMOR}.bed

# Bring in conservative Manta calls
cat ../structural/${BCBATCH}_manta_regions_merged.bed >> ${BCTUMOR}.bed

# Sort and merge
bedtools sort -i ${BCTUMOR}.bed > ${BCTUMOR}_sorted.bed
bedtools merge -i ${BCTUMOR}_sorted.bed > ${BCTUMOR}_sorted_merged.bed 

sambamba view -f bam -L ./${BCTUMOR}_sorted_merged.bed ${BCRESULT}/${BCNORMAL}/${BCNORMAL}-ready.bam -t 28 -o ${BCNORMAL}_mini.bam
samtools index ${BCNORMAL}_mini.bam

sambamba view -f bam -L ./${BCTUMOR}_sorted_merged.bed {BCRESULT}/${BCTUMOR}/${BCTUMOR}-ready.bam -t 28 -o ${BCTUMOR}_mini.bam
samtools index ${BCTUMOR}_mini.bam

cd ..
```

## Additional information

```
mkdir log
cd log
ln ${BCRESULT}${BCFINAL}/data_versions.csv .
ln ${BCRESULT}${BCFINAL}/programs.txt .
ln ${BCRESULT}/config/* .
```



