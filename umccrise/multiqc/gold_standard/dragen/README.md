## Preparing MultiQC background samples for DRAGEN runs

First, we are running 5 TES tasks to process the samples with Dragen
(see below "### TES definitions and launch jsons").

Then, we connect to gadi, pull the Dragen results, and extract the QC files into a new folder:

```
# ssh rjn
cd /g/data/gx8/extras/umccrise_multiqc_background
iap files download gds://"umccr-primary-data-dev/results/*" .

DIR=umccrised
mkdir ${DIR}
umccrise Alice  -c -j30 -o ${DIR}/Alice  multiqc
umccrise Bob    -c -j30 -o ${DIR}/Bob    multiqc
umccrise Chen   -c -j30 -o ${DIR}/Chen   multiqc
umccrise Dakota -c -j30 -o ${DIR}/Dakota multiqc
umccrise Elon   -c -j30 -o ${DIR}/Elon   multiqc

DIR_QCONLY=${DIR}.qconly_hg38
mkdir ${DIR_QCONLY}
echo "" > ${DIR_QCONLY}/background_multiqc_filelist.txt
for batch in Alice Bob Chen Dakota Elon ; do
    # copy the QC files themselves too:
    mkdir ${DIR_QCONLY}/${batch}
    cd ${DIR}/${batch}
    for fpath in $(cat work/${batch}/multiqc_data/filelist.txt | grep -v gold_standard) ; do
        cp -r ${fpath} ../../${DIR_QCONLY}/${batch}
        echo $(readlink -e ../../${DIR_QCONLY}/${batch}/$(basename $fpath)) >> ../../${DIR_QCONLY}/background_multiqc_filelist.txt
    done
    cd ../../
done
```

Then we add the umccrised.qconly folder to the repo (this directory).


### TES definitions and launch jsons

https://aps2.platform.illumina.com/v1/tasks
```
{
    "id": "tdn.6f07a2cfa86941c494bb39c5477b8095",
    "name": "MultiQCBackgroundSamples",
    "href": "https://aps2.platform.illumina.com/v1/tasks/tdn.6f07a2cfa86941c494bb39c5477b8095",
    "description": "MultiQCBackgroundSamples",
    "taskVersions": [],
    "acl": [
        "tid:YXdzLXVzLXBsYXRmb3JtOjEwMDAwNTM3OjBiYTU5YWUxLWZkYWUtNDNiYS1hM2I1LTRkMzY3YTQzYWJkNQ",
        "wid:e4730533-d752-3601-b4b7-8d4d2f6373de"
    ],
    "tenantId": "YXdzLXVzLXBsYXRmb3JtOjEwMDAwNTM3OjBiYTU5YWUxLWZkYWUtNDNiYS1hM2I1LTRkMzY3YTQzYWJkNQ",
    "subTenantId": "wid:e4730533-d752-3601-b4b7-8d4d2f6373de",
    "createdBy": "590dfb6c-6e4f-3db8-9e23-2d1039821653",
    "timeCreated": "2020-03-16T04:54:54.6454838Z",
    "modifiedBy": "590dfb6c-6e4f-3db8-9e23-2d1039821653",
    "timeModified": "2020-03-16T04:54:54.6454838Z"
}
```

https://aps2.platform.illumina.com/v1/tasks/tdn.6f07a2cfa86941c494bb39c5477b8095/versions
```
{
    "version": "2",
    "description": "Somatic tumor-normal analysis through DRAGEN",
    "execution": {
        "image": {
            "name": "699120554104.dkr.ecr.us-east-1.amazonaws.com/public/dragen",
            "tag": "3.5.7"
        },
        "command": "bash",
        "args": [
            "-c",
            "/opt/edico/bin/dragen --partial-reconfig DNA-MAPPER --ignore-version-check true; mkdir -p /ephemeral/ref; tar -C /ephemeral/ref -xvf /mount/index/hg38/hg38_dragen_ht.tar; /opt/edico/bin/dragen --lic-instance-id-location /opt/instance-identity -f --ref-dir /ephemeral/ref --tumor-fastq1 /mount/fastqs/t1.fastq.gz --tumor-fastq2 /mount/fastqs/t2.fastq.gz -1 /mount/fastqs/n1.fastq.gz -2 /mount/fastqs/n2.fastq.gz --RGID {{normal_name}} --RGSM {{normal_name}} --RGID-tumor {{tumor_name}} --RGSM-tumor {{tumor_name}} --output-directory /output/{{prefix}} --output-file-prefix {{prefix}} --sv-reference /mount/index/hg38/hg38.fa --enable-map-align true --output-format BAM --enable-map-align-output true --enable-variant-caller true --enable-duplicate-marking true --vc-enable-clustered-events-filter true --vc-enable-triallelic-filter true --enable-sv true"
        ],
        "inputs": [
            {
                "Path": "/mount/index/hg38/hg38_dragen_ht.tar",
                "Url": "gds://umccr-refdata-dev/dragen/genomes/hg38/hg38_dragen_ht.tar",
                "mode": "download"
            },
            {
                "Path": "/mount/index/hg38/hg38.fa",
                "Url": "gds://umccr-refdata-dev/dragen/genomes/hg38/hg38.fa",
                "mode": "download"
            },
            {
                "Path": "/mount/index/hg38/hg38.fa.fai",
                "Url": "gds://umccr-refdata-dev/dragen/genomes/hg38/hg38.fa.fai",
                "mode": "download"
            },
            {
                "Path": "/mount/fastqs/t1.fastq.gz",
                "Url": "{{fq_t1}}",
                "mode": "download"
            },
            {
                "Path": "/mount/fastqs/t2.fastq.gz",
                "Url": "{{fq_t2}}",
                "mode": "download"
            },
            {
                "Path": "/mount/fastqs/n1.fastq.gz",
                "Url": "{{fq_n1}}",
                "mode": "download"
            },
            {
                "Path": "/mount/fastqs/n2.fastq.gz",
                "Url": "{{fq_n2}}",
                "mode": "download"
            }
        ],
        "outputs": [
            {
                "path": "/output/{{prefix}}",
                "url": "{{results_folder}}/{{prefix}}"
            },
            {
                "path": "/var/log/tessystemlogs",
                "url": "{{results_folder}}/{{prefix}}/dragenlogs"
            }
        ],
        "environment": {
            "resources": {
                "Type": "fpga",
                "Size": "medium"
            }
        }
    }
}
```

//Alice_T
//gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E194_PRJ180506_Missing_R1_001.fastq.gz
//gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E194_PRJ180507_Missing_R1_001.fastq.gz
https://aps2.platform.illumina.com/v1/tasks/tdn.db193bad9280403fae32ed18cd5872cf/versions/tvn.777ca152d6a0481e877d3f711817edc8:launch
```
{
    "name": "Alice",
    "arguments": {
        "prefix": "Alice",
        "results_folder": "gds://umccr-primary-data-dev/multiqc_controls",
        "fq_t1": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E194_PRJ180506_Missing_R1_001.fastq.gz",
        "fq_t2": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E194_PRJ180506_Missing_R2_001.fastq.gz",
        "fq_n1": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E194_PRJ180507_Missing_R1_001.fastq.gz",
        "fq_n2": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E194_PRJ180507_Missing_R2_001.fastq.gz",
        "normal_name": "Alice_N",
        "tumor_name": "Alice_T"
    }
 }
```

//Bob_T
//gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E201_PRJ180492_Missing_R1_001.fastq.gz
//gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E201_PRJ180493_Missing_R1_001.fastq.gz
https://aps2.platform.illumina.com/v1/tasks/tdn.db193bad9280403fae32ed18cd5872cf/versions/tvn.777ca152d6a0481e877d3f711817edc8:launch
```
{
    "name": "Bob",
    "arguments": {
        "prefix": "Bob",
        "results_folder": "gds://umccr-primary-data-dev/multiqc_controls",
        "fq_t1": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E201_PRJ180492_Missing_R1_001.fastq.gz",
        "fq_t2": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E201_PRJ180492_Missing_R2_001.fastq.gz",
        "fq_n1": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E201_PRJ180493_Missing_R1_001.fastq.gz",
        "fq_n2": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E201_PRJ180493_Missing_R2_001.fastq.gz",
        "normal_name": "Bob_N",
        "tumor_name": "Bob_T"
    }
 }
```

//Chen
//gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E199_PRJ180494_Missing_R1_001.fastq.gz
//gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E199_PRJ180495_Missing_R1_001.fastq.gz
https://aps2.platform.illumina.com/v1/tasks/tdn.db193bad9280403fae32ed18cd5872cf/versions/tvn.777ca152d6a0481e877d3f711817edc8:launch
```
{
    "name": "Chen",
    "arguments": {
        "prefix": "Chen",
        "results_folder": "gds://umccr-primary-data-dev/multiqc_controls",
        "fq_t1": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E199_PRJ180494_Missing_R1_001.fastq.gz",
        "fq_t2": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E199_PRJ180494_Missing_R2_001.fastq.gz",
        "fq_n1": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E199_PRJ180495_Missing_R1_001.fastq.gz",
        "fq_n2": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180718_A00130_0068_BH5M73DSXX_E199_PRJ180495_Missing_R2_001.fastq.gz",
        "normal_name": "Chen_N",
        "tumor_name": "Chen_T"
    }
 }
```

//Dakota
//gds://umccr-fastq-data-dev/PD/multiqc_controls/180521_A00130_0060_AH5CHKDSXX_E190_PRJ180253_Missing_R1_001.fastq.gz
//gds://umccr-fastq-data-dev/PD/multiqc_controls/180521_A00130_0060_AH5CHKDSXX_E190_PRJ180254_Missing_R1_001.fastq.gz
https://aps2.platform.illumina.com/v1/tasks/tdn.db193bad9280403fae32ed18cd5872cf/versions/tvn.777ca152d6a0481e877d3f711817edc8:launch
```
{
    "name": "Dakota",
    "arguments": {
        "prefix": "Dakota",
        "results_folder": "gds://umccr-primary-data-dev/multiqc_controls",
        "fq_t1": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180521_A00130_0060_AH5CHKDSXX_E190_PRJ180253_Missing_R1_001.fastq.gz",
        "fq_t2": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180521_A00130_0060_AH5CHKDSXX_E190_PRJ180253_Missing_R2_001.fastq.gz",
        "fq_n1": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180521_A00130_0060_AH5CHKDSXX_E190_PRJ180254_Missing_R1_001.fastq.gz",
        "fq_n2": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180521_A00130_0060_AH5CHKDSXX_E190_PRJ180254_Missing_R2_001.fastq.gz",
        "normal_name": "Dakota_N",
        "tumor_name": "Dakota_T"
    }
 }

//Elon
//gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E202_PRJ180499_Missing_R1_001.fastq.gz
//gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E202_PRJ180500_Missing_R1_001.fastq.gz
https://aps2.platform.illumina.com/v1/tasks/tdn.db193bad9280403fae32ed18cd5872cf/versions/tvn.777ca152d6a0481e877d3f711817edc8:launch
{
    "name": "Elon",
    "arguments": {
        "prefix": "Elon",
        "results_folder": "gds://umccr-primary-data-dev/multiqc_controls",
        "fq_t1": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E202_PRJ180499_Missing_R1_001.fastq.gz",
        "fq_t2": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E202_PRJ180499_Missing_R2_001.fastq.gz",
        "fq_n1": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E202_PRJ180500_Missing_R1_001.fastq.gz",
        "fq_n2": "gds://umccr-fastq-data-dev/PD/multiqc_controls/180720_A00130_0070_AH5TY5DSXX_E202_PRJ180500_Missing_R2_001.fastq.gz",
        "normal_name": "Elon_N",
        "tumor_name": "Elon_T"
    }
 }
```






