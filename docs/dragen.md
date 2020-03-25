# Processing validation samples

```
{
    "name": "P033",
    "arguments": {
        "prefix": "P033",
        "results_folder": "gds://umccr-primary-data-dev/FiveValidationSamples",
        "fq_t1": "gds://umccr-fastq-data-dev/FiveValidationSamples/P033/181102_A00130_0081_AHFFWLDSXX_CCR170105_MH17B001P033_R1_001.fastq.gz",
        "fq_t2": "gds://umccr-fastq-data-dev/FiveValidationSamples/P033/181102_A00130_0081_AHFFWLDSXX_CCR170105_MH17B001P033_R2_001.fastq.gz",
        "fq_n1": "gds://umccr-fastq-data-dev/FiveValidationSamples/P033/181102_A00130_0081_AHFFWLDSXX_CCR170115b_MH17T002P033_R1_001.fastq.gz",
        "fq_n2": "gds://umccr-fastq-data-dev/FiveValidationSamples/P033/181102_A00130_0081_AHFFWLDSXX_CCR170115b_MH17T002P033_R2_001.fastq.gz",
        "tumor_name": "P033_T",
        "normal_name": "P033_N"
    }
 }

cd /g/data/gx8/projects/Hofmann_Cromwell/
# P025 tumor - normal
iap files upload 2019-02-01T0241_Cromwell_WGS_2016.249.18.WH.P025/data/180920_A00130_0075_AHC7NCDSXX_CCR180149_VPT_WH025_E_R1_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/P025/P025_T_R1.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_2016.249.18.WH.P025/data/180920_A00130_0075_AHC7NCDSXX_CCR180149_VPT_WH025_E_R2_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/P025/P025_T_R2.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_2016.249.18.WH.P025/data/180920_A00130_0075_AHC7NCDSXX_CCR180135_WH18B002P025_R1_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/P025/P025_T_R1.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_2016.249.18.WH.P025/data/180920_A00130_0075_AHC7NCDSXX_CCR180135_WH18B002P025_R2_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/P025/P025_N_R2.fastq.gz

# P033
iap files upload 2019-02-01T0241_Cromwell_WGS_2016.249.17.MH.P033/data/181102_A00130_0081_AHFFWLDSXX_CCR170115b_MH17T002P033_R1_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/P033/P033_T_R1.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_2016.249.17.MH.P033/data/181102_A00130_0081_AHFFWLDSXX_CCR170115b_MH17T002P033_R2_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/P033/P033_T_R2.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_2016.249.17.MH.P033/data/181102_A00130_0081_AHFFWLDSXX_CCR170105_MH17B001P033_R1_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/P033/P033_T_R1.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_2016.249.17.MH.P033/data/181102_A00130_0081_AHFFWLDSXX_CCR170105_MH17B001P033_R2_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/P033/P033_N_R2.fastq.gz

# FFPE
iap files upload 2019-02-01T0241_Cromwell_WGS_CUP-Pairs8/data/181211_A00130_0084_BHFTGGDSXX_PRJ180660_8_DNA009529_FFPE_R1_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/CUP_Pairs8/CUP_Pairs8_T_R1.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_CUP-Pairs8/data/181211_A00130_0084_BHFTGGDSXX_PRJ180660_8_DNA009529_FFPE_R2_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/CUP_Pairs8/CUP_Pairs8_T_R2.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_CUP-Pairs8/data/181211_A00130_0084_BHFTGGDSXX_PRJ180661_8_DNA008678_Blood_R1_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/CUP_Pairs8/CUP_Pairs8_T_R1.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_CUP-Pairs8/data/181211_A00130_0084_BHFTGGDSXX_PRJ180661_8_DNA008678_Blood_R2_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/CUP_Pairs8/CUP_Pairs8_N_R2.fastq.gz

# SFRC
iap files upload 2019-02-01T0241_Cromwell_WGS_SFRC01073/data/181015_A00130_0078_BHF7T5DSXX_PRJ180598_SFRC01073_S2T_R1_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/SFRC/SFRC_T_R1.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_SFRC01073/data/181015_A00130_0078_BHF7T5DSXX_PRJ180598_SFRC01073_S2T_R2_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/SFRC/SFRC_T_R2.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_SFRC01073/data/181015_A00130_0078_BHF7T5DSXX_PRJ180599_SFRC01073_S2N_R1_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/SFRC/SFRC_T_R1.fastq.gz
iap files upload 2019-02-01T0241_Cromwell_WGS_SFRC01073/data/181015_A00130_0078_BHF7T5DSXX_PRJ180599_SFRC01073_S2N_R2_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/SFRC/SFRC_N_R2.fastq.gz

# B_All
iap files upload 2019-09-11T0258_All_WGS_B_ALL_Case_10/data/190418_A00130_0101_BHKJT3DSXX_MDX190025_B_ALL_Case_10T_R1_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/B_ALL_Case_10/B_ALL_Case_10_T_R1.fastq.gz
iap files upload 2019-09-11T0258_All_WGS_B_ALL_Case_10/data/190418_A00130_0101_BHKJT3DSXX_MDX190025_B_ALL_Case_10T_R2_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/B_ALL_Case_10/B_ALL_Case_10_T_R2.fastq.gz
iap files upload 2019-09-11T0258_All_WGS_B_ALL_Case_10/data/190418_A00130_0101_BHKJT3DSXX_MDX190026_B_ALL_Case_10G_R1_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/B_ALL_Case_10/B_ALL_Case_10_T_R1.fastq.gz
iap files upload 2019-09-11T0258_All_WGS_B_ALL_Case_10/data/190418_A00130_0101_BHKJT3DSXX_MDX190026_B_ALL_Case_10G_R2_001.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/B_ALL_Case_10/B_ALL_Case_10_N_R2.fastq.gz

# SEQC
iap files upload 2019-12-31T0000_SEQC_SEQC50/data/T_SRR7890936_50pc_1.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/SEQC50/SEQC50_T_R1.fastq.gz
iap files upload 2019-12-31T0000_SEQC_SEQC50/data/T_SRR7890936_50pc_2.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/SEQC50/SEQC50_T_R2.fastq.gz
iap files upload 2019-12-31T0000_SEQC_SEQC50/data/N_SRR7890889_1.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/SEQC50/SEQC50_T_R1.fastq.gz
iap files upload 2019-12-31T0000_SEQC_SEQC50/data/N_SRR7890889_2.fastq.gz gds://umccr-fastq-data-dev/FiveValidationSamples/SEQC50/SEQC50_N_R2.fastq.gz
```