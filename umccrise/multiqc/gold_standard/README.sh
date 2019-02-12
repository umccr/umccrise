#############################
### Prepare gold standard ###
#############################
# ssh spa

cd /data/cephfs/punim0010/data/Results/Tothill-A5/2018-08-11/making_background_qc

DIR=umccrised_2019
# Run recent version of umccrise:
umccrise ../final -c -j30 -o ${DIR} multiqc -s E201,E199,E194,E190,E202

# Subset filelist to QC files only:
DIR_QCONLY=${DIR}.qconly
mkdir ${DIR_QCONLY}
for batch in E201 E199 E194 E190 E202 ; do
    cat ${DIR}/work/${batch}_*/multiqc_data/filelist.txt | grep -v gold_standard | grep ${batch}\
        >> ${DIR_QCONLY}/background_multiqc_filelist.txt
done

# Subset output to QC files only:
for fpath in $(cat ${DIR_QCONLY}/background_multiqc_filelist.txt) ; do
    new_path=${DIR_QCONLY}/${fpath}
    mkdir -p $(dirname ${new_path})
    cp -r ${DIR}/${fpath} ${new_path}
done

# Rename sample IDs in all files and filenames into target nice names (Alice, Bob, ...)
python rename.py ${DIR_QCONLY} ${DIR_QCONLY}.renamed




############## OLD ##############

#### Inititate the subset final directory
#mkdir final.subset
#for f in $(cat final/2018-08-11_2018-07-31T0005_Tothill-A5_WGS-merged/multiqc/list_files_final.txt | grep -v trimmed | grep -v target_info.yaml) ; do
#    mkdir -p final.subset/`dirname $f`
#    cp final/$f final.subset/$f
#done
#
#cat final/2018-08-11_2018-07-31T0005_Tothill-A5_WGS-merged/multiqc/list_files_final.txt |\
#    grep -P "E201|E199|E194|E190|E202" |\
#    grep -v "indexcov.tsv" |\
#    grep -v "Per_base_N_content.tsv" |\
#    grep -v "Per_base_sequence_content.tsv" |\
#    grep -v "Per_base_sequence_quality.tsv" |\
#    grep -v "Per_sequence_GC_content.tsv" |\
#    grep -v "Per_sequence_quality_scores.tsv" |\
#    grep -v "Per_tile_sequence_quality.tsv" |\
#    grep -v "Sequence_Length_Distribution.tsv" |\
#    grep -v "sort-chr.qsig.vcf" |\
#    grep -v "ped_check.rel-difference.csv" |\
#    grep -v ".html" |\
#    grep -v "verifybamid" |\
#    grep -v "qsignature" |\
#    > final.subset/list_files_final.txt
#
#### Clean up ###
#cd final.subset
#
## Remove sample dirs:
#ls | grep -v -P "E201|E199|E194|E190|E202|2018-08-11_2018-07-31T0005_Tothill-A5_WGS-merged|list_files_final.txt" | xargs rm -rf
#
## Clean up qc/coverage
#find . -name "*-indexcov.tsv" -delete
## Clean up qc/fastqc
#find . -name "fastqc_report.html" -delete
#find . -name "Per_base_N_content.tsv" -delete
#find . -name "Per_base_sequence_content.tsv" -delete
#find . -name "Per_base_sequence_quality.tsv" -delete
#find . -name "Per_sequence_GC_content.tsv" -delete
#find . -name "Per_sequence_quality_scores.tsv" -delete
#find . -name "Per_tile_sequence_quality.tsv" -delete
#find . -name "Sequence_Length_Distribution.tsv" -delete
## Clean up qc/qsignature
#rm -rf */qc/qsignature */mixup_check
## Clean up qc/peddy -  TODO: remove peddy and VerifyBAMID
#find . -path "*qc/peddy/*.ped_check.rel-difference.csv" -delete
#find . -path "*qc/peddy/*.html" -delete
## Clean up bcbio metrics
#cd 2018-08-11_2018-07-31T0005_Tothill-A5_WGS-merged/multiqc/report/metrics
#ls | grep -v -P "E201|E199|E194|E190|E202" | xargs rm
#cd ../../../../
#
#cd ..
#### Rename ###
#python rename.py final.subset final.subset.renamed
#
#### Run conpair to complete MultiQC report ###
#conpair -T final/*E194-T*/*-ready.bam -N final/*E194-B*/*-ready.bam -tn Alice_T  -nn Alice_B  -o final.subset.renamed/conpair &
#conpair -T final/*E201-T*/*-ready.bam -N final/*E201-B*/*-ready.bam -tn Bob_T    -nn Bob_B    -o final.subset.renamed/conpair &
#conpair -T final/*E199-T*/*-ready.bam -N final/*E199-B*/*-ready.bam -tn Chen_T   -nn Chen_B   -o final.subset.renamed/conpair &
#conpair -T final/*E190-T*/*-ready.bam -N final/*E190-B*/*-ready.bam -tn Dakota_T -nn Dakota_B -o final.subset.renamed/conpair &
#conpair -T final/*E202-T*/*-ready.bam -N final/*E202-B*/*-ready.bam -tn Eugene_T -nn Eugene_B -o final.subset.renamed/conpair &
#cd final.subset.renamed
#rm conpair/*.pileup
#ls conpair/concordance/* >> list_files_final.txt
#ls conpair/contamination/* >> list_files_final.txt

