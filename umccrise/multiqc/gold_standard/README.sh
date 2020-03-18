#####################################
### Prepare gold standard (bcbio) ###
#####################################
# ssh spa

cd /data/cephfs/punim0010/data/Results/Tothill-A5/2018-08-11

DIR=umccrised
# Run recent version of umccrise:
umccrise final -c -j30 -o ${DIR} multiqc -s E201,E199,E194,E190,E202

cd making_background_qc

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
python rename.py ${DIR_QCONLY} ${DIR_QCONLY}.renamed ${DIR_QCONLY}.renamed_hg38



######################################
### Prepare gold standard (DRAGEN) ###
######################################
# ssh rjn

cd /g/data/gx8/extras/umccrise_multiqc_background
iap files download gds://"umccr-primary-data-dev/results/*" .

umccrise Alice -c -j30 -o umccrised/Alice multiqc
umccrise Bob -c -j30 -o umccrised/Bob multiqc
umccrise Chen -c -j30 -o umccrised/Chen multiqc
umccrise Dakota -c -j30 -o umccrised/Dakota multiqc
umccrise Elon -c -j30 -o umccrised/Elon multiqc
















