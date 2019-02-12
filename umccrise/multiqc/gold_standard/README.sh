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
