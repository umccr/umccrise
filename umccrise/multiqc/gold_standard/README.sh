#####################################
### Prepare gold standard (bcbio) ###
#####################################
# ssh rjn
cd /g/data/gx8/extras/umccrise_multiqc_background/bcbio

# Run the recent version of umccrise:
DIR=umccrised
source /g/data/gx8/extras/umccrise/load_umccrise.sh
umccrise final -c -j30 -o ${DIR} -T multiqc -s E201 -s E199 -s E194 -s E190 -s E202

# Subset filelist to QC files only:
DIR_QCONLY=${DIR}.qconly
mkdir ${DIR_QCONLY}
touch ${DIR_QCONLY}/background_multiqc_filelist.txt
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

See `dragen/README.sh`















