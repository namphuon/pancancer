#First filter out non-good reads
samtools flagstat ${BAM} > ${OUTPUT_DIR}/${TYPE}.flagstat
python /home/read_count.py -b ${BAM} -s ${BED} -o ${OUTPUT_DIR}/${TYPE}.csv