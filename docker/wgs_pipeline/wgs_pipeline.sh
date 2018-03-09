mkdir $OUT_DIR/normal
mkdir $OUT_DIR/tumor
mkdir $OUT_DIR/tmp/
ln -s $DATA_FOLDER /home/
ln -s $NORMAL_BAM /home/input/normal.bam
ln -s $TUMOR_BAM /home/input/tumor.bam
ln -s $CONFIG_NORMAL /home/input/bic-seq-normal.config
ln -s $CONFIG_TUMOR /home/input/bic-seq-tumor.config
ln -s $CONFIG /home/input/bic-seq.config
cd /home/
echo "Printing paths"
ls /home/ 
ls /home/input/
perl -p -i -e "s/REPLACE/$OUT_DIR/" /home/input/bic-seq-normal.config
perl -p -i -e "s/REPLACE/$OUT_DIR/" /home/input/bic-seq-tumor.config
perl -p -i -e "s/REPLACE/$OUT_DIR/" /home/input/bic-seq.config
/home/samtools-0.1.7a_getUnique-0.1.3/samtools view -U BWA,$OUT_DIR/normal/,N,N /home/input/normal.bam
/home/samtools-0.1.7a_getUnique-0.1.3/samtools view -U BWA,$OUT_DIR/tumor/,N,N /home/input/tumor.bam
/home/NBICseq-norm_v0.2.4/NBICseq-norm.pl --tmp=$OUT_DIR/tmp/ /home/input/bic-seq-normal.config  /home/output/normal/gam.model.txt
/home/NBICseq-norm_v0.2.4/NBICseq-norm.pl --tmp=$OUT_DIR/tmp/ /home/input/bic-seq-tumor.config  $OUT_DIR/tumor/gam.model.txt
/home/NBICseq-seg_v0.7.2/NBICseq-seg.pl --tmp=$OUT_DIR/tmp/ --control /home/input/bic-seq.config  $OUT_DIR/cnv.csv
tar -zcf $OUTPUT_FILE $OUT_DIR
cp $OUTPUT_FILE /home/output
mv /home/output $OUT_DIR