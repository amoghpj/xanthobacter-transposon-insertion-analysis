echo "V0C4" 
infile="V0C4_S2_R1_processed"
readsdir="processed"

bowtie2\
    -x reference/V0C4\
    -f -U $readsdir/$infile.fasta\
    -S $readsdir/"$infile"_mapped.sam --no-unal


echo "VSF"
infile="VSF_S3_R1_processed"
bowtie2\
    -x reference/V0C4\
    -f -U $readsdir/$infile.fasta\
    -S $readsdir/"$infile"_mapped.sam --no-unal
