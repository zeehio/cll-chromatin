
function make_windows() {
    CHROM=$1
    COUNT=`grep -P "^chr${CHROM}\t" resources/genomes/hg19/hg19.chrom.sizes | cut -f 2`
    echo "chr$CHROM 1 $COUNT" | sed 's/ /\t/g' > ~/windows/tmp_${CHROM}.bed
}
export -f make_windows


function windows_position() {
    CHROM=$1
    bedtools makewindows -b ~/windows/tmp_${CHROM}.bed -w 1 | cut -f 1,2,3 > ~/windows/${CHROM}_position.bed
}
export -f windows_position


function get_coverage() {
    CHROM=$1
    BAM=$2
    NAME=`basename ${BAM/.trimmed.bowtie2.filtered.bam/}`

    echo "#! /bin/bash" > job.sh
    echo "#SBATCH --partition=shortq" >> job.sh
    echo "#SBATCH --ntasks=1" >> job.sh

    echo "#SBATCH --cpus-per-task=2" >> job.sh
    echo "#SBATCH --mem-per-cpu=4000" >> job.sh

    echo "#SBATCH --job-name=${NAME}_chr${CHROM}" >> job.sh
    echo "#SBATCH --output=/home/arendeiro/coverage/${NAME}_${CHROM}.log" >> job.sh

    echo "hostname" >> job.sh
    echo "date" >> job.sh

    echo "bedtools coverage -d -abam $BAM -b ~/windows/tmp_${CHROM}.bed | cut -f 5 > ~/coverage/${NAME}_${CHROM}.coverage.bed" >> job.sh

    echo "date" >> job.sh
    sbatch job.sh
}

export -f get_coverage

# Start
# make 1bp windows for each chromosome
parallel make_windows ::: {1..22}

# get bam files
readarray BAMS < ~/bam_list.txt

# get coverage in each chromosome for each sample
for BAM in ${BAMS[*]}
do
    for CHROM in {1..22}
    do
        get_coverage $CHROM $BAM
    done
done

# while you wait for those jobs,
parallel windows_position ::: {1..22}


# concatenate columns across samples for each chromosome
for CHROM in {1..22}
do
    # add "chrom start end" to begining
    paste ~/windows/${CHROM}_position.bed ${BAMS[*]/.trimmed.bowtie2.filtered.bam/_${CHROM}.coverage.bed} | bgzip > ~/coverage/all_chr_${CHROM}.coverage.bgzip
done

# index with tabix
for CHROM in {1..22}
do
    tabix -p bed ~/coverage/all_chr_${CHROM}.coverage.bgzip
done

# try:
# load into R with datatable
# or
# get sections out with tabix, calculate stats, write to file, paste position columns again
