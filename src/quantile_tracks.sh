
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


function paste_together() {
    CHROM=$1

    #load bgzip
    module load htslib

    # read bams
    readarray BAMS < ~/bam_list.txt
    
    # get sample and file names
    i=0
    for BAM in ${BAMS[*]}
    do
        NAMES[i]=`basename ${BAM/.trimmed.bowtie2.filtered.bam/}`
        FILENAMES[i]=~/coverage/${NAMES[i]}_${CHROM}.coverage.bed 
        let i=i+1
    done

    paste ~/windows/${CHROM}_position.bed ${FILENAMES[*]} | bgzip > ~/coverage/all_chr_${CHROM}.bed.gz
}

export -f paste_together


function quantilize() {
    CHROM=$1

    zcat all_chr_${CHROM}.bed.gz | \
    python ~/quantilize.py ~/coverage/all_chr_${CHROM}
}

export -f quantilize


function concatenate_chroms() {
    SUFFIX=$1
    cat ~/coverage/all_chr_{1..22}_${SUFFIX}.bed > all_${SUFFIX}.bedgraph
}

export -f concatenate_chroms


function bedgraph_to_bigwig() {
    BEDGRAPH=$1
    BIGWIG=${BEDGRAPH/bedgraph/bigwig}
    
    bedGraphToBigWig $BEDGRAPH ~/resources/genomes/hg19/hg19.chromSizes $BIGWIG
}

export -f bedgraph_to_bigwig


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
parallel paste_together ::: {1..22}


# get quantiles
parallel quantilize ::: {1..22}


# concat chromosomes
parallel concatenate_chroms ::: quant5 quant25 mean quant75 quant95


# make bigwigs
bedgraph_to_bigwig ::: all_{quant5,quant25,mean,quant75,quant95}.bedgraph
