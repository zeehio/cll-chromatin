
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
    RATIO=$3

    echo "#! /bin/bash" > job.sh
    echo "#SBATCH --partition=shortq" >> job.sh
    echo "#SBATCH --ntasks=1" >> job.sh

    echo "#SBATCH --cpus-per-task=2" >> job.sh
    echo "#SBATCH --mem-per-cpu=4000" >> job.sh

    echo "#SBATCH --job-name=${NAME}_chr${CHROM}" >> job.sh
    echo "#SBATCH --output=/home/arendeiro/coverage/${NAME}_${CHROM}.log" >> job.sh

    echo "hostname" >> job.sh
    echo "date" >> job.sh

    echo "bedtools coverage -d -abam $BAM -b ~/windows/tmp_${CHROM}.bed | cut -f 5 | awk '{print \$1/${RATIO}}' > ~/coverage/${NAME}_${CHROM}.coverage.bed" >> job.sh

    echo "date" >> job.sh
    sbatch job.sh
}

export -f get_coverage


function paste_together() {
    CHROM=$1

    #load bgzip
    module load htslib

    # read bams
    readarray BAMS < /home/arendeiro/cll-patients/data/bamfiles.txt
    
    # get sample and file names
    i=0
    for BAM in ${BAMS[*]}
    do
        NAMES[i]=`basename ${BAM/.trimmed.bowtie2.filtered.bam/}`
        FILENAMES[i]=~/coverage/${NAMES[i]}_${CHROM}.coverage.bed 
        let i=i+1
    done

    echo "#! /bin/bash" > job.sh
    echo "#SBATCH --partition=shortq" >> job.sh
    echo "#SBATCH --ntasks=1" >> job.sh

    echo "#SBATCH --cpus-per-task=2" >> job.sh
    echo "#SBATCH --mem-per-cpu=4000" >> job.sh

    echo "#SBATCH --job-name=${NAME}_chr${CHROM}" >> job.sh
    echo "#SBATCH --output=/home/arendeiro/coverage/${NAME}_${CHROM}.log" >> job.sh

    echo "hostname" >> job.sh
    echo "date" >> job.sh

    echo "paste ~/windows/${CHROM}_position.bed ${FILENAMES[*]} | bgzip > ~/coverage/all_chr_${CHROM}.bed.gz" >> job.sh

    echo "date" >> job.sh
    sbatch job.sh
}

export -f paste_together


function split_files() {
    FILE=$1
    NAME=${FILE/.bed.gz/}
    
    echo "#! /bin/bash" > job.sh
    echo "#SBATCH --partition=shortq" >> job.sh
    echo "#SBATCH --ntasks=1" >> job.sh

    echo "#SBATCH --cpus-per-task=1" >> job.sh
    echo "#SBATCH --mem-per-cpu=4000" >> job.sh

    echo "#SBATCH --job-name=${NAME}_chr${CHROM}" >> job.sh
    echo "#SBATCH --output=/scratch/users/arendeiro/logs/split_files_${NAME}.log" >> job.sh

    echo "hostname" >> job.sh
    echo "date" >> job.sh

    echo "zcat $FILE | split -a 5 -d -l 1000000 - ${NAME}_" >> job.sh

    echo "date" >> job.sh
    sbatch job.sh
}

export -f split_files


function compress(){
    CHROM=$1

    echo "#! /bin/bash" > job.sh
    echo "#SBATCH --partition=shortq" >> job.sh
    echo "#SBATCH --ntasks=1" >> job.sh

    echo "#SBATCH --cpus-per-task=1" >> job.sh
    echo "#SBATCH --mem-per-cpu=2000" >> job.sh

    echo "#SBATCH --job-name=chr${CHROM}" >> job.sh

    echo "hostname" >> job.sh
    echo "date" >> job.sh

    echo "gzip all_chr_${CHROM}_* ">> job.sh

    echo "date" >> job.sh
    sbatch job.sh
}

export -f compress


function quantilize() {
    FILE=$1
    
    echo "#! /bin/bash" > job.sh
    echo "#SBATCH --partition=shortq" >> job.sh
    echo "#SBATCH --ntasks=1" >> job.sh

    echo "#SBATCH --cpus-per-task=1" >> job.sh
    echo "#SBATCH --mem-per-cpu=2000" >> job.sh

    echo "#SBATCH --job-name=quantilize_${FILE}" >> job.sh
    echo "#SBATCH --output=/scratch/users/arendeiro/logs/quantilize_files_${FILE}.log" >> job.sh

    echo "hostname" >> job.sh
    echo "date" >> job.sh

    echo "zcat $FILE | python ~/quantilize.py /fhgfs/groups/lab_bock/arendeiro/coverage/$FILE" >> job.sh

    echo "date" >> job.sh
    sbatch job.sh
}

export -f quantilize


function concatenate_back() {
    SUFFIX=$1
    
    echo "#! /bin/bash" > job.sh
    echo "#SBATCH --partition=longq" >> job.sh
    echo "#SBATCH --ntasks=1" >> job.sh

    echo "#SBATCH --cpus-per-task=4" >> job.sh
    echo "#SBATCH --mem-per-cpu=8000" >> job.sh

    echo "#SBATCH --job-name=concatenate_back_${FILE}" >> job.sh
    echo "#SBATCH --output=/scratch/users/arendeiro/logs/concatenate_back_files_${FILE}.log" >> job.sh

    echo "hostname" >> job.sh
    echo "date" >> job.sh

    echo "cd /fhgfs/groups/lab_bock/arendeiro/coverage/" >> job.sh

    echo "INPUT=(\`ls -v | grep -e all_chr_.*_${SUFFIX} | sort -z\`)" >> job.sh

    echo "cat \${INPUT[@]} > all_${SUFFIX}.bedgraph" >> job.sh

    echo "date" >> job.sh
    sbatch job.sh
}

export -f concatenate_back


function bedgraph_to_bigwig() {
    BEDGRAPH=$1
    BIGWIG=${BEDGRAPH/bedgraph/bigwig}
    
    echo "#! /bin/bash" > job.sh
    echo "#SBATCH --partition=longq" >> job.sh
    echo "#SBATCH --ntasks=1" >> job.sh

    echo "#SBATCH --cpus-per-task=2" >> job.sh
    echo "#SBATCH --mem-per-cpu=100000" >> job.sh

    echo "#SBATCH --job-name=bedgraph_to_bigwig_${BEDGRAPH}" >> job.sh
    echo "#SBATCH --output=/scratch/users/arendeiro/logs/bedgraph_to_bigwig_${BEDGRAPH}.log" >> job.sh

    echo "hostname" >> job.sh
    echo "date" >> job.sh

    echo "cd /fhgfs/groups/lab_bock/arendeiro/coverage/" >> job.sh

    echo "bedGraphToBigWig all_${BEDGRAPH} ~/resources/genomes/hg19/hg19.chromSizes ~/coverage/all_${BIGWIG}" >> job.sh

    echo "date" >> job.sh
    sbatch job.sh
}

export -f bedgraph_to_bigwig


# Start
# make 1bp windows for each chromosome
parallel make_windows ::: {1..22}

# get bam files
readarray BAMS < /home/arendeiro/cll-patients/data/bamfiles.txt
readarray RATIOS < /home/arendeiro/cll-patients/data/bamsizes.txt

# get coverage in each chromosome for each sample
for INDEX in ${!BAMS[@]}
do
    for CHROM in {1..22}
    do
        get_coverage $CHROM ${BAMS[INDEX]} ${RATIOS[INDEX]}
    done
done
# while you wait for those jobs,
parallel windows_position ::: {1..22}


# concatenate columns across samples for each chromosome
for CHROM in {1..22}
do
    paste_together $CHROM
done

# split files in chunks
for CHROM in {1..22}
do
    split_files all_chr_$CHROM.bed.gz
done

# compress file chunks
for CHROM in {1..22}
do
    compress $CHROM
done


# get quantiles and mean
for FILE in `ls ~/coverage/ | grep -v bed | grep -v job`
do
    quantilize $FILE
done


# concat chuncks across chromosomes
for SUFFIX in quant5 quant25 mean quant75 quant95
do
    concatenate_back $SUFFIX
done


# make bigwigs
for SUFFIX in quant5 quant25 mean quant75 quant95
do
    bedgraph_to_bigwig ${SUFFIX}.bedgraph
done
