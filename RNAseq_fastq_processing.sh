#!/bin/bash
#SBATCH -n 3                    # Number of cores requested
#SBATCH -t 120:00:00             # Runtime in minutes or 
#SBATCH -p medium               # Partition (queue) to submit to
#SBATCH --mem-per-cpu=64G       # 32 GB memory needed (memory PER CORE)
#SBATCH --open-mode=append      # append adds to outfile, truncate deletes first
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file


module load gcc/6.2.0
module load python/2.7.12
module load samtools/1.3.1
module load bedtools/2.26.0
module load star/2.5.4a
module load subread/1.6.0

#alignment of RNAseq fastq files 
for file in *-5_NA_sequence.fastq; do
	out=${file%-5_NA_sequence.fastq}
	STAR --runMode alignReads --genomeDir /n/groups/shared_databases/star_reference/mm10 --readFilesIn ${out}-5_NA_sequence.fastq,${out}-6_NA_sequence.fastq --outFileNamePrefix ${out} --runThreadN 3 --outSAMtype BAM SortedByCoordinate --clip3pNbases 0 --clip5pNbases 3
	samtools index ${out}Aligned.sortedByCoord.out.bam
done

#featurescount to count reads within gene ORFs 
for file in *Aligned.sortedByCoord.out.bam; do 
	out=${file%.bam}
	featureCounts -t gene -g gene_id -s 1 -O -a /n/data2/bch/medicine/kourembanas/rmr15/Reference_Files/Mus_musculus.GRCm38.91.chr.gtf -o GENEfeatures_${out} $file
done 
