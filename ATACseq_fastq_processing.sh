#!/bin/bash
#SBATCH -n 6                    # Number of cores requested
#SBATCH -t 120:00:00             # Runtime in minutes or 
#SBATCH -p medium               # Partition (queue) to submit to
#SBATCH --mem-per-cpu=32G       # 32 GB memory needed (memory PER CORE)
#SBATCH --open-mode=append      # append adds to outfile, truncate deletes first
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file


module load gcc/6.2.0
module load python/2.7.12
module load fastqc/0.11.5
module load bowtie2/2.2.9
module load samtools/1.3.1
module load bedtools/2.26.0
module load macs2/2.1.1.20160309



#bamtofastq
for file in *_mm9_bestmap.bam; do
	out=${file%_mm9_bestmap.bam}
	samtools sort -n $file -o ${out}_mm9.nsorted.bam
	bedtools bamtofastq -i ${out}_mm9.nsorted.bam -fq ${out}_1_sequence.fastq -fq2 ${out}_2_sequence.fastq
	rm ${out}_mm9.nsorted.bam
done

#cleaning of bam files 
for file in *_1_sequence.fastq; do
	out=${file%_1_sequence.fastq}
	bowtie2 -p 6 -x /home/rmr15/mm10_Bowtie2Index/mm10 -X 4000 -1 ${out}_1_sequence.fastq -2 ${out}_2_sequence.fastq -S ${out}.sam
	samtools view -b ${out}.sam | samtools sort - -o ${out}.sorted.bam
	samtools index ${out}.sorted.bam
    samtools view -f 0x2 -b ${out}.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > ${out}.sorted.noMt.bam
	samtools index ${out}.sorted.noMt.bam
	samtools rmdup ${out}.sorted.noMt.bam ${out}.rmdupNOS.sorted.noMt.bam
	samtools index ${out}.rmdupNOS.sorted.noMt.bam
	samtools sort -n ${out}.rmdupNOS.sorted.noMt.bam -o ${out}.nsorted.rmdupNOS.noMt.bam
	bedtools bamtobed -bedpe -i ${out}.nsorted.rmdupNOS.noMt.bam 2>/dev/null | awk 'BEGIN {OFS="\t"}; {if($9 =="+") print $1,$2+4,$3+4,$4,$5-5,$6-5,$7,$8,$9,$10 ; else print $1,$2-5,$3-5,$4,$5+4,$6+4,$7,$8,$9,$10  }' > ${out}.shifted.rmdupNOS.noMt.bed
	bedtools bedpetobam -i ${out}.shifted.rmdupNOS.noMt.bed -g /home/rmr15/UCSCscripts/mm10.chrom.sizes > ${out}.shifted.rmdupNOS.noMt.bam
	samtools sort ${out}.shifted.rmdupNOS.noMt.bam -o ${out}.sorted.shifted.rmdupNOS.noMt.bam
	samtools index ${out}.sorted.shifted.rmdupNOS.noMt.bam
	samtools flagstat ${out}.sorted.shifted.rmdupNOS.noMt.bam > ${out}.sorted.shifted.rmdupNOS.noMt.bam.FLAGSTAT.txt
	rm ${out}.sam
	rm ${out}.sorted.bam
	rm ${out}.sorted.bam.bai
	rm ${out}.sorted.noMt.bam
	rm ${out}.sorted.noMt.bam.bai
	#rm ${out}.rmdupNOS.sorted.noMt.bam
	#rm ${out}.rmdupNOS.sorted.noMt.bam.bai
	rm ${out}.nsorted.rmdupNOS.noMt.bam
	rm ${out}.shifted.rmdupNOS.noMt.bed
	rm ${out}.shifted.rmdupNOS.noMt.bam
done

#ATACseq peak calling
for file in *.sorted.shifted.rmdupNOS.noMt.bam; do 
	out=${file%.sorted.shifted.rmdupNOS.noMt.bam}
	macs2 callpeak -t $file -f BAMPE -g mm -q 0.05 --nomodel --nolambda --keep-dup all -n ${out}_macs2 -B
	macs2 bdgcmp -t ${out}_macs2_treat_pileup.bdg -c ${out}_macs2_control_lambda.bdg -o ${out}_macs2_FE.bdg -m FE
	rm ${out}_macs2_treat_pileup.bdg
	rm ${out}_macs2_control_lambda.bdg
done 

