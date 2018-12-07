#!/bin/bash

#################################################################################################################################################
#This is a quick ATAC pipeline to get your analysis going, it's ideal if your just starting out with ATAC data and do not know what to do!!	#
#It will align, filter and merge your .fastq files, you can then view the aligned .bams in a genome browser.					#
#I've set it up so that it should run reasnobly well on an average laptop. (Quad core, at least 8gb RAM)					#
#It will take all fastq files from 4 lanes of an illumina nextseq, for paired end data (I plan to add an option for non-paied end)		#
#It works well just running in the background by setting the threads to run at around 4. (On a quad core processor)				#
#Or if you want it to run as quick as possble, run it on all the threads you have and leave it alone for a bit					#
#Dependencies: The following must be installed: Bowtie2, Samtools, Bamtools									#
#Input: Takes compressed .fastq files, and path to bowtie2 genome										#
#Make sure the file is set to executable													#
#################################################################################################################################################

echo -n "Enter the path where the compressed .fastq files are located>"
read path

echo -n "Enter lane 1 forward alignment file name>"
read forward1

echo -n "Enter lane 1 reverse alignment file name>"
read reverse1

echo -n "Enter lane 2 forward alignment file name>"
read forward2

echo -n "Enter lane 2 reverse alignment file name>"
read reverse2

echo -n "Enter lane 3 forward alignment file name>"
read forward3

echo -n "Enter lane 3 reverse alignment file name>"
read reverse3

echo -n "Enter lane 4 forward alignment file name>"
read forward4

echo -n "Enter lane 4 reverse alignment file name>"
read reverse4

echo -n "Enter number of threads to run >"
read threads

echo -n "Enter file path for bowtie2 genome >"
read genome

echo -n "Enter output file name >"
read output

workingdir=$(pwd)

bowtie2 -p $threads -x $genome -1 $path/$forward1 -2 $path/$reverse1 | samtools view -bS - > alignment_temp.bam

echo "Alignment 1 completed, will now filter alignment 1 and start alignment 2......"

bowtie2 -p $threads -x $genome -1 $path/$forward2 -2 $path/$reverse2 | samtools view -bS - > alignment2_temp.bam & 

bamtools index -in alignment_temp.bam

bamtools sort -in alignment_temp.bam -out sorted_temp.bam

bamtools index -in sorted_temp.bam

mkdir ATAC_raw_alignment
mv alignment_temp.bam $workingdir/ATAC_raw_alignment

samtools view -b sorted_temp.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chX chY> filtered_temp.bam

rm sorted_temp.bam

echo "Mitochondrial DNA alignment 1 filtered"

bamtools filter -isDuplicate false -in filtered_temp.bam -out filtered2_temp.bam

echo "Duplicates Mitochondrial DNA alignment 1 filtered"

bamtools filter -isMapped true -in filtered2_temp.bam -out filtered_temp.bam

echo "non mapped alignment 1 filtered"  

bamtools filter -isPaired true -in filtered_temp.bam -out filtered2_temp.bam

echo "non paired alignment 1 filtered"

bamtools filter -isPrimaryAlignment true -in filtered2_temp.bam -out filteredcomplete1.bam

rm filtered2_temp.bam
rm filtered_temp.bam

wait

echo "Alignment 2 Complete and alignment 1 has been filtered, will now filter alignment 2 and start alignment 3"

bowtie2 -p $threads -x $genome -1 $path/$forward3 -2 $path/$reverse3 | samtools view -bS - > alignment3_temp.bam &

bamtools index -in alignment2_temp.bam

bamtools sort -in alignment2_temp.bam -out sorted_temp.bam

bamtools index -in sorted_temp.bam

mv alignment2_temp.bam $workingdir/ATAC_raw_alignment

samtools view -b sorted_temp.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chX chY> filtered_temp.bam

rm sorted_temp.bam

echo "Mitochondrial DNA alignment 2 filtered"

bamtools filter -isDuplicate false -in filtered_temp.bam -out filtered2_temp.bam

echo "Duplicates Mitochondrial DNA alignment 2 filtered"

bamtools filter -isMapped true -in filtered2_temp.bam -out filtered_temp.bam

echo "non mapped alignment 2 filtered"  

bamtools filter -isPaired true -in filtered_temp.bam -out filtered2_temp.bam

echo "non paired alignment 2 filtered"

bamtools filter -isPrimaryAlignment true -in filtered2_temp.bam -out filteredcomplete2.bam

wait

echo "Alignment 3 complete and alignment 2 has been filtered, will now filter alignment 3 and start alignment 4"

bowtie2 -p $threads -x $genome -1 $path/$forward4 -2 $path/$reverse4 | samtools view -bS - > alignment4_temp.bam &

bamtools index -in alignment3_temp.bam

bamtools sort -in alignment3_temp.bam -out sorted_temp.bam

bamtools index -in sorted_temp.bam

mv alignment3_temp.bam $workingdir/ATAC_raw_alignment

samtools view -b sorted_temp.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chX chY > filtered_temp.bam

rm sorted_temp.bam

echo "Mitochondrial DNA alignment 3 filtered"

bamtools filter -isDuplicate false -in filtered_temp.bam -out filtered2_temp.bam

echo "Duplicates Mitochondrial DNA alignment 3 filtered"

bamtools filter -isMapped true -in filtered2_temp.bam -out filtered_temp.bam

echo "non mapped alignment 3 filtered"  

bamtools filter -isPaired true -in filtered_temp.bam -out filtered2_temp.bam

echo "non paired alignment 3 filtered"

bamtools filter -isPrimaryAlignment true -in filtered2_temp.bam -out filteredcomplete3.bam

wait

echo "Alignment 4 complete and Alignment 3 has been filtered, will now filter alignment 4"

bamtools index -in alignment4_temp.bam

bamtools sort -in alignment4_temp.bam -out sorted_temp.bam

bamtools index -in sorted_temp.bam

mv alignment4_temp.bam $workingdir/ATAC_raw_alignment

samtools view -b sorted_temp.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chX chY> filtered_temp.bam

rm sorted_temp.bam

echo "Mitochondrial DNA alignment 4 filtered"

bamtools filter -isDuplicate false -in filtered_temp.bam -out filtered2_temp.bam

echo "Duplicates Mitochondrial DNA alignment 4 filtered"

bamtools filter -isMapped true -in filtered2_temp.bam -out filtered_temp.bam

echo "non mapped alignment 4 filtered"  

bamtools filter -isPaired true -in filtered_temp.bam -out filtered2_temp.bam

echo "non paired alignment 4 filtered"

bamtools filter -isPrimaryAlignment true -in filtered2_temp.bam -out filteredcomplete4.bam

echo "Alignment 4 complete, will now Merge bam files...."

bamtools merge -in filteredcomplete1.bam -in filteredcomplete2.bam -in filteredcomplete3.bam -in filteredcomplete4.bam -out $output

bamtools index -in $output

rm filteredcomplete1.bam
rm filteredcomplete2.bam
rm filteredcomplete3.bam
rm filteredcomplete4.bam

