#!/bin/bash

# make a samples.txt
for f in *_R1_001.fastq.gz; do echo ${f%_R*}; done > samples.txt

# run the loop to trim off parts of reads that go past reverse complements of opposite primers
for sample in $(cat samples.txt)
do
    echo "On sample: $sample"
    cutadapt -a ^TCGATGAAGAACGCAGCG...GCATATCAATAAGCGGAGGA \
    -A ^TCCTCCGCTTATTGATATGC...CGCTGCGTTCTTCATCGA \
    -m 215 -M 285 --discard-untrimmed \
    -o ${sample}_R1_001_t.fastq.gz -p ${sample}_R2_001_t.fastq.gz \
    ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz \
    >> cutadapt_primer_trimming_stats.txt 2>&1
done

# output another file showing reads and bp that passed the trimming step
paste samples.txt <(grep "passing" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") > trim_info.txt
sed -i  '1i sample_name,reads_retained,bp_retained' trim_info.txt

# put raw reads in a separate directory
mkdir raw_reads
mv *_001.fastq.gz raw_reads

 
