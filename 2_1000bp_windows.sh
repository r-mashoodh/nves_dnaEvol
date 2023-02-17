#!/bin/bash

# create pileup
samtools mpileup -B ../map/F1.qf.rmd.sort.bam > F1.pileup
samtools mpileup -B ../map/N1.qf.rmd.sort.bam > N1.pileup
samtools mpileup -B ../map/F2.qf.rmd.sort.bam > F2.pileup
samtools mpileup -B ../map/N2.qf.rmd.sort.bam > N2.pileup


 filtering indels
perl ../software/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --indel-window 5 --min-count 2 --input F1.pileup --output F1.indels.gtf
perl ../software/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --indel-window 5 --min-count 2 --input N1.pileup --output N1.indels.gtf
perl ../software/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --indel-window 5 --min-count 2 --input F2.pileup --output F2.indels.gtf
perl ../software/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --indel-window 5 --min-count 2 --input N2.pileup --output N2.indels.gtf

perl ../software/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input F1.pileup --gtf F1.indels.gtf --output F1.idf.pileup
perl ../software/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input N1.pileup --gtf N1.indels.gtf --output N1.idf.pileup
perl ../software/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input F2.pileup --gtf F2.indels.gtf --output F2.idf.pileup
perl ../software/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input N2.pileup --gtf N2.indels.gtf --output N2.idf.pileup


##CALCULATING TAJIMA’S / PI DIVERSITY
## –min-coverage –max-coverage: for subsampled files not important; should contain target coverage, i.e.: 10
## –min-covered-fraction minimum percentage of sites having sufficient coverage in the given window
## –min-count minimum occurrence of allele for calling a SNP
## –measure which population genetics measure should be computed (pi/theta/D)
## –pool-size number of chromosomes (thus number of diploids times two)
## –region compute the measure only for a small region; default is the whole genome
## –output a file containing the measure () for the windows
## –snp-output a file containing for every window the SNPs that have been used for computing the measure (e.g.)
## –window-size –step-size control behaviour of sliding window; if step size is smaller than window size than the windows will be overlapping.

perl ../software/popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N1.idf.pileup --min-count 2 --min-coverage 40 --max-coverage 500 --min-covered-fraction 0.25 --pool-size 104 --window-size 1000 --step-size 1000 --output N1.1000.pi --snp-output N1.1000.snps
perl ../software/popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input F2.idf.pileup --min-count 2 --min-coverage 40 --max-coverage 500 --min-covered-fraction 0.25 --pool-size 104 --window-size 1000 --step-size 1000 --output F2.1000.pi --snp-output F2.1000.snps
perl ../software/popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N2.idf.pileup --min-count 2 --min-coverage 40 --max-coverage 500 --min-covered-fraction 0.25 --pool-size 118 --window-size 1000 --step-size 1000 --output N2.1000.pi --snp-output N2.1000.snps
perl ../software/popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input F1.idf.pileup --min-count 2 --min-coverage 40 --max-coverage 500 --min-covered-fraction 0.25 --pool-size 82 --window-size 1000 --step-size 1000 --output F1.1000.pi --snp-output F1.1000.snps

perl ../software/popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure theta --input N1.idf.pileup --min-count 2 --min-coverage 40 --max-coverage 500 --min-covered-fraction 0.25 --pool-size 104 --window-size 1000 --step-size 1000 --output N1.1000.theta --snp-output N1.1000.snps
perl ../software/popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure theta --input F2.idf.pileup --min-count 2 --min-coverage 40 --max-coverage 500 --min-covered-fraction 0.25 --pool-size 104 --window-size 1000 --step-size 1000 --output F2.1000.theta --snp-output F2.1000.snps
perl ../software/popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure theta --input N2.idf.pileup --min-count 2 --min-coverage 40 --max-coverage 500 --min-covered-fraction 0.25 --pool-size 118 --window-size 1000 --step-size 1000 --output N2.1000.theta --snp-output N2.1000.snps
perl ../software/popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure theta --input F1.idf.pileup --min-count 2 --min-coverage 40 --max-coverage 500 --min-covered-fraction 0.25 --pool-size 82 --window-size 1000 --step-size 1000 --output F1.1000.theta --snp-output F1.1000.snps


