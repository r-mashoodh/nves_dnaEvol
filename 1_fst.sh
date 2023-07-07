###### Script to go from raw reads to Fst windows

#### Trim fq reads using TrimGalore
./software/TrimGalore/trim_galore --fastqc --illumina --paired -o ./trimmed_reads raw_data/F1_H7HC7DSX2_L2_1.fq raw_data/F1_H7HC7DSX2_L2_2.fq
./software/TrimGalore/trim_galore --fastqc --illumina --paired -o ./trimmed_reads raw_data/F2_H7HC7DSX2_L2_1.fq raw_data/F2_H7HC7DSX2_L2_2.fq
./software/TrimGalore/trim_galore --fastqc --illumina --paired -o ./trimmed_reads raw_data/N1_H7HC7DSX2_L2_1.fq raw_data/N1_H7HC7DSX2_L2_2.fq
./software/TrimGalore/trim_galore --fastqc --illumina --paired -o ./trimmed_reads raw_data/N2_H7HC7DSX2_L2_1.fq raw_data/N2_H7HC7DSX2_L2_2.fq


#### Align to vespilloides genome using bwa
cd trimmed_reads

for i in *_H7HC7DSX2_L2_1_val_1.fq.gz;
do
filename=`echo $i | awk -F '_H7HC7DSX2_L2_1_val_1.fq.gz' '{print $1}'`
# align
bwa mem -t 8 ../ref/nves.fa $filename"_H7HC7DSX2_L2_1_val_1.fq.gz" $filename"_H7HC7DSX2_L2_2_val_2.fq.gz" > ../map/$filename"_aligned.sam"
done

cd ../map

#### Process alignments
for i in *_aligned.sam;
do
filename=`echo $i | awk -F '_aligned.sam' '{print $1}'`
# convert sam to bam
samtools view -Sb map/$filename"_aligned.sam" > map/$filename".bam"
# sort sam -- needs to be via picard for MarkDuplicates
java -jar software/picard.jar SortSam -I map/$filename".bam" -O map/$filename".sort.bam" -VALIDATION_STRINGENCY SILENT -SO coordinate
# mark duplicates
java -jar software/picard.jar MarkDuplicates -I map/$filename".sort.bam" -O map/$filename".rmd.sort.bam" -M map/$filename".dupstat.txt" -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true
# filter poor qual alignments from bam, remove duplicates
samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b map/$filename".rmd.sort.bam" > map/$filename".qf.rmd.sort.bam"
# index bam
samtools index map/$filename".qf.rmd.sort.bam"
done

#### Create pileup
### note that for both glm SNP analysis and popoolation, cols should be in this order F1,N1,F2,N2
samtools mpileup -Q 0 -f ../ref/nves.fa -B F1.qf.rmd.sort.bam N1.qf.rmd.sort.bam F2.qf.rmd.sort.bam N2.qf.rmd.sort.bam > all_pops.pileup
# filtering indels
perl ../software/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --indel-window 5 --min-count 2 --input all_pops.mpileup --output indels.gtf
perl ../software/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input all_pops.pileup --gtf indels.gtf --output all_pops.idf.pileup

# filtering repeats
perl ../software/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input all_pops.idf.pileup --gtf ../ref/n_vespilloides.filteredRepeats.gff --output all_pops.idf.tef.pileup
perl ../software/popoolation2_1201/fst-sliding.pl --input all_pops.idf.tef.sync --output windows_500bp.maxDP600.mcnt6.idf.tef.fst --suppress-noninformative --min-count 6 --min-coverage 40 --max-coverage 2% --min-covered-fraction 0.5 --window-size 500 --step-size 250 --pool-size 82:104:104:118

#### Create sync file
java -jar ../software/popoolation2_1201/mpileup2sync.jar --fastq-type sanger --threads 16 --min-qual 20 --input all_pops.idf.tef.pileup --output all_pops.idf.tef.sync


#### Calculate genetic divergence between populations
# calculate FST for 500bp sliding windows
perl ../software/popoolation2_1201/fst-sliding.pl --input all_pops.idf.tef.sync --output windows_500bp.maxDP600.idf.tef.fst --suppress-noninformative --min-count 4 --min-coverage 40 --max-coverage 2% --min-covered-fraction 0.5 --window-size 500 --step-size 250 --pool-size 82:104:104:118

perl ../software/popoolation2_1201/snp-frequency-diff.pl --input all_pops.idf.tef.sync --output-prefix all_freq --min-count 6 --min-coverage 30 --max-coverage 700

# 500 bp windows FET
perl ../software/popoolation2_1201/fisher-test.pl --input all_pops.idf.tef.sync --output all_pops.idf.500bp.tef.fet --min-count 3 --min-coverage 30 --max-coverage 700 --suppress-noninformative --window-size 500 --step-size 250
