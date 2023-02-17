# Pooled Sequencing data analysis

# trim reads
./software/TrimGalore/trim_galore --fastqc --illumina --paired -o ./trimmed_reads raw_data/F1_H7HC7DSX2_L2_1.fq raw_data/F1_H7HC7DSX2_L2_2.fq
./software/TrimGalore/trim_galore --fastqc --illumina --paired -o ./trimmed_reads raw_data/F2_H7HC7DSX2_L2_1.fq raw_data/F2_H7HC7DSX2_L2_2.fq
./software/TrimGalore/trim_galore --fastqc --illumina --paired -o ./trimmed_reads raw_data/N1_H7HC7DSX2_L2_1.fq raw_data/N1_H7HC7DSX2_L2_2.fq
./software/TrimGalore/trim_galore --fastqc --illumina --paired -o ./trimmed_reads raw_data/N2_H7HC7DSX2_L2_1.fq raw_data/N2_H7HC7DSX2_L2_2.fq

# align to vespilloides genome
cd trimmed_reads

for i in *_H7HC7DSX2_L2_1_val_1.fq.gz;
do
filename=`echo $i | awk -F '_H7HC7DSX2_L2_1_val_1.fq.gz' '{print $1}'`
# align
bwa mem -t 8 ../ref/nves.fa $filename"_H7HC7DSX2_L2_1_val_1.fq.gz" $filename"_H7HC7DSX2_L2_2_val_2.fq.gz" > ../map/$filename"_aligned.sam"
done

cd ../map

# for all mapped files

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

# create pileup
samtools mpileup -B F1.qf.rmd.sort.bam F2.qf.rmd.sort.bam N1.qf.rmd.sort.bam N2.qf.rmd.sort.bam > all_pops.pileup

# filtering indels
perl ../software/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --indel-window 5 --min-count 2 --input all_pops.mpileup --output indels.gtf
perl ../software/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input all_pops.pileup --gtf indels.gtf --output all_pops.idf.pileup

# create sync file
#java -jar ../software/popoolation2_1201/mpileup2sync.jar --fastq-type sanger --threads 16 --min-qual 15 --input all_pops.idf.norep.pileup --output all_pops.idf.sync

# calculate FST for 500bp sliding windows
perl ../software/popoolation2_1201/fst-sliding.pl --input all_pops.idf.sync --output maxdp400.windows_500bp.idf.fst --suppress-noninformative --min-count 4 --min-coverage 40 --max-coverage 400 --min-covered-fraction 1 --window-size 500 --step-size 250 --pool-size 82:104:104:118
