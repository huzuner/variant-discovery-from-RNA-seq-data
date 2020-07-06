# variant-discovery-from-RNA-seq-data
Variant discovery from RNA-seq data using GATK

#Variant identification from RNA-Seq data#
#by huzuner#

1) ALIGNMENT - STAR TWO-PASS

Download human genome assembly 38 contained in GATK Bundle:

Homo_sapiens_assembly38.fasta.gz file and associated files (index etc.) at ftp://ftp.broadinstitute.org/bundle/hg38 

*p.s. Other annotation files are listed in Base Recalibration step

There are four steps of alignment using STAR:

*STAR uses genome index files that must be saved in unique directories. The human genome index was built from the FASTA file hg19.fa as follows:

genomeDir=/path/to/genome
mkdir $genomeDir
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles genome.fa\ --runThreadN <n>

*Alignments are executed as follows:

runDir=/path/to/1pass
mkdir $runDir
cd $runDir
STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN

*A new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

genomeDir=/path/to/hg19_2pass
mkdir $genomeDir
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles genome.fa
--sjdbFileChrStartEnd /path/to/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN

*The resulting index is then used to produce the final alignments as follows:

runDir=/path/to/2pass
mkdir $runDir
cd $runDir
STAR --genomeDir $genomeDir --readFilesIn mate1.fq mate2.fq --runThreadN 

2) PICARD - ADDING READ GROUP INFORMATION

java -jar picard.jar AddOrReplaceReadGroups I=$DATADIR/data.bam  O=data_rg_added_sorted.bam SO=coordinate RGID=flowcell1.lane1 RGLB=Samplelib RGPL=illumina RGPU=unit1 RGSM=SampleID

3) PICARD - MARK DUPLICATES

java -jar picard.jar MarkDuplicates I=$DATADIR/data_rg_added_sorted.bam O=data_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 

4) GATK3 - Split'N'Trim

java -jar -Xms128G -Xmx128G GenomeAnalysisTK.jar -T SplitNCigarReads -R genome.fa -I $data_reordered_after_mark_duplicates.bam -o data_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

5) GATK4 - Base Recalibration

Download the known sites data first:

*dbSNP 144: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_144.hg38.vcf.gz AND ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_144.hg38.vcf.gz.tbi
*Mills and 1000G indels: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz AND ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
*Known indels (in BETA): ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz AND ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

Then perform base recalibration.
Base recalibration consists of two steps:

gatk BaseRecalibrator -I $DATADIR/split.bam -R genome.fa --known-sites Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites dbsnp_146.hg38.vcf.gz -O recal_data.table

gatk ApplyBQSR -I $DATADIR/split.bam -R genome.fa --bqsr-recal-file recal_data.table -O recal.bam

6) VARIANT CALLING - GATK3 HAPLOTYPECALLER

java -jar -Xms128G -Xmx128G  GenomeAnalysisTK.jar -T HaplotypeCaller -I $DATADIR/data_split.bam -R genome.fasta -dontUseSoftClippedBases -stand_call_conf 20.0 -o data_raw.vcf --num_cpu_threads_per_data_thread 48

p.s. you may need to use --output-mode EMIT_ALL_ACTIVE_SITES to reveal all genotypes, may be required for comparison

7) VARIANT FILTRATION

*GATK3
java -jar Xmx128G GenomeAnalysisTK.jar -T Variant Filtration -V data.vcf -R genome.fasta -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o variants_filtered.vcf

*GATK4
gatk Variant Filtration -V data.vcf -R genome.fa -window 35 -cluster 3 --filter-expression "FS > 30.0 || QD < 2.0" --filter-name "FS, QD" -O variants_filtered.vcf --num_cpu_threads_per_data_thread 48

8) SNPEFF - VARIANT ANNOTATION

java -Xmx16g -jar snpEff.jar -c snpEff.config -v hg38 -stats annotated-hg38.html $DATADIR/variants.vcf > annotated-hg38.ann.vcf
