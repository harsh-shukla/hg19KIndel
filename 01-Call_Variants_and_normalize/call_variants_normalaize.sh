date

####
echo Step1 Alignment Using Bowtie and Conversions

bowtie2 -p 10 -x /home/cancer/Harsh/hg19K_and_hg19Kindel/hg19/hg19 -1 ../EBI_FTP_FASTQ/ERR194147_1.fastq.gz -2 ../EBI_FTP_FASTQ/ERR194147_2.fas
tq.gz -S ERR194147.sam 1>>ERR194147.stdout 2>>ERR194147.stderr

samtools view -bS ERR194147.sam -o ERR194147.bam
samtools sort -@ 4 -m 2G ERR194147.bam ERR194147_sorted
samtools index ERR194147_sorted.bam

rm ERR194147.bam
rm ERR194147.sam

####
echo Step2 remove duplicates using picard

~/Harsh/Software/jdk1.8.0_151/bin/java -Xmx12g -jar ~/Harsh/Software/picard/picard.jar MarkDuplicates I=ERR194147_sorted.bam O=ERR194147_sorted.dupRem.bam M=ERR194147_sorted.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true


####
echo Step3 running mpileup to generate binary vcf bcf files

samtools mpileup -ugf /home/cancer/Harsh/hg19K_and_hg19Kindel/hg19/hg19.fa  ERR194147_sorted.dupRem.bam > ERR194147_sorted.dupRem.bam.mpileup.bcf


####
echo Step4 converting and extrating only variants from vcf files

bcftools call -vmO z -o ERR194147_sorted.dupRem.bam.mpileup.bcf.call.vcf.gz ERR194147_sorted.dupRem.bam.mpileup.bcf

####
echo Step5 indexing the vcf files

tabix -p vcf ERR194147_sorted.dupRem.bam.mpileup.bcf.call.vcf.gz

####
echo Step6 generating stats of variants from vcf files

bcftools stats -F /home/cancer/Harsh/hg19K_and_hg19Kindel/hg19/hg19.fa -s - ERR194147_sorted.dupRem.bam.mpileup.bcf.call.vcf.gz > ERR194147_sorted.dupRem.bam.mpileup.bcf.call.vcf.gz.stats

####
echo Step7 Filtering 

bcftools filter -i'%QUAL>10 & DP>3' ERR194147_sorted.dupRem.bam.mpileup.bcf.call.vcf.gz > ERR194147_sorted.dupRem.bam.mpileup.bcf.call.qualGT10.dpGT3.vcf

rm ERR194147_sorted.dupRem.bam.mpileup.bcf


### Normalizing the variants

# STEP 1
vcfallelicprimitives ERR194147_sorted.dupRem.bam.mpileup.bcf.call.qualGT10.dpGT3.vcf > step1.vcf

# STEP 2
vt normalize step1.vcf -r /home/cancer/Harsh/hg19K_and_hg19Kindel/hg19/hg19.fa | vt uniq - -o ERR194147_sorted.dupRem.bam.mpileup.bcf.call.qualGT10.dpGT3.normalized.uniq.vcf

rm step1.vcf

date


