
compare.sh -- shell script that takes two vcfs as input and outputs comparision statistics

In the first step,
It first sort and compares the two vcf using vcftools --diff-site option.

In second step <find_stats.py> python script parses the output from the step 1 to generate comaparison statistics

Eg : ./compare.sh <S1_hg19.vcf> <S1_hg19Kindel.vcf>

Here S1_hg19.vcf       is the vcf called for sample 1 using hg19       as reference
     S1_hg19Kindel.vcf is the vcf called for sample 1 using hg19Kindel as reference 
