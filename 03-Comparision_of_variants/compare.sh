# Sort
vcf-sort $1 > $1.sort.vcf
vcf-sort $2 > $2.sort.vcf

vcftools --vcf $1.sort.vcf --out vcftools_out --diff $2.sort.vcf --diff-site

rm $1.sort.vcf
rm $2.sort.vcf

python find_stats.py $1 $2 > final_out.txt

