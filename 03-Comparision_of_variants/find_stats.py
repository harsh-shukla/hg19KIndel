import re
import sys


## usage   python find_stats.py <V1.vcf> <V2.vcf>
## example python find_stats.py hg19_varaiant_file hg19K_varaint_file   
 

f1  = open(sys.argv[1],"r")  # hg19
f2  = open(sys.argv[2],"r")  # hg19K


SNP_f1=0
Indel_f1=0

for line in f1:

    if line[0]=="#":
        continue

    col = re.split("\t|\n",line)

    ref = col[3]
    alt = col[4]

    if len(ref)==1 and len(alt)==1:
        SNP_f1 = SNP_f1 + 1

    else:
 
        if "," in alt:
            aa = re.split(",",alt)

            if (len(aa) == 2 and len(alt) == 3) or (len(aa)==3 and len(alt)==5):
                SNP_f1 = SNP_f1 + 1
    
            else:
                Indel_f1 = Indel_f1 + 1 

        else: 

            Indel_f1 = Indel_f1 + 1

f1.close()

SNP_f2=0
Indel_f2=0

for line in f2:

    if line[0]=="#":
        continue

    col = re.split("\t|\n",line)

    ref = col[3]
    alt = col[4]

    if len(ref)==1 and len(alt)==1:
        SNP_f2 = SNP_f2 + 1

    else:
 
        if "," in alt:
            aa = re.split(",",alt)

            if (len(aa) == 2 and len(alt) == 3) or (len(aa)==3 and len(alt)==5):
                SNP_f2 = SNP_f2 + 1 

            else:      
                Indel_f2 = Indel_f2 + 1                 

        else:    
            Indel_f2 = Indel_f2 + 1  

f2.close()

#################################################################################################################################


f1 = open("vcftools_out.diff.sites_in_files","r")
fout = open("mixed.txt","w")


for line in f1:
    fout.write(line)
    break

common=0
unique_hg19=0
unique_hg19K=0


SNP_same        = 0    # B
SNP_hetrozygous = 0    # O
SNP_unq_1 = 0
SNP_unq_2 = 0
SNP_non_match_O = 0

Indel_same = 0 #B
Indel_hetrozygous = 0 #O
Indel_unq_1 = 0
Indel_unq_2 = 0
Indel_non_match = 0

for line in f1:

    cols = re.split("\t|\n",line)

    present = cols[3]
    ref1    = cols[4]
    ref2    = cols[5]
    alt1    = cols[6]
    alt2    = cols[7]


    if (len(ref1)==1 and len(alt1)==1) and (len(ref2)==1 and len(alt2)==1):    # SNP only strictly

        if present == "B":
            SNP_same = SNP_same + 1

        elif present == "1":
            SNP_unq_1 = SNP_unq_1 + 1

        elif present == "2":
            SNP_unq_2 = SNP_unq_2 + 1 

        else:

            if ref1==alt2 and ref2==alt1:   # site is hetrozygous 
                SNP_hetrozygous = SNP_hetrozygous + 1 

            else:
                SNP_non_match_O = SNP_non_match_O + 1

    else:

        if (("," in alt1 and len(alt1)==3) or ("," in alt2 and len(alt2)==3)) and (len(ref1)==1 and len(ref2)==1):

            #print line,         
            if present == "B":
                SNP_same = SNP_same + 1

            elif present == "1":
                SNP_unq_1 = SNP_unq_1 + 1

            elif present == "2":
                SNP_unq_2 = SNP_unq_2 + 1 

            else:

                if ref1==alt2 and ref2==alt1:   # site is hetrozygous 
                    SNP_hetrozygous = SNP_hetrozygous + 1 
                else:
                    SNP_non_match_O = SNP_non_match_O + 1


        elif ( (len(ref1)>1 or len(alt1)>1) or (ref1=="." and alt1==".") ) and ( (len(ref2)>1 or len(alt2)>1) or (ref2=="." and alt2==".") ):        # Both variants are indels strictly  

            if present== "B":
                Indel_same = Indel_same + 1
  
            elif present == "1":
                Indel_unq_1 = Indel_unq_1 + 1

            elif present == "2":
                Indel_unq_2 = Indel_unq_2 + 1

            else:  

                if ref1==alt2 and ref2==alt1:   # site is hetrozygous 
                    Indel_hetrozygous = Indel_hetrozygous + 1 
                    #print line,

                else:
                    Indel_non_match = Indel_non_match + 1                          # same varaints
                    #print line,

        else:								       # One is SNP other is Indel and are overlapping "O"

            fout.write(line) 

            if (len(ref1) == 1 and len(alt1) == 1 ):
                SNP_unq_1 = SNP_unq_1 + 1
                Indel_unq_2 = Indel_unq_2 + 1

            elif (len(ref2) == 1 and len(alt2) == 1 ):
                SNP_unq_2 = SNP_unq_2 + 1 
                Indel_unq_1 = Indel_unq_1 + 1
  
            else:
                pass
                print line,

print "\n ****    SNP analyses  ****** "

print "SNP Same               : ",SNP_same
print "SNP_hetrozygous        : ",SNP_hetrozygous
print "SNP Unique to hg19     : ",SNP_unq_1
print "SNP_Unique to hg19K    : ",SNP_unq_2
print "SNP not match overlapp : ",SNP_non_match_O

print "\n ****    Indel analyses  ****** "

print "Indel Same               : ",Indel_same
print "Indel_hetrozygous        : ",Indel_hetrozygous
print "Indel Unique to hg19     : ",Indel_unq_1
print "Indel_Unique to hg19K    : ",Indel_unq_2
print "Indel not match overlapp : ",Indel_non_match


total_SNPs  = SNP_same   + SNP_hetrozygous   + SNP_unq_1   + SNP_unq_2    + SNP_non_match_O
total_Indel = Indel_same + Indel_hetrozygous + Indel_unq_1 + Indel_unq_2  + Indel_non_match


common_SNP = SNP_same + SNP_hetrozygous + SNP_non_match_O
fp_SNP     = (float(SNP_unq_1) / SNP_f1) * 100
fn_SNP     = (float(SNP_unq_2) / SNP_f2) * 100


common_Indel = Indel_same + Indel_hetrozygous + Indel_non_match
fp_Indel     = (float(Indel_unq_1) / Indel_f1) * 100
fn_Indel     = (float(Indel_unq_2) / Indel_f2) * 100



print "\n\n\n\n  ***** FINAL SNP ANALYSES  ********"

print " Common  SNP            : ",common_SNP
print "SNP Unique to hg19      : ",SNP_unq_1
print "SNP Unique to hg19K     : ",SNP_unq_2
print "False Positives         : ",fp_SNP
print "False Negatives         : ",fn_SNP 


print "\n\n  ***** FINAL Indel ANALYSES  ********"

print " Common Indels            : ",common_Indel
print "Indel Unique to hg19      : ",Indel_unq_1
print "Indel Unique to hg19K     : ",Indel_unq_2
print "False Positives           : ",fp_Indel
print "False Negatives           : ",fn_Indel 


print "\nStats for hg19"
print "Number of SNPs    :",SNP_f1
print "Number of Indels  :",Indel_f1


print"\n\nStats for hg19K"
print "Number of SNPs    :",SNP_f2
print "Number of Indels  :",Indel_f2


fout.close()








