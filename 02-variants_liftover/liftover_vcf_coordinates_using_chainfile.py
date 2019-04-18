import re
import sys
from pyliftover import LiftOver

###############################################################################################################################  

# Usage  :  python liftover_vcf_coordinates_using_chainfile.py hg19Kindeltohg19.over.chain <INPUT_VCF_FILE.vcf>

# Input  :  1) <INPUT_VCF_FILE.vcf>              - variants called using hg19Kindel as reference
#           2) hg19Kindeltohg19.over.chain       - chain file used for lifting cordinates from hg19Kindel coordinate space to hg19 coordinate space 

# Output :  1) INPUT_VCF_FILE.vcf.CONVERTED.vcf  - all the variants with positions lifted to hg19 coordinate frame
#           2) INPUT_VCF_FILE.vcf.UNMAPPED.vcf   - the variants that could not be mapped. (That particualr position is deleted in hg19)

# Dependency : PyLiftover Library; Ref:https://pypi.org/project/pyliftover/
             
###############################################################################################################################  


## INPUT  ##

chain_file = sys.argv[1]               # chain file for converting
vcf_file   = open(sys.argv[2],"r")     # vcf file - variants which have to be lifted Over


##  Output  ###

coll = re.split("/",sys.argv[2])
fname = coll[-1]

filename1 = fname + ".CONVERTED.vcf"
fout = open(filename1,"w")
filename2 = fname + ".UNMAPPED.vcf"
funmapped = open(filename2,"w")


lo = LiftOver(chain_file)

for line in vcf_file:

    if line[0:3]!="chr":
        fout.write(line)
        continue

    cols = re.split("\t|\n",line)
    chr_name = cols[0]
    pos = int(cols[1])

    pos_0_based = pos - 1          # Since vcf is 1-based and the conversion works on 0-based coordinate                              
    mapped_list = lo.convert_coordinate(chr_name,pos_0_based,'+')

    ## Sanity check : convert cordinate in theory should return only one mapped coordinate (for this particular scenario)   
    if len(mapped_list) > 1:
        print " More than one mapping detected at ",chr_name," : ",str(pos) 
        continue

    if len(mapped_list) == 0:   # The position cannot be mapped to hg19 (that particular position is deleted in hg19)
        funmapped.write(line)
        continue

    mapped_chr = mapped_list[0][0]

    ## Sanity check : if the position maps to a different chromosome and not the same chromosome
    if mapped_chr != chr_name:
        print " Different chromosome mapping detected at",chr_name," : ",str(pos)
        sys.exit(0)

    mapped_pos = mapped_list[0][1] + 1 # Converting 0-based to 1-based coordinates 

    write_line = chr_name + "\t" + str(mapped_pos) + "\t" + cols[2] + "\t" + cols[3] + "\t" + cols[4] + "\t" + cols[5] + "\t" + cols[6] + "\t" + cols[7] + "\t" + cols[8] + "\t" + cols[9] + "\n"
    fout.write(write_line)  
 



