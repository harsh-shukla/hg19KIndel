############################################################## USAGE ###########################################################
#
#  Run   : python replace_variants.py <INPUT_1.fasta> <INPUT_2.vcf> <OUTPUT_1.fasta>
#   
#  INPUT :  
#  <INPUT_1.fasta> is the fasta file in which you want to replace the indels
#  <INPUT_2.vcf>   is the vcf file containing the variants to be replaced

#  OUTPUT :
#  <OUTPUT_1.fasta> - The output fasta file with the replaced variants 
#
#  DESCRIPTION :
#  This python script takes a fasta file and a vcf file as input and replaces for all the records (variants) the REF with the ALT field to give a modified/updated fasta file. 
##################################################################################################################################


def Make_Substitutions_Current_CHROM(current_chr_local,current_chr_name_l):

    len_curr_chr = len(current_chr_local)
    #print len_curr_chr

    # Here we are working with a local copy of entire chromosome string ; manipulating and returning it; cannot pass an immutable object through pass by reference ; also using the globla keyword was not viable 
    

    print "\nProcessing ",current_chr_name_l," for variants................."

    vcf_file.seek(0, 0)               # bring the file pointer to the starting of the vcf file

    displaced=0                       # keep track of iteratively changing REF POS value ; since you make some changes the REF POS will shift

    for line in vcf_file:

        if line[0:3]!="chr":          # ignore the header lines 
            continue

        col = re.split("\t|\n",line)

        current_line_chrom = col[0]

        if current_line_chrom != current_chr_name_l:    # variants in rest of chromosomes is ignored
            continue

        POS = int(col[1])
        REF = col[3]
        ALT = col[4]

        len_REF = len(REF)
        len_ALT = len(ALT) 

        NEW_POS = POS + displaced

        ref_start_index = NEW_POS - 1                     # since a string in python starts from 0 
        ref_end_index   = ref_start_index +  len_REF
    

        #print POS,REF,ALT,current_chr_name_l
        #print current_chr_local[ref_start_index:ref_end_index]       

        current_chr_local = current_chr_local[:ref_start_index] + ALT + current_chr_local[ref_end_index:]

        #print current_chr_local[(ref_start_index-10):(ref_start_index+len_ALT+10)]

        displaced = displaced - len_REF + len_ALT

        #break 

    print "    : Done"
        
    return current_chr_local


def Write_Modified_CHROM():

    

    len_curr_chr = len(current_chr)
    #print len_curr_chr
     
    print "Writing modified",current_chr_name,"to new fasta file ................."

    track_print=0
    
    fasta_header = ">" + current_chr_name
    fout.write(fasta_header)
    fout.write("\n")

    
    while(track_print < len_curr_chr):
  
        fout.write(current_chr[track_print:(track_print + BASES_PER_LINE)])
        fout.write("\n")

        track_print = track_print + BASES_PER_LINE  
  

    fout.write(current_chr[track_print:])


    if last==False:
        fout.write("\n") 
        fout.write("\n")    

    print "    : Done\n\n"
    


import re
import sys
import os


## Process the input files #

fasta_file = open(sys.argv[1],"r")          # The fasta file in which replaement should happen #
vcf_file   = open(sys.argv[2],"r")          # The vcf file that contains indel for 


## Sub section to create the output file

#file_ff = re.split("\.",sys.argv[1])
#file_vf = re.split("\.",sys.argv[2])
#fout_filename = file_ff[0] + "_" + file_vf[0] + ".fa" 

fout_filename = sys.argv[3]
fout       = open(fout_filename,"w")


## Paramters to control the structure of written output fasta file

BASES_PER_LINE = 60   # DEFAULT IS 60


## Find out all chromosome present in the vcf file

chr_present=[]

for line in vcf_file:

    if line[0:3]!="chr":           
        continue

    col1 = re.split("\t|\n",line)

    chromosome = col1[0]

    if chromosome not in chr_present:

        chr_present.append(chromosome)


#print chr_present 

## Process one chromosome at a time

current_chr = ""
first=True
last=False

for line in fasta_file:


    if line[0]==">":            # starting of a new chromosome


        if first==False:         # process the previous chromosome for substitutions

            if current_chr_name in chr_present:
                current_chr = Make_Substitutions_Current_CHROM(current_chr,current_chr_name)

            Write_Modified_CHROM()

        del current_chr          
 
        current_chr = ""
      
        col1 = re.split(">|\n",line)
        current_chr_name = col1[1]

        first=False
        continue 

    col1 = re.split("\n",line)

    current_chr = current_chr + col1[0]

    #print current_chr

   
#Process the last fasta sequence :  chromosome

#len_curr_chr = len(current_chr)
#print len_curr_chr

if current_chr_name in chr_present:
    current_chr = Make_Substitutions_Current_CHROM(current_chr,current_chr_name)

last=True
Write_Modified_CHROM()



# Clean up
del current_chr


# Closing all the files

 
fasta_file.close()
vcf_file.close()
fout.close()
