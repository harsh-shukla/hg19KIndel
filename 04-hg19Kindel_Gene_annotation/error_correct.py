import re
import sys


## Input ##


f_refSeq         = open(sys.argv[1],"r")

fout = open("hg19Kindel_refGene.txt","w")

include_chr=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]

for line in f_refSeq: 
    fout.write(line)
    break

for line in f_refSeq:

    cols = re.split("\t|\n",line)

    chrm_refSeq   = cols[2]

    if chrm_refSeq not in include_chr:
        fout.write(line) 
        continue

    bin_rec       = cols[0]
    name          = cols[1]
#   chrm_refSeq   = cols[2]
    strand        = cols[3]
 
    txStart       = int(cols[4])    
    txEnd         = int(cols[5])
    cdsStart      = int(cols[6])
    cdsEnd        = int(cols[7])
    exonCount     = int(cols[8])

    exonStarts    = cols[9]
    exonEnds      = cols[10]
    score         = cols[11]
    name2         = cols[12]
    cdsStartStat  = cols[13]
    cdsEndStat    = cols[14]
    exonFrames    = cols[15]

    #print cols

    cols_exonStarts = re.split(",",exonStarts)
    cols_exonEnds   = re.split(",",exonEnds)
    cols_exonFrames = re.split(",",exonFrames)


    #print cols_exonFrames

    flag=0

    start  = []
    end    = []
    frames = []


    

    for j in xrange(0,exonCount):

        if j==0:

            current_starts = int(cols_exonStarts[j])
            current_ends   = int(cols_exonEnds[j])
            current_frame  = int(cols_exonFrames[j])
  
            start.append(current_starts)
            end.append(current_ends)
            frames.append(current_frame)  

            previous_starts = current_starts
            previous_ends   = current_ends

        else:

            current_starts = int(cols_exonStarts[j])
            current_ends   = int(cols_exonEnds[j])
            current_frame  = int(cols_exonFrames[j])
        
            if current_starts == previous_ends:
                end[-1]   = current_ends 
                flag=1

                if strand == "-":
                    frames[-1] = current_frame
  
                previous_starts =  previous_starts 
                previous_ends   =  current_ends
 
            else: 
                start.append(current_starts)
                end.append(current_ends)
                frames.append(current_frame)  

                previous_starts = current_starts
                previous_ends   = current_ends



    if flag==1:

        exonStarts_converted = ""
        exonEnds_converted   = ""
        exonFrames_converted = "" 

        if len(start) != len(end):
            print "Erorr ",line    

        exonCount_converted = len(start)

        for i in xrange(0,len(start)):

            exonStarts_converted = exonStarts_converted + str(start[i]) + ","
            exonEnds_converted   = exonEnds_converted + str(end[i]) + ","
            exonFrames_converted = exonFrames_converted + str(frames[i]) + ","


        write_line = bin_rec + "\t" + name + "\t" + chrm_refSeq + "\t" + strand + "\t" + str(txStart) + "\t" + str(txEnd) + "\t" + str(cdsStart) + "\t" + str(cdsEnd) + "\t" + str(exonCount_converted) + "\t" + exonStarts_converted + "\t" + exonEnds_converted + "\t" + score + "\t" + name2 + "\t" +  cdsStartStat + "\t" + cdsEndStat + "\t" + exonFrames_converted + "\n"


        fout.write(write_line)        

        print name

    else:

        fout.write(line)




 




