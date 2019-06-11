#USAGE
#
# ./1_FASTQalign.sh yourinput.fastq bowtie_built_reference 
#
#    INPUT:  program will use the bowtie_built_reference file and carry out alignment of all input sequences
#    OUTPUT: program will output a new file name youinput.fastq.map, consisting of mapped reads
#
#EXAMPLE USE
#
#    for i in `ls ../FASTQout`; do ./1_FASTQalign.sh ../FASTQout/$i URA3; done
#
#       program takes each fastq file in `ls ../FASTQout` and 
#       carries out alignment on the reference sequence URA3 (note this is built using bowtie build)

#print name of file being processed
echo $1

#bowtie mapping
bowtie $2 -5 4 -3 4 -n 3 -e 450 -l 16 $1 > $1.map 
