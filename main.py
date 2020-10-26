# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


from Bio import SeqIO
import array
import csv
import os
import pandas as pd

#directory to look into to find multiple FASTA files
directory = "100nM_FASTAs"
file_count = 0
max_seq_length = 1500


position_row = [0]*max_seq_length
ii = 0
while ii < max_seq_length:
    position_row[ii] = ii
    ii += 1

rows = [position_row, position_row]
all_data = pd.DataFrame(position_row)

for filename in os.listdir(directory):
#This section opens the file and iterates through the sequencing to count how many sequences there are total and it
#records this as the size and also records the length of the sequence (should be the same for every sequence).
    file_count += 1
    aa = SeqIO.parse(directory + "/" + filename, "fasta")
    size = 0
    seq_length = 0
    for record in aa:
        size +=1
        seq_length = len(record.seq)

    print(seq_length)

#Initiating a lot of variables here. These count how many of each base shows up at each position. The lists are all of
#the same length of the sequence.
    A_count = [0]*max_seq_length
    T_count = [0]*max_seq_length
    C_count = [0]*max_seq_length
    G_count = [0]*max_seq_length
    Blank_count = [0]*max_seq_length

    #Loop through every sequence in the FASTA file
    for record in SeqIO.parse(directory + "/" + filename, "fasta"):

        #Initiate the counter for the position at the start of the sequence
        pos = 0

        #Loop through each position of the sequence and add 1 to the value corresponding to that position in the correct
        #list
        while pos < seq_length:
            if record.seq[pos] == "a":
                A_count[pos] += 1
            if record.seq[pos] == "t":
                T_count[pos] += 1
            if record.seq[pos] == "c":
                C_count[pos] += 1
            if record.seq[pos] == "g":
                G_count[pos] += 1
            if record.seq[pos] == "-":
                Blank_count[pos] += 1
            pos += 1


    #Initiate some variables to convert all of the above values into percentages
    A_percent = []
    T_percent = []
    C_percent = []
    G_percent = []
    Blank_percent = []
#
    #Loop through each array and convert to percentages
    for number in A_count:
        A_percent.append(number/size*100)
    for number in T_count:
        T_percent.append(number/size*100)
    for number in C_count:
        C_percent.append(number/size*100)
    for number in G_count:
        G_percent.append(number/size*100)
    for number in Blank_count:
        Blank_percent.append(number/size*100)


    #Need a section where we will start with the percent lists and remove every index where the Blanks are there over
    #99% of the time and add an additional item to the end of the lists
    ii = 0
    while ii < 1390:
        if Blank_percent[ii] > 50:
            #remove that index from every list
            del A_percent[ii]
            del T_percent[ii]
            del C_percent[ii]
            del G_percent[ii]
            del Blank_percent[ii]

            #add 0s to the end of every list to keep them all the same size
            A_percent.append(0)
            T_percent.append(0)
            C_percent.append(0)
            G_percent.append(0)
            Blank_percent.append(0)

            ii -= 1
        ii += 1



    if file_count == 1:
        percents = [A_percent, T_percent, C_percent, G_percent, Blank_percent]
        all_data = pd.DataFrame(percents, index=["File 1 A", "File 1 T", "File 1 C", "File 1 G", "File 1 Blank"])
    else:
        all_data.loc["File " + str(file_count) + " A"] = A_percent
        all_data.loc["File " + str(file_count) + " T"] = T_percent
        all_data.loc["File " + str(file_count) + " C"] = C_percent
        all_data.loc["File " + str(file_count) + " G"] = G_percent
        all_data.loc["File " + str(file_count) + " Blank"] = Blank_percent

    print(file_count)


#outside the loop, we make the fakeCSV file with all of the different percentages given
all_data.to_csv("100nM_Output_v2.csv")
