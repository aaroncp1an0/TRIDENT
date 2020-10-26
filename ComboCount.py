# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


from Bio import SeqIO
import array
import csv
import os
import pandas as pd

#directory to look into to find multiple FASTA files
directory = "2nM_FASTAs"
file_count = 0
max_seq_length = 1500
seq_count = 0


position_row = [0]*max_seq_length
ii = 0
while ii < max_seq_length:
    position_row[ii] = ii
    ii += 1

rows = [position_row, position_row]
all_data = pd.DataFrame(position_row)

offset = 7
offsets = []

D54 = 0
Y57 = 0
Q237 = 0
N277 = 0


D54N = 0
Y57H = 0
Q237S = 0
N277D = 0
D54N_Y57H = 0
D54N_Q237S = 0
D54N_N277D = 0

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



    #Loop through every sequence in the FASTA file
    for record in SeqIO.parse(directory + "/" + filename, "fasta"):

        #Initiate the counter for the position at the start of the sequence
        seq = ""
        pos = 0
        seq_count += 1

        while pos < seq_length:
            if record.seq[pos] != "-":
                seq += record.seq[pos]
            pos +=1

        pos = 0


        while pos < len(seq):
            if seq[pos: pos + 3] == "atg":
                offset = pos
                break
            pos += 1

        if offset > 20:
            offset = 7

        if seq[160 + offset:160 + offset + 4] == "aaca" or seq[160 + offset + 1:160 + offset + 5] == "aaca" or seq[160 + offset - 1:160 + offset + 3] == "aaca" or seq[160 + offset + 2:160 + offset + 6] == "aaca" or seq[160 + offset - 2:160 + offset + 2] == "aaca":
            D54N += 1
            if seq[169 + offset:169 + offset + 4] == "cact" or seq[169 + offset + 1:169 + offset + 5] == "cact" or seq[169 + offset + 2:169 + offset + 6] == "cact" or seq[169 + offset - 1:169 + offset + 3] == "cact" or seq[169 + offset - 2:160 + offset + 2] == "cact":
                D54N_Y57H += 1
            if seq[709 + offset:709 + offset + 4] == "taaa" or seq[709 + offset + 1:709 + offset + 5] == "taaa" or seq[709 + offset + 2:709 + offset + 6] == "taaa" or seq[709 + offset - 1:709 + offset + 3] == "taaa" or seq[709 + offset - 2:709 + offset + 2] == "taaa":
                D54N_Q237S += 1
            if seq[829 + offset:829 + offset + 4] == "gact" or seq[829 + offset + 1:829 + offset + 5] == "gact" or seq[829 + offset + 2:829 + offset + 6] == "gact" or seq[829 + offset - 1:829 + offset + 3] == "gact" or seq[829 + offset - 2:829 + offset + 2] == "gact":
                D54N_N277D += 1
        if seq[160 + offset:160 + offset + 4] == "gaca" or seq[160 + offset + 1:160 + offset + 5] == "gaca" or seq[160 + offset - 1:160 + offset + 3] == "gaca" or seq[160 + offset + 2:160 + offset + 6] == "gaca" or seq[160 + offset - 2:160 + offset + 2] == "gaca":
            D54 += 1
        if seq[169 + offset:169 + offset + 4] == "cact" or seq[169 + offset + 1:169 + offset + 5] == "cact" or seq[169 + offset + 2:169 + offset + 6] == "cact" or seq[169 + offset - 1:169 + offset + 3] == "cact" or seq[169 + offset - 2:169 + offset + 2] == "cact":
            Y57H +=1
        if seq[169 + offset:169 + offset + 4] == "tact" or seq[169 + offset + 1:169 + offset + 5] == "tact" or seq[169 + offset + 2:169 + offset + 6] == "tact" or seq[169 + offset - 1:169 + offset + 3] == "tact" or seq[169 + offset - 2:169 + offset + 2] == "tact":
            Y57 +=1
        if seq[709 + offset:709 + offset + 4] == "taaa" or seq[709 + offset + 1:709 + offset + 5] == "taaa" or seq[709 + offset + 2:709 + offset + 6] == "taaa" or seq[709 + offset - 1:709 + offset + 3] == "taaa" or seq[709 + offset - 2:709 + offset + 2] == "taaa":
            Q237S += 1
        if seq[709 + offset:709 + offset + 4] == "caaa" or seq[709 + offset + 1:709 + offset + 5] == "caaa" or seq[709 + offset + 2:709 + offset + 6] == "caaa" or seq[709 + offset - 1:709 + offset + 3] == "caaa" or seq[709 + offset - 2:709 + offset + 2] == "caaa":
            Q237 += 1
        if seq[829 + offset:829 + offset + 4] == "gact" or seq[829 + offset + 1:829 + offset + 5] == "gact" or seq[829 + offset + 2:829 + offset + 6] == "gact" or seq[829 + offset - 1:829 + offset + 3] == "gact" or seq[829 + offset - 2:829 + offset + 2] == "gact":
            N277D += 1
        if seq[829 + offset:829 + offset + 4] == "aact" or seq[829 + offset + 1:829 + offset + 5] == "aact" or seq[829 + offset + 2:829 + offset + 6] == "aact" or seq[829 + offset - 1:829 + offset + 3] == "aact" or seq[829 + offset - 2:829 + offset + 2] == "aact":
            N277 += 1
    print(file_count)

print("D54N Count: " + str(D54N))
print("D54 Count: " + str(D54))
print("Total Sequence Count: " + str(seq_count))
print("Percent D54 Accounted for: " + str(round(100*(D54 + D54N)/seq_count,2)) + "%")

print("Y57H Count: " + str(Y57H))
print("Y57 Count: " + str(Y57))
print("Total Sequence Count: " + str(seq_count))
print("Percent Y57 Accounted for: " + str(round(100*(Y57 + Y57H)/seq_count,2)) + "%")

print("Q237* Count: " + str(Q237S))
print("Q237 Count: " + str(Q237))
print("Total Sequence Count: " + str(seq_count))
print("Percent Q237 Accounted for: " + str(round(100*(Q237 + Q237S)/seq_count,2)) + "%")

print("N277D Count: " + str(N277D))
print("N277 Count: " + str(N277))
print("Total Sequence Count: " + str(seq_count))
print("Percent Q237 Accounted for: " + str(round(100*(N277 + N277D)/seq_count,2)) + "%")

print("D54N + Y57H: " + str(D54N_Y57H))
print("D54N + Q237*: " + str(D54N_Q237S))
print("D54N + N277D: " + str(D54N_N277D))