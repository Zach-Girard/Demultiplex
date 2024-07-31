#!/usr/bin/env python

import bioinfo

import argparse

# arguments to pass file names
def get_args():
    parser = argparse.ArgumentParser(description="A program used to calculate and plot the mean quality score per base position.")
    parser.add_argument("-r1", "--ForwardRead", help="The name of the R1 fastq file.", required=True, type=str)
    parser.add_argument("-r4", "--ReverseRead", help="The name of the R4 fastq file.", required=True, type=str)
    parser.add_argument("-i1", "--index1", help="The name of the index 1 fastq file.", required=True, type=str)
    parser.add_argument("-i2", "--index2", help="The name of the index 2 fastq file.", required=True, type=str)
    return parser.parse_args()
	
args = get_args()






Read1 = args.ForwardRead


# Initialize function for sequence read quality scores
import copy
def init_list101(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    for i in range(101):
        new_value = copy.copy(value)
        lst.append(new_value)
    return lst   


# Initialize list for forward reads
Read1_list = []
init_list101(Read1_list)


# Read forward read file and sum the quality scores for each base position
import gzip
counter = 0
with gzip.open(Read1, 'rt') as fh:
    for line in fh:
        line = line.strip()
        counter += 1
        if counter%4 == 0:
            index = 0
            for phred_score in line:
                Read1_list[index] += bioinfo.convert_phred(phred_score) # type: ignore
                index += 1 



# Take the sum in each position and calculate the mean
index = 0
for sum in Read1_list:
    Read1_list[index] = sum/ (counter / 4)
    index += 1






Read4 = args.ReverseRead

# Initialize list for reverse reads
Read4_list = []
init_list101(Read4_list)

# Read reverse read file and sum the quality scores for each base position
counter = 0
with gzip.open(Read4, 'rt') as fh:
    for line in fh:
        line = line.strip()
        counter += 1
        if counter%4 == 0:
            index = 0
            for phred_score in line:
                Read4_list[index] += bioinfo.convert_phred(phred_score) # type: ignore
                index += 1 

# Take the sum in each position and calculate the mean
index = 0
for sum in Read4_list:
    Read4_list[index] = sum/ (counter / 4)
    index += 1


# Initialize function for index read quality scores
def init_list8(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 8 values of 0.0.'''
    for i in range(8):
        new_value = copy.copy(value)
        lst.append(new_value)
    return lst   







Index1 = args.index1

# Initialize list for first index
Index1_list = []
init_list8(Index1_list)

# Read first index file and sum the quality scores for each base position
counter = 0
with gzip.open(Index1, 'rt') as fh:
    for line in fh:
        line = line.strip()
        counter += 1
        if counter%4 == 0:
            index = 0
            for phred_score in line:
                Index1_list[index] += bioinfo.convert_phred(phred_score) # type: ignore
                index += 1 

# Take the sum in each position and calculate the mean
index = 0
for sum in Index1_list:
    Index1_list[index] = sum/ (counter / 4)
    index += 1






Index2 = args.index2

# Initialize list for second index
Index2_list = []
init_list8(Index2_list)

# Read first index file and sum the quality scores for each base position
counter = 0
with gzip.open(Index2, 'rt') as fh:
    for line in fh:
        line = line.strip()
        counter += 1
        if counter%4 == 0:
            index = 0
            for phred_score in line:
                Index2_list[index] += bioinfo.convert_phred(phred_score) # type: ignore
                index += 1 

# Take the sum in each position and calculate the mean
index = 0
for sum in Index2_list:
    Index2_list[index] = sum/ (counter / 4)
    index += 1





# Histograms for mean quality scores per base position

import matplotlib.pyplot as plt # type: ignore

# Plot for the forward read
plt.plot(Read1_list, color = "red")
plt.title("Mean Quality Score per Base - Forward Read")
plt.xlabel("Base Position")
plt.ylabel("Mean Quality Score")
plt.ylim(0,41)
plt.savefig("ForwardRead.png")
plt.close()


# Plot for the reverse read
plt.plot(Read4_list, color = "red")
plt.title("Mean Quality Score per Base - Reverse Read")
plt.xlabel("Base Position")
plt.ylabel("Mean Quality Score")
plt.ylim(0,41)
plt.savefig("ReverseRead.png")
plt.close()


# Plot for the first index 
plt.plot(Index1_list, color = "red")
plt.title("Mean Quality Score per Base - Index 1")
plt.xlabel("Base Position")
plt.ylabel("Mean Quality Score")
plt.ylim(0,41)
plt.savefig("Index1.png")
plt.close()


# Plot for the second index
plt.plot(Index2_list, color = "red")
plt.title("Mean Quality Score per Base - Index 2")
plt.xlabel("Base Position")
plt.ylabel("Mean Quality Score")
plt.ylim(0,41)
plt.savefig("Index2.png")
plt.close()

