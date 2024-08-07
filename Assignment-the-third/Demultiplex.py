#!/usr/bin/env python

import bioinfo

import argparse

# arguments to pass file names
def get_args():
    parser = argparse.ArgumentParser(description="A program used to demultiplex fastq files. Creates a forward and reverse file for each dual-matched index. Also creates a 1 forward and 1 reverse file for index-hopping reads and unknown index reads.")
    parser.add_argument("-r1", "--ForwardRead", help="The name of the R1 fastq file.", required=True, type=str)
    parser.add_argument("-r4", "--ReverseRead", help="The name of the R4 fastq file.", required=True, type=str)
    parser.add_argument("-i1", "--index1", help="The name of the index 1 fastq file.", required=True, type=str)
    parser.add_argument("-i2", "--index2", help="The name of the index 2 fastq file.", required=True, type=str)
    parser.add_argument("-ids", "--indexes", help="The name of the file containing the list of indexes used", required=True, type=str)
    return parser.parse_args()
	
args = get_args()


# Stats to Report
MatchedIndexCount = 0
PercentMatched = 0.0

UnknownIndexCount = 0 
PercentUnknown = 0.0

LowQualCount = 0
PercentLowQual = 0.0

HoppedCount = 0
PercentHopped = 0.0

TotalRecords = 0


# Extract indexes from indexes.txt and store as keys and set the value to 0, we will use this as a counter.
index_dict = {}


with open(args.indexes, "r") as indexes:
    next(indexes)
    for line in indexes:
        line = line.strip()
        line = line.split("\t")
        index = line[4]
        index_dict[index] = 0
indexes.close()


# Dictionary that relates the sample number to the index used
index_sample = {}

# Creating a list of indexes for open/closing matched files
index_list = []

with open(args.indexes, "r") as indexes:
    next(indexes)
    for line in indexes:
        line = line.strip()
        line = line.split("\t")
        index = line[4]
        sample = line[0]
        index_sample[index] = sample
        index_list.append(index)
indexes.close()



# Define reverse complement function
def reverse_complement(seq):
    '''This function returns the reverse complement of DNA'''
    # complement
    seq_comp = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq_comp = seq_comp.upper()
         
    # Reverse strand
    reverse_comp = seq_comp[::-1]
    return reverse_comp


# Open files to write to 
Unknown_Forward = open("UnknownIndexes_Forward.fq", "a")
Unknown_Reverse = open("UnknownIndexes_Reverse.fq", "a")
Hopped_Forward = open("HoppedIndexes_Forward.fq", "a")
Hopped_Reverse = open("HoppedIndexes_Reverse.fq", "a")


# Loop to open matched files
matchedfile_dict = {}

for i in index_list:
    fh1 = open(f"{i}_{i}_Forward_Matched.fq", "a")
    fh2 = open(f"{i}_{i}_Reverse_Matched.fq", "a")
    matchedfile_dict[i] = (fh1, fh2)


# Read each record from all 4 files into memory
import gzip
with gzip.open(args.ForwardRead, "rt") as R1, gzip.open(args.index1, "rt") as R2, gzip.open(args.index2, "rt") as R3, gzip.open(args.ReverseRead, "rt") as R4:
    while True:
        headerR1 = R1.readline().strip()
        seqR1 = R1.readline().strip()
        spacer = R1.readline().strip()
        qualityR1 = R1.readline().strip()
        

        headerR2 = R2.readline().strip()
        seqR2 = R2.readline().strip()
        spacer = R2.readline().strip()
        qualityR2 = R2.readline().strip()

        headerR3 = R3.readline().strip()
        seqR3 = R3.readline().strip()
        seqR3_revcomp = reverse_complement(seqR3)
        spacer = R3.readline().strip()
        qualityR3 = R3.readline().strip()

        headerR4 = R4.readline().strip()
        seqR4 = R4.readline().strip()
        spacer = R4.readline().strip()
        qualityR4 = R4.readline().strip()

        TotalRecords += 1

        # set new headers that include the indexes
        new_headerR1 = headerR1 + " " + seqR2 + "-" + seqR3_revcomp # type: ignore
       
        new_headerR4 = headerR4 + " " + seqR2 + "-" + seqR3_revcomp # type: ignore
        

        # break while True loop once at EOF
        if headerR1 == "":
            TotalRecords -= 1
            break
    


        # Check for unknown indexes
        if seqR2 not in index_dict or seqR3_revcomp not in index_dict:
            UnknownIndexCount += 1
            # Write to unknown file
            Unknown_Forward.write(f"{new_headerR1}\n{seqR1}\n{spacer}\n{qualityR1}\n")
            Unknown_Reverse.write(f"{new_headerR4}\n{seqR4}\n{spacer}\n{qualityR4}\n")

        # Check for index average quality scores below threshold
        elif bioinfo.qual_score(qualityR2) < 20 or bioinfo.qual_score(qualityR3) < 20: # type: ignore
            LowQualCount += 1
            # Write to unknown file
            Unknown_Forward.write(f"{new_headerR1}\n{seqR1}\n{spacer}\n{qualityR1}\n")
            Unknown_Reverse.write(f"{new_headerR4}\n{seqR4}\n{spacer}\n{qualityR4}\n")


            # check if they match
        elif seqR2 == seqR3_revcomp:
            MatchedIndexCount += 1
            index_dict[seqR2] += 1
            # write to matched files with the index in file name
            forward_matched = matchedfile_dict[seqR2][0]
            forward_matched.write(f"{new_headerR1}\n{seqR1}\n{spacer}\n{qualityR1}\n")
                
            reverse_matched = matchedfile_dict[seqR2][1]
            reverse_matched.write(f"{new_headerR4}\n{seqR4}\n{spacer}\n{qualityR4}\n")
                
        else:
            # They are index hopped
            HoppedCount += 1
            # write to hopped file
            Hopped_Forward.write(f"{new_headerR1}\n{seqR1}\n{spacer}\n{qualityR1}\n")
            Hopped_Reverse.write(f"{new_headerR4}\n{seqR4}\n{spacer}\n{qualityR4}\n")


# Close files opened before loop 
Unknown_Forward.close()
Unknown_Reverse.close()
Hopped_Forward.close()
Hopped_Reverse.close()

for file in matchedfile_dict:
    matchedfile_dict[file][0].close()
    matchedfile_dict[file][1].close()

# dictionary with counts for each sample
sample_count = {}

for index in index_sample:
    sample_count[index_sample[index]] = index_dict[index]
print(f"Sample counts are: {sample_count}")



# convert counts to percentage
sample_percentage = {}

for sample in sample_count:
    sample_percentage[sample] = sample_count[sample] / TotalRecords * 100
print(f"Sample percentages are: {sample_percentage}")

print(f"Matched index counts are: {index_dict}")

print(f"Number of matched indexes: {MatchedIndexCount}") 
PercentMatched = MatchedIndexCount / TotalRecords * 100
print(f"Percent dual matched: {PercentMatched}")


print(f"Number of unknown indexes: {UnknownIndexCount}") 
PercentUnknown = UnknownIndexCount / TotalRecords * 100
print(f"Percent unknown: {PercentUnknown}")


print(f"Number of low quality indexes (mean < 20) that don't contain N: {LowQualCount}")
PercentLowQual = LowQualCount / TotalRecords * 100
print(f"Percent Low Quality: {PercentLowQual}")


print(f"Number of hopped indexes: {HoppedCount}")
PercentHopped = HoppedCount / TotalRecords * 100
print(f"Percent Hopped: {PercentHopped}")


print(f"Total number of records: {TotalRecords}")