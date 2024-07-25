# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz |  |  |  |
| 1294_S1_L008_R2_001.fastq.gz |  |  |  |
| 1294_S1_L008_R3_001.fastq.gz |  |  |  |
| 1294_S1_L008_R4_001.fastq.gz |  |  |  |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem
    
    We need to de-multiplex fastq files in order to identify the level of index swapping and undetermined index-pairs.
   
3. Describe output

    We will output 48 FASTQ files of matching index pairs (read1 and read2 for 24 different index pairs), two FASTQ files with index-hopped reads-pairs, and two FASTQ files for undetermined (non-matching or low quality) index-pairs.
   
5. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
6. Pseudocode
7. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
