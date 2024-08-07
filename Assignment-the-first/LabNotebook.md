# Demulitplexing Lab Notebook


## Path to fastq files
```/projects/bgmp/shared/2017_sequencing```

# First Assignment Part 1 Questions
Length of read in 1294_S1_L008_R1_001.fastq.gz:

```zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | awk 'NR == 2 || NR % 4 == 2' | head -50 | awk '{ print length }'```

This gives me a charcter length of 101 in each line.


Length of read in 1294_S1_L008_R2_001.fastq.gz:

```zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | awk 'NR == 2 || NR % 4 == 2' | head -50 | awk '{ print length }'```

This gives me a charcter length of 8 in each line.


Length of read in 1294_S1_L008_R3_001.fastq.gz:


```zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | awk 'NR == 2 || NR % 4 == 2' | head -50 | awk '{ print length }'```

This gives me a charcter length of 8 in each line.


Length of read in 1294_S1_L008_R4_001.fastq.gz:

```zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | awk 'NR == 2 || NR % 4 == 2' | head -50 | awk '{ print length }'```

This gives me a charcter length of 101 in each line.


Phred encoding:

```zcat 1294_S1_L008_R1_001.fastq.gz | head -4 | tail -1```
```A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ```

The ```<``` indicates a +33 phred encoding.


```zcat 1294_S1_L008_R2_001.fastq.gz | head -4 | tail -1```
```#AA<FJJJ```

The ```<``` indicates a +33 phred encoding.


```zcat 1294_S1_L008_R3_001.fastq.gz | head -4 | tail -1```
```#AAAAJJF```

The ```#``` indicates a +33 phred encoding.


```zcat 1294_S1_L008_R4_001.fastq.gz | head -4 | tail -1```
```#AAFAFJJ-----F---7-<FA-F<AFFA-JJJ77<FJFJFJJJJJJJJJJAFJFFAJJJJJJJJFJF7-AFFJJ7F7JFJJFJ7FFF--A<A7<-A-7--```

The ```#``` and ```-``` indicate a +33 phred encoding.



## Table

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | forward read | 101 | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index 1 | 8 | +33 |
| 1294_S1_L008_R3_001.fastq.gz | index 2 | 8 | +33 |
| 1294_S1_L008_R4_001.fastq.gz | reverse read | 101 | +33 |


## Questions Continued...
1. What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.

I would suggest a quality score cutoff of 20 for the average qscore of the index reads. We are first filtering out all index reads with N's, so highly unconfident base calls will already be removed. The likelihood of including incorrect dual-matched index reads to your output files is very low. I would not filter the biological reads by quality score, because it would not be included after alignment. 


2. How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command(s) you used. CHALLENGE: use a one-line command)


```zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | awk 'NR == 2 || NR % 4 == 2' | grep "N" | wc -l```
Count: 3976613

```zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | awk 'NR == 2 || NR % 4 == 2' | grep "N" | wc -l```
Count: 3328051

Total: 7304664


## Script for part 1 

```FirstAssignment.py```



## Script for sbatch

```Part1.sh```



### Copied bioinfo.py module to this directory 

```scp bioinfo.py Talapas:/projects/bgmp/zgirard/bioinfo/Bi622/Demultiplex/Assignment-the-first```



## Test run to calculate mean quality score per base only for Read1 file

```
# Base Pair     Mean Quality Score
0       31.264301940663003
1       31.44230414073784
2       35.32827246196721
3       35.92871739094916
4       36.09835260322436
5       39.573363997339165
6       39.73660445426991
7       39.816290257915185
8       40.00983631690454
9       40.066305108014255
10      40.10437861747057
11      40.11889414780287
12      40.064250105372594
13      40.07103662198092
14      40.073171105034156
15      40.06607201851381
16      39.9829783604249
17      40.01358335127224
18      40.02123504565017
19      40.01391455314801
20      39.98320344985344
21      40.04378032744052
22      39.973750976729356
23      40.018210123760646
24      39.971674256067296
25      39.95728093192634
26      39.89541670347016
27      39.92127150709283
28      39.86349854459119
29      39.90736658101001
30      39.81043722526508
31      39.870599522388
32      39.852375548537275
33      39.895397397033726
34      39.91179441984523
35      39.904167636909385
36      39.8940528068339
37      39.86885067528549
38      39.88070280934528
39      39.86989397165538
40      39.85233052129154
41      39.837741423883685
42      39.82973229201909
43      39.825489189875306
44      39.74388511984836
45      39.78127342011759
46      39.79134226492084
47      39.78587638234381
48      39.780806010548176
49      39.77567073796273
50      39.77657128287746
51      39.76766332944465
52      39.687864192365005
53      39.72259013422378
54      39.72944116896192
55      39.726058457208154
56      39.67179063839349
57      39.68900449717738
58      39.676117080584355
59      39.66556348813431
60      39.64553211194038
61      39.64080384370144
62      39.589305976831426
63      39.60720985145262
64      39.59399951660956
65      39.581014640090295
66      39.56906488643318
67      39.55893994477335
68      39.544322313041576
69      39.44859068864032
70      39.46707217616147
71      39.46282866933408
72      39.45364293226201
73      39.432105023600556
74      39.40741819468797
75      39.383661906830355
76      38.64629549113497
77      39.05744404006825
78      39.17772193878081
79      39.22938594341392
80      39.32101735202107
81      39.32748463107315
82      39.30458273217514
83      39.28915057694875
84      39.271086665101066
85      39.26378823749097
86      39.22070142488686
87      39.0747698310351
88      39.11058770287364
89      39.11503113716906
90      39.09731839984742
91      39.05424496382604
92      39.01740024173927
93      38.96179684312923
94      38.93702517381195
95      38.94044496779854
96      38.92241449603119
97      38.701608888514855
98      38.753279268979526
99      38.83143660190091
100     37.66828871290474
```

Command being timed: "./FirstAssignment.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
        User time (seconds): 7423.05
        System time (seconds): 6.10
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 2:04:15
        Exit status: 0

Success!!!



Now Running the script on all 4 files with plots

sbatch Part1.sh 
Submitted batch job 7785413
Command being timed: "./FirstAssignment.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
        User time (seconds): 16485.78
        System time (seconds): 4.73
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 4:36:13
        Exit status: 0

Code ran successfully, but plots add new line each time. Added plt.close() between each plot to correct this.


New run:
sbatch Part1.sh 
Submitted batch job 7791841
Command being timed: "./FirstAssignment.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
        User time (seconds): 16778.78
        System time (seconds): 5.88
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 4:40:49
        Exit status: 0

Success! All plots are correct. Going to run one more time to remove plt.xlim for the index file plots. 



sbatch Part1.sh 
Submitted batch job 7794835


# Third Assignment: Demultiplexing script

Created a python script ```Demultiplex.py``` and an sbatch script ```Demultiplex.sh``` 

Running python script on test input files


```./Demultiplex.py -r1 ../TEST-input_FASTQ/read1_test.fq.gz -r4 ../TEST-input_FASTQ/read4_test.fq.gz -i1 ../TEST-input_FASTQ/read2_test.fq.gz  -i2 ../TEST-input_FASTQ/read3_test.fq.gz -ids /projects/bgmp/shared/2017_sequencing/indexes.txt ```



```@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1 GTAGCGTA-GTAGCGTA
@K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1 GTAGCGTA-GTAGCGTA
@K00337:83:HJKJNBBXX:8:1101:1489:1191 1:N:0:1 GTAGCGTA-TCGGATTC
@K00337:83:HJKJNBBXX:8:1101:1489:1191 4:N:0:1 GTAGCGTA-TCGGATTC
@K00337:83:HJKJNBBXX:8:1101:1671:1191 1:N:0:1 TACCGGAT-ATCCGGTN
@K00337:83:HJKJNBBXX:8:1101:1671:1191 4:N:0:1 TACCGGAT-ATCCGGTN
@K00337:83:HJKJNBBXX:8:1101:1286:1191 1:N:0:1 AAAAAAAA-TTTTTTTT
@K00337:83:HJKJNBBXX:8:1101:1286:1191 4:N:0:1 AAAAAAAA-TTTTTTTT
@K00337:83:HJKJNBBXX:8:1101:1448:1191 1:N:0:1 AGGATAGC-AGGATAGC
@K00337:83:HJKJNBBXX:8:1101:1448:1191 4:N:0:1 AGGATAGC-AGGATAGC
 -
 -
sample counts are: {'1': 1, '2': 0, '3': 0, '4': 0, '6': 0, '7': 0, '8': 0, '10': 0, '11': 0, '14': 0, '15': 0, '16': 0, '17': 0, '19': 0, '21': 0, '22': 0, '23': 0, '24': 0, '27': 0, '28': 0, '29': 0, '31': 0, '32': 0,'34': 1}
sample percentages are: {'1': 0.2, '2': 0.0, '3': 0.0, '4': 0.0, '6': 0.0, '7': 0.0, '8': 0.0, '10': 0.0, '11': 0.0, '14': 0.0, '15': 0.0, '16': 0.0, '17': 0.0, '19': 0.0, '21': 0.0, '22': 0.0, '23': 0.0, '24': 0.0, '27': 0.0, '28': 0.0, '29': 0.0, '31': 0.0, '32': 0.0, '34': 0.2}
index counts are: {'GTAGCGTA': 1, 'CGATCGAT': 0, 'GATCAAGG': 0, 'AACAGCGA': 0, 'TAGCCATG': 0, 'CGGTAATC': 0, 'CTCTGGAT': 0, 'TACCGGAT': 0, 'CTAGCTCA': 0, 'CACTTCAC': 0, 'GCTACTCT': 0, 'ACGATCAG': 0, 'TATGGCAC': 0, 'TGTTCCGT': 0, 'GTCCTAAG': 0, 'TCGACAAG': 0, 'TCTTCGAC': 0, 'ATCATGCG': 0, 'ATCGTGGT': 0, 'TCGAGAGT': 0, 'TCGGATTC': 0, 'GATCTTGC': 0, 'AGAGTCCA': 0, 'AGGATAGC': 1}
Number of matched indexes: 2
Number of unknown indexes: 2
Number of low quality indexes that don't contain N: 0
Number of hopped indexes: 1
Total number of records: 5```



Q score threshold is set to 20

Running sbatch script on actual fq files:

sbatch Demultiplex.sh 
Submitted batch job 7941864

1 min in, I have 2 x Hopped, 2 x Unknown, and 2 x 24 matched files!

8/7/2024
Script has been running for over 8 hours, with no failures.
It's working but obviously needs some optimization.
Probably due to opening each matched file every time I need to write to a file. I'm going to let the code run overnight, and see if it finishes. 
Will discuss with others how to avoid this tomorrow. 

Final Output:
Sample counts are: {'1': 8119243, '2': 5604966, '3': 6587100, '4': 8872034, '6': 10629633, '7': 5064906, '8': 34976387, '10': 76363857, '11': 17332036, '14': 4191388, '15': 7416557, '16': 7942853, '17': 11184304, '19': 15733007, '21': 8830276, '22': 3853350, '23': 42094112, '24': 10087503, '27': 6887592, '28': 11741547, '29': 4611350, '31': 3641072, '32': 11316780, '34': 8673180}
Sample percentages are: {'1': 0.02235186780137198, '2': 0.015430189620286607, '3': 0.0181339551475941, '4': 0.024424263579409737, '6': 0.02926284526686799, '7': 0.013943431590651461, '8': 0.09628823504772864, '10': 0.2102258592909307, '11': 0.04771422377684964, '14': 0.011538680450906186, '15': 0.020417408569412192, '16': 0.021866274999003087, '17': 0.030789826644966264, '19': 0.043312177327622776, '21': 0.02430930590470414, '22': 0.010608078831045791, '23': 0.11588297414428239, '24': 0.027770388631297677, '27': 0.018961194516999583, '28': 0.03232388860976273, '29': 0.01269481472421218, '31': 0.01002368816886957, '32': 0.03115452641301786, '34': 0.023876828514370542}
Matched index counts are: {'GTAGCGTA': 8119243, 'CGATCGAT': 5604966, 'GATCAAGG': 6587100, 'AACAGCGA': 8872034, 'TAGCCATG': 10629633, 'CGGTAATC': 5064906, 'CTCTGGAT': 34976387, 'TACCGGAT': 76363857, 'CTAGCTCA': 17332036, 'CACTTCAC': 4191388, 'GCTACTCT': 7416557, 'ACGATCAG': 7942853, 'TATGGCAC': 11184304, 'TGTTCCGT': 15733007, 'GTCCTAAG': 8830276, 'TCGACAAG': 3853350, 'TCTTCGAC': 42094112, 'ATCATGCG': 10087503, 'ATCGTGGT': 6887592, 'TCGAGAGT': 11741547, 'TCGGATTC': 4611350, 'GATCTTGC': 3641072, 'AGAGTCCA': 11316780, 'AGGATAGC': 8673180}
Number of matched indexes: 331755033
Number of unknown indexes: 30783962
Number of low quality indexes that don't contain N: 1044899
Number of hopped indexes: 7806306
Total number of records: 363246735


I seem to have read every record.
I need to multiply these percentage by 100. 

Command being timed: "./Demultiplex.py -ids /projects/bgmp/shared/2017_sequencing/indexes.txt -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
        User time (seconds): 9750.24
        System time (seconds): 10696.60
        Percent of CPU this job got: 53%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 10:42:39
        Exit status: 0


Job took almost 11 hours to run...

Talked to Lauren who showed me how to loop over an index list to open and close the matched files outside of the main loop.
Running new version of script

sbatch Demultiplex.sh 
Submitted batch job 8370024
Sample counts are: {'1': 8093048, '2': 5590177, '3': 6566013, '4': 8840103, '6': 10599469, '7': 5037387, '8': 34857831, '10': 76119286, '11': 17272877, '14': 4183306, '15': 7382032, '16': 7924413, '17': 11150072, '19': 15695547, '21': 8803030, '22': 3840626, '23': 42006477, '24': 10053967, '27': 6872523, '28': 11698720, '29': 4595129, '31': 3633844, '32': 11279057, '34': 8643481}
Sample percentages are: {'1': 2.227975428326974, '2': 1.5389476246772047, '3': 1.8075903696698057, '4': 2.433635914167267, '6': 2.9179805291298764, '7': 1.3867673167110506, '8': 9.596185633987874, '10': 20.95525676232162, '11': 4.755136202394221, '14': 1.1516431111211503, '15': 2.0322362980082946, '16': 2.1815510606034767, '17': 3.0695587669907067, '19': 4.32090518308444, '21': 2.423429903643869, '22': 1.0573050298717757, '23': 11.564171939494514, '24': 2.7678065709248565, '27': 1.8919710317561425, '28': 3.2205988031798825, '29': 1.2650159126688365, '31': 1.0003789848241857, '32': 3.1050676890461246, '34': 2.3795068660424437}
Matched index counts are: {'GTAGCGTA': 8093048, 'CGATCGAT': 5590177, 'GATCAAGG': 6566013, 'AACAGCGA': 8840103, 'TAGCCATG': 10599469, 'CGGTAATC': 5037387, 'CTCTGGAT': 34857831, 'TACCGGAT': 76119286, 'CTAGCTCA': 17272877, 'CACTTCAC': 4183306, 'GCTACTCT': 7382032, 'ACGATCAG': 7924413, 'TATGGCAC': 11150072, 'TGTTCCGT': 15695547, 'GTCCTAAG': 8803030, 'TCGACAAG': 3840626, 'TCTTCGAC': 42006477, 'ATCATGCG': 10053967, 'ATCGTGGT': 6872523, 'TCGAGAGT': 11698720, 'TCGGATTC': 4595129, 'GATCTTGC': 3633844, 'AGAGTCCA': 11279057, 'AGGATAGC': 8643481}
Number of matched indexes: 330738415
Percent dual matched: 91.05062293264659
Number of unknown indexes: 30783962
Percent unknown: 8.474669978795541
Number of low quality indexes (mean < 20) that don't contain N: 1044899
Percent Low Quality: 0.28765544169309604
Number of hopped indexes: 679459
Percent Hopped: 0.1870516468647681
Total number of records: 363246735
        Command being timed: "./Demultiplex.py -ids /projects/bgmp/shared/2017_sequencing/indexes.txt -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
        User time (seconds): 2833.96
        System time (seconds): 84.37
        Percent of CPU this job got: 97%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 50:07.88
        Exit status: 0


My counts match Lauren's!!!



