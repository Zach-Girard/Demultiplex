#!/usr/bin/env python

# Author: Zach Girard zgirard@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "1.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning



DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')





def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33


if __name__ == "__main__":
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Passed Convert Phred Score Tests!")



def qual_score(phred_score: str) -> float:
    """qual_score calculates the average quality score for an entire sequence"""
    sum = 0
    for score in phred_score:
        sum += convert_phred(score)
    return sum / len(phred_score)


    
if __name__ == "__main__":
    assert qual_score("A") == 32.0, "wrong average phred score for 'A'"
    assert qual_score("AC") == 33.0, "wrong average phred score for 'AC'"
    assert qual_score("@@##") == 16.5, "wrong average phred score for '@@##'"
    assert qual_score("EEEEAAA!") == 30.0, "wrong average phred score for 'EEEEAAA!'"
    assert qual_score("$") == 3.0, "wrong average phred score for '$'"
    print("Passed Quality Score Tests!")




def validate_base_seq(seq: str,RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)


if __name__ == "__main__":
	assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
	assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
	assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
	assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
	print("Passed DNA and RNA tests")






def gc_content(seq: str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(seq), "String contains invalid characters - Are you sure you used a nucleotide sequence?"
    seq = seq.upper()
    return (seq.count("G")+seq.count("C"))/len(seq)


if __name__ == "__main__":
    assert gc_content("AATAGCAT") == 0.25, "Wrong GC content calculated"
    assert gc_content("AATAAT") == 0.0, "GC Content does not work on sequences without G and C"
    assert gc_content("GCGCGC") == 1.0, "GC Content does not work if all "
    print("Passed GC Content Tests!")



def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    mid = (len(lst) - 1) // 2
    if len(lst) % 2 == 0:
        median = (lst[mid] + lst[mid + 1]) / 2
    else:
        median = lst[mid]
    return median

if __name__ == "__main__":
    assert calc_median([1,2,3]) == 2, "Calc Median does not work on odd number lists"
    assert calc_median([1,2,3,4]) == 2.5, "Calc Median does not work on even number lists"
    assert calc_median([42]) == 42, "Calc Median does not work on lists of length 1." 





def oneline_fasta(filename: str, output: str):
    '''This function takes a multiline fasta and converts each read into one line. It takes an input and output file name.'''
    with open(filename, "r") as fasta:
        with open(output, "w") as output:
            first_time = True
            for line in fasta:
                if first_time:
                    output.write(line)
                elif ">" in line:
                    output.write(f"\n{line}")
                else:
                    output.write(line.strip())
   




