#!/usr/bin/env python

import argparse, sys, os
from argparse import RawTextHelpFormatter
from collections import Counter

__author__ = "Swati Mishra"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: $"

# --------------------------------------
# This script calculates GC % N50 statistics , shortest and longest sequences from each fasta file and 300 bp detection in multifasta files
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
    Reporting sequences\n\
    author: " + __author__ + "\n\
    version: " + __version__ + "\n\
    description: report number of fasta sequences greater than 300 bp in multifasta files recursively in a directory\n\
                 report 10 shortest and 10 longest sequences in each fasta file\n\
                 report GC% and N50 statistics for each fasta file")
    parser.add_argument('-d', '--directory', type=str, default=None, help='Directory to search for fasta files')
    
    # parse the arguments
    args = parser.parse_args()

    if args.directory == None:
        parser.print_help()
        exit(1)
    # send back the user input
    return args

class Record(object):
    def __init__(self,seq_array):
        self.head = seq_array[0]
        self.seq = seq_array[1]
        self.length = len(seq_array[1])      
        
class Multifasta(object):        
    def __init__(self):
        fname = ''
        fastaList = list()
    def add_items(self,fname, fastaList):
        self.fname = fname
        self.fastaList = fastaList
    
    def GC(self):
        gcSum = sum([self.GCcounter(fas.seq) for fas in self.fastaList])
        Sum = sum([fas.length for fas in self.fastaList])
        try:
            return float(gcSum) / float(Sum) * 100
        except ZeroDivisionError:
            return 0
    def GCcounter(self, string):
        counter = Counter(string)
        return (counter['G'] + counter['C'])
                
    def N50(self):
        seqL = sorted([fas.length for fas in self.fastaList])
        Sum = sum(seqL)
        inSum = 0
        for l in seqL:
            inSum+=l
            if inSum >= float(Sum)/float(2):
                return l            


    def getCount_300_seq(self):
         return sum([1 for fasta in self.fastaList if fasta.length >=300])  
    def long10(self):
        if len(self.fastaList)  >= 10:
            return sorted(self.fastaList,key=lambda fasta: fasta.length,reverse=True)[0:10]
        else:
            return sorted(self.fastaList,key=lambda fasta: fasta.length,reverse=True)    
    def short10(self):
        
        if len(self.fastaList)  >= 10:
            return sorted(self.fastaList, key=lambda fasta: fasta.length)[0:10]
        else:
            return sorted(self.fastaList, key=lambda fasta: fasta.length)
        
        
def ProcessRecord(directory,output):
    if not os.path.exists(os.path.join(output,'longest10')):
        os.makedirs(os.path.join(output,'longest10'))
    if not os.path.exists(os.path.join(output,'shortest10')):    
        os.makedirs(os.path.join(output,'shortest10'))
    results = open(os.path.join(output,'results.txt'),'w')
    results.write('\t'.join(['sequenceName','CountOfSequences300bp','GC%','N50 Statistics'])+ '\n')
    for root, dirs, files in os.walk(directory):
        
        for file in files:
            if file.endswith(".fas") or file.endswith(".fa") \
              or file.endswith(".seq") or file.endswith(".fna") \
              or file.endswith(".ffn") or file.endswith(".faa") \
              or file.endswith(".frn"): 
                    fastafile = os.path.join(root, file)
                    linecount = 0
                    multifasta = Multifasta()
                    fasta_collection=[]
                    seq = ''
                    longest10 = open(os.path.join(output,'longest10',os.path.basename(fastafile)),'w')
                    shortest10 = open(os.path.join(output,'shortest10',os.path.basename(fastafile)),'w')
                    fp = open(fastafile,'r')
                    for line in fp:
                        linecount+=1
                        if linecount == 1:
                           if line[0] == '>':
                               item = []
                               seq = ''
                               item.append(''.join(line[1:]))
                        else:
                           if line[0] == '>':
                               if item:
                                    item.append(seq)
                                    fasta_collection.append(Record(item))

                               item = []
                               seq = ''
                               item.append(''.join(line[1:]))
                           else: 
                               seq = seq + line.upper()
                    
                    multifasta.add_items(os.path.basename(fastafile), fasta_collection)
                    results.write('\t'.join([multifasta.fname, str(multifasta.getCount_300_seq()), str(multifasta.GC()),str(multifasta.N50())]) + '\n')
                    longest10.write(''.join([''.join(['>'+fasta.head, fasta.seq]) for fasta in multifasta.long10()]))
                    shortest10.write(''.join([''.join(['>'+fasta.head, fasta.seq]) for fasta in multifasta.short10()]))
                    longest10.close()
                    shortest10.close()
    results.close()            

#
# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    ProcessRecord(args.directory,'../Results')
    
# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

