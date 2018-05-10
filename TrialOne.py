# -*- coding: utf-8 -*-
#!/usr/bin/python3
import vcf
import glob
import os
import pandas as pd
import numpy as np
import gffutils
import pprint
from Bio import SeqIO
from collections import namedtuple
import urllib.request, urllib.parse, urllib.error

InNucs = 2 #number of nucleotides on the boundaries of introns
RecNucs = 2 #number of nucleotides on the boundaries of introns that will be written to the output text file.
# Note that this must be a number greater than or equal to InNucs. Default 5.

#Initialized GeneInfo named tuple. Note: namedtuple is immutable
gffInfoFields = ["Chromosome", "Source", "SeqType", "LeftPos", "RightPos", "Score", "Orientation", "Frame", "Name"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)

VCFDir = '/Users/Amelia/Desktop/lab' #directory containing the VCFs
VarEffDir = '/Users/Amelia/Desktop/lab' #directory of output file containing information on each variant
RefGenomeFile = '/Users/Amelia/Desktop/lab/Cr53C.fa' #Reference genome from which variant calling was done.
RefGenesFile = '/Users/Amelia/Desktop/lab/Creinhardtii_281_v5.5.gene_exons.gff3' #gff3 file containing the gene location and structure information requried for characterizing variants.


def fa_file(fa_file):
    '''
    read a fasta file and return  a dictionary with
    chromosome as key and sequence in the form of string as value
    :param fa_file: fasta file
    :return: dictionary
    '''
    fasta_sequences = SeqIO.parse(open(fa_file), 'fasta')
    mydict = {}
    for seq_record in fasta_sequences:
        mydict[seq_record.id] = str(seq_record.seq).upper()
    return mydict

def gff_file(gff_file):
    '''
    parse a gff file and return a dictionary
    :param gff_file: gff file
    :return: dictionary with name as key and everything else as value
    '''
    ret = {}
    with open(gff_file) as infile:
        for line in infile:
            parts = line.strip().split("\t")
            if(len(parts)==9):
            # Normalize data
                normalizedInfo = [
                    urllib.parse.unquote(parts[0]),
                    urllib.parse.unquote(parts[1]),
                    urllib.parse.unquote(parts[2]),
                    int(parts[3]),
                    int(parts[4]),
                    urllib.parse.unquote(parts[5]),
                    urllib.parse.unquote(parts[6]),
                    urllib.parse.unquote(parts[7]),
                    urllib.parse.unquote(parts[8])
                ]
                ret[parts[8]] = normalizedInfo
    return ret


def vcf_file(vcf_files):
    '''
    read a vfc file and perform decision making for each variant file
    :param vcf_files: vfc file
    :return: cvs
    '''
    variantOutputFile = cvs_files()
    for f in vcf_files:
        with open(f, mode='r') as vcf:
            vcf.read()
            #do something here with each vcf file


def vcf_files(vcf_dir):
    '''
    read the directory containing vfc files
    :param vcf_dir: vfc directory
    :return: list of vfc files
    '''
    files = []
    for file in glob.glob(os.path.join(vcf_dir, '*.vcf')):
        files.append(file)
    return vcf_file(files)


def cvs_files():
    '''
    create a cvs file and store the variant information. Can be mutated
    :param not sure if needed
    :return: DataFrame
    '''
    d = {'Chromosome':[], 'Position':[],'Ref':[],'Alt':[],'GenomeLocation':[],'GeneName':[],
         'mRNAName':[],'GeneLocation':[],'MutationType':[],'RefCodonOrSpliceJun':[],
         'AltCodonOrSpliceJun':[],'RefAA':[],'AltAA':[]}
    df = pd.DataFrame(data=d)
    return df


def main():
    fa_file(RefGenomeFile)
    gff_file(RefGenesFile)
    vcf_files(VCFDir)



if __name__ == '__main__':
    main()
