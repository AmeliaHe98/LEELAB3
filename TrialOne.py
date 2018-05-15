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

RefGenes = []
RefGenesOnly = []
variantOutputFile = pd.DataFrame()
InNucs = 2 #number of nucleotides on the boundaries of introns
RecNucs = 2 #number of nucleotides on the boundaries of introns that will be written to the output text file.



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

def vcf_file(vcf_file):
    '''
    read a vfc file and perform decision making for each variant file
    :param vcf_files: vfc file
    :return: cvs
    '''
    for f in vcf_file:
        vcf_reader = vcf.Reader(open(f, 'r'))
        for record in vcf_reader:
            findPositions(record)
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
    :return:
    '''
    d = {'Chromosome':[], 'Position':[],'Ref':[],'Alt':[],'GenomeLocation':[],'GeneName':[],
         'mRNAName':[],'GeneLocation':[],'MutationType':[],'RefCodonOrSpliceJun':[],
         'AltCodonOrSpliceJun':[],'RefAA':[],'AltAA':[]}
    global variantOutputFile
    variantOutputFile = pd.DataFrame(data=d)

def gff_file(gff_file):
    '''
    parse a gff file and return a dictionary of refGenes if (seqType==gene)
    :param gff_file: gff file
    :return: none
    '''
    with open(gff_file) as infile:
        for line in infile:
            parts = line.strip().split("\t")
            if(len(parts)==9):
                #Normalize data
                normalizedInfo = [
                        urllib.parse.unquote(parts[0]), #Chromosome
                        urllib.parse.unquote(parts[1]), #Source
                        urllib.parse.unquote(parts[2]), #SeqType
                        int(parts[3]), #Left.Pos
                        int(parts[4]), #Right.Pos
                        urllib.parse.unquote(parts[5]), #Score
                        urllib.parse.unquote(parts[6]), #Orientation
                        urllib.parse.unquote(parts[7]), #Frame
                        urllib.parse.unquote(parts[8]) #Name
                ]
                # store it in a list
                RefGenes.append(normalizedInfo)

    # gene = np.array(RefGenes)
    # print(np.where(gene=='gene'))

def findPositions(record):
    '''
    Determine which positions in the genome are affected by the given variant
    :param Record
    :return: list
    '''
    start = record.POS
    end = len(record.REF)-1+record.POS
    VarPos = []

    while start < end:
        VarPos.append(start)
        start = start + 1
    global variantOutputFile
    variantOutputFile.Chromosome = record.CHROM
    variantOutputFile.Position = record.POS
    variantOutputFile.Ref = record.REF
    variantOutputFile.Alt = record.ALT

    RefGene = SeqType(record,VarPos)
    if RefGene is None:
        return
    s = RefGene[8]
    startmRNAName = s.find('Name=') + 5
    startgeneName = s.find('ID=') + 3
    endgeneName = s.find(';')
    variantOutputFile.GeneName = s[startgeneName:endgeneName]
    # Returns from Name= to the end of the string
    variantOutputFile.mRNAName = s[startmRNAName:]
    if (RefGene[2] == "five_prime_UTR"):
        if((VarPos[1] > RefGene[4] and RefGene[6] == "-") or (VarPos[len(VarPos)] < RefGene[3] and RefGene[6]==  "+")):
#            global variantOutputFile
            variantOutputFile.GeneLocation = "five_prime_untranscribed"
        else:
#            global variantOutputFile
            variantOutputFile.GeneLocation="five_prime_UTR"
    elif (RefGene[2] == "three_prime_UTR"):
        if((VarPos[1] > RefGene[4] and RefGene[6] == "+") or (VarPos[len(VarPos)] < RefGene[3] and RefGene[6]==  "-")):

            variantOutputFile.GeneLocation = "five_prime_untranscribed"
        else:
            variantOutputFile.GeneLocation = "three_prime_UTR"
    #only cds and mRNA left






def SeqType(record, VarPos):
    '''
    if any genes that the given variant position is found within.
    :param Record, list
    :return: gene
    '''
    for gene in RefGenes:
        # return gene
        global variantOutputFile
        if (len(VarPos) != 0 and record.CHROM == gene[0]
                and gene[4] >= VarPos[0] >= gene[3]
                and gene[4] >= VarPos[len(VarPos) - 1] >= gene[3]):
            variantOutputFile.GenomeLocation = "Genic"
            return gene
        else:
            variantOutputFile.GenomeLocation = "Intergenic"
            return None




def main():
    fa_file(RefGenomeFile)
    gff_file(RefGenesFile)
    vcf_files(VCFDir)



if __name__ == '__main__':
    main()
