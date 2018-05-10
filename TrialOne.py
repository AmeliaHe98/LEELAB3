# -*- coding: utf-8 -*-
#!/usr/bin/python3
import vcf
from BCBio import GFF
import gffutils
from Bio import SeqIO
from collections import namedtuple
import urllib.request, urllib.parse, urllib.error

#Set the number of nucleotides on the boundaries of introns that should be genotyped
# (i.e. that should be checked for mutations). If nucleotides around the boundaries
# of introns are mutated, this will affect splice sites and therefore very strongly
# affect the protein structure of the gene in question. Such mutations would be equivalent
# to nonsense mutations in terms of their effect. The default is 2.
InNucs = 2

#Set the number of nucleotides on the boundaries of introns that will be written to the
# output text file. Note that this must be a number greater than or equal to InNucs. Default 5.
RecNucs = 2

##Declare directories that need to be manually set - These are directories that contain either files or programs required for the code

#Set the directory containing the VCFs for which the genomic effect of their variants will be determined.
VCFDir = '/Users/Amelia/Desktop/lab'

#Set the directory into which the text file containing information on each variant will be written.
VarEffDir = '/Users/Amelia/Desktop/lab'
##Declare reference files that need to be manually set

#Reference genome from which variant calling was done.
RefGenomeFile = '/Users/Amelia/Desktop/lab/Cr53C.fa'

#gff3 file containing the gene location and structure information requried for characterizing variants.
RefGenesFile = '/Users/Amelia/Desktop/lab/Creinhardtii_281_v5.5.gene_exons.gff3'

# gffInfoFields = ["Chromosome", "Source", "SeqType", "LeftPos", "RightPos", "Score", "Orientation", "Frame", "Name"]
# GFFRecord = namedtuple("GFFRecord", gffInfoFields)

def fa_file(fa_file):
    fasta_sequences = SeqIO.parse(open(fa_file), 'fasta')
    mydict = {}
    for seq_record in fasta_sequences:
        #chromosome: sequence(as a string)
        mydict[seq_record.id] = str(seq_record.seq).upper()
    return mydict

# def parseGFFAttributes(attributeString):
#     ret = {}
#     for attribute in attributeString.split(";"):
#         key, value = attribute.split("=")
#         ret[urllib.parse.unquote(key)] = urllib.parse.unquote(value)
#     return ret

def gff_file(gff_file):
    open(gff_file).read()
    # I was going to write a dictionary to store gff data before I
    # found its library. I'm keeping this code for now in case we might
    # use this later
    # with open(gff_file) as infile:
    #     for line in infile:
    #         parts = line.strip().split("\t")
    #         assert len(parts) == len(gffInfoFields)
    #         # Normalize data
    #         normalizedInfo = {
    #             "Chromosome": urllib.parse.unquote(parts[0]),
    #             "Source": urllib.parse.unquote(parts[1]),
    #             "SeqType": urllib.parse.unquote(parts[2]),
    #             "LeftPos": int(parts[3]),
    #             "RightPos": int(parts[4]),
    #             "Score": float(parts[5]),
    #             "Orientation": urllib.parse.unquote(parts[6]),
    #             "Frame": urllib.parse.unquote(parts[7]),
    #             "Name": parseGFFAttributes(parts[8])
    #         }
    #         #    yield normalizedInfo
    #         yield GFFRecord(**normalizedInfo)


def main():
    fa_file(RefGenomeFile)
    gff_file(RefGenesFile)


if __name__ == '__main__':
    main()
