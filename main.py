from Bio import SeqIO

from vcf import parser


def main():
    ### This code is for determining the genomic effect of variants present in a vcf file using an
    ### input gff3 file containing information on the location and structure of genes.
    ### Based on this gff3 file and the vcf, the code goes through each line of the vcf and determines the genomic effect
    ### of the given variant. This includes whether or not the variant is found within a gene, and if it is found
    ### within a gene, whether or not it is found within a protein coding sequence (and what sort of change it makes
    ### to the protein-coding sequence).

    ## Begin section containing variables/directories/files that need to be manually set.

    # Declare strings that need to be manually set - These strings are typically required for naming files

    # Set the aligner through which the sequencing was aligned
    Al = "Bowtie"

    # Set the name of the reference file used for aligning
    RefName = "Cr53C"

    # Set the number of nucleotides on the boundaries of introns that should be genotyped (i.e. that should be checked for mutations). If nucleotides around the boundaries of introns are mutated, this will affect splice sites and therefore very strongly affect the protein structure of the gene in question. Such mutations would be equivalent to nonsense mutations in terms of their effect. The default is 2.
    InNucs = 2

    # Set the number of nucleotides on the boundaries of introns that will be written to the output text file. Note that this must be a number greater than or equal to InNucs. Default 5.
    RecNucs = 2



if __name__ == "__main__":
    main()