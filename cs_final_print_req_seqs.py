# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 10:46:29 2017

@author: Zachary_Roga
"""

import numpy as np
from numpy import log1p
import os
import matplotlib.pyplot as plt
import math
from math import log
import argparse
import random
#note: check that the sequence names don't leave a space
def main():
    src_gene_files = r"C:\Users\Zachary_Roga\Documents\past semesters\cs4775\Final_project_data\\"

    if not os.path.exists(src_gene_files):
        src_gene_files = input("please link the Final_project_data folder")
    gene_results_path = src_gene_files + "Final_project_datatest"    

    strain_genome_locations=[src_gene_files+ "AWRI796\AWRI796_ADVS01000000_cds.fsa",
                             src_gene_files+ "AWRI1631\AWRI1631_ABSV01000000_cds.fsa",
                             src_gene_files+ "BC187\BC187_JRII00000000_cds.fsa",
                             src_gene_files+ "CBS7960\CBS7960_AEWL01000000_cds.fsa",
                             src_gene_files+ "CLIB382\CLIB382_AFDG01000000_cds.fsa",
                             src_gene_files+ "EC1118\EC1118_PRJEA37863_cds.fsa",
                             src_gene_files+ "Fosters_B\FostersB_AEHH01000000_cds.fsa",
                             src_gene_files+ "Fosters_O\FostersO_AEEZ01000000_cds.fsa",
                             src_gene_files+ "KyoKai7\Kyokai7_BABQ01000000_cds.fsa",
                             #"C:\Users\Zachary_Roga\Documents\cs4775\Final_project_data\LalvinQA23\LalvinQA23_ADVV01000000_cds.fsa",
                             src_gene_files+ "M22\M22_ABPC01000000_cds.fsa",
                             src_gene_files+ "T73\T73_AFDF01000000_cds.fsa",
                             src_gene_files+ "VL3\VL3_AEJS01000000_cds.fsa"] #list of all the genome file locations
    seqs=[]
    genes_src=["YKR069W_", "YGR155W_"] #"YFR030W_",, "YBR213W_", "YMR307W_","YAL022C_", "YLR289W_" #"YOL030W_", "YAL022C_", "YLR289W_" /// "YBR213W_", "YGR155W_"
    #for genes in genes_src:
    for lmao in range(0,1):
        seqs=[] 
        for p1 in strain_genome_locations:
            with open(p1, "r") as infile:
                checker=0
                named=0
                k=[]
                for line in infile: #YAL022C, YAL020C#line[1:1+len("YLR289W_")]=="YLR288W_" or  line[1:1+len("YLR211W_")]=="YLR289W_" or line[1:1+len("YLR211W_")]=="YLR293C_"
                    for genes in genes_src:
                    #for lmao in range(0,1):
                        if  line[1:1+len(genes)]==genes:
                                if named==0:
                                    names=['>']
                                    iters=1+len("YAL027W_")
                                    while line[iters]!=' ' and iters < min(20, len(line)):
                                        names.append(line[iters])
                                        iters+=1
                                    seqs.append(''.join(names))
                                    named=1
                                    
                                checker=1
                                #print line[0:min(20, len(line))]
                        elif checker!=0:
                                if line[0]!='>':
                                    k.append(line.strip())
                                else:
                                    checker=0
                if len(k)>0:
                    seqs.append(''.join(k))
                k=[]
                checker=0
#                            i=i+len("YLR289W")
#                k=''.join(k)
#                seqs.append(k)
        file_pathsout=gene_results_path+ "_tree_data"+".txt"
        with open(file_pathsout, "w+") as outfile:
            for row in seqs:
                outfile.write("%s\n" % row)
        outfile.close()
    
    #print seqs
    
if __name__ == "__main__":
    main()