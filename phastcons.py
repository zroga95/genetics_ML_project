# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 22:46:46 2017

@author: Zachary_Roga
"""

import numpy as np
from  numpy import log1p
import os
import textwrap
import matplotlib.pyplot as plt
import math
from math import log
import argparse
import random
import re
#create a binary tree with branch lengths that are scaled
def sumLogProb(a, b):
    if a > b: 
        return a + log1p(np.exp(b - a))
    else: 
        return b + log1p(np.exp(a - b))

def genBaseTree(scale, B):
    return BinaryTree("root", BinaryTree("CLIB382-FostersO-FostersB-VL3-CBS7960-AWRI1631-Kyokai7-M22",
        BinaryTree("CLIB382-FostersO-FostersB-VL3-CBS7960-AWRI1631-Kyokai7",
        BinaryTree('CLIB382-FostersO-FostersB-VL3-CBS7960-AWRI1631', 
        BinaryTree('CLIB382-FostersO-FostersB-VL3-CBS7960',BinaryTree('CLIB382-FostersO-FostersB-VL3',
        BinaryTree('CLIB382-FostersO-FostersB', BinaryTree('CLIB382'),BinaryTree('FostersO-FostersB', BinaryTree('FostersO'),
        BinaryTree('FostersB'), scale*0.001010, scale*0.000001), scale*0.000001, scale*0.000001),
        BinaryTree("VL3"), scale*0.000001, scale*0.000001), BinaryTree("CBS7960"),scale*0.000001,scale*0.000001),
        BinaryTree("AWRI1631"), scale*0.000001, scale*0.000001), BinaryTree("Kyokai7"),scale*0.000001, scale*0.001907), 
        BinaryTree("M22"), scale*0.000001, scale*0.000953),
        BinaryTree("EC1118-BC187-AWRI796", BinaryTree("EC1118"),BinaryTree("BC187-AWRI796",BinaryTree("BC187"),
        BinaryTree("AWRI796"),scale*0.000001,scale*0.001010), scale*0.000001, scale*0.000001), scale*0.000953, scale*0.000001)

def genTestTree(scale, B):
    return BinaryTree('root', BinaryTree('HD', BinaryTree('H'), BinaryTree('D'), scale*0.07517, scale*0.11761),
        BinaryTree('MR', BinaryTree('M'), BinaryTree('R'), scale*0.03059, scale*0.03161), scale*0,scale*0.14289)
    
#create a post-traversal order for the tree, with the lengths of the branch after the entry
def Transverse_tree(tree, arrays):
    if tree.getNodeValue() is None: 
        return
    if tree.getLeftChild() is None: 
        arrays.append(tree)
        return
    if tree.getRightChild() is None: 
        arrays.append(tree)
        return
    Transverse_tree(tree.getLeftChild(), arrays)
    arrays.append(tree.getlChildLength())
    Transverse_tree(tree.getRightChild(),  arrays)
    arrays.append(tree.getrChildLength())
    arrays.append(tree)

def pruneTree(tree, names):
    checks=[0,0]
    if (tree.getLeftChild()!=None and tree.getRightChild()!=None):
        #check if the branch is invalid for current genes
        if tree.getLeftChild().getNodeValue() not in names and "-" not in tree.getLeftChild().getNodeValue() :
            #if one of the children of the branch is invalid, make this branch the valid child
            tree.setLeftChild(None)
            if tree.getRightChild().getRightChild()!=None:#the other part of the tree is not a leaf    
                tree.setlChildLength(tree.getRightChild().getlChildLength()+tree.getrChildLength())
                tree.setrChildLength(tree.getRightChild().getrChildLength()+tree.getrChildLength())
                tree.setLeftChild(tree.getRightChild().getLeftChild())
                tree.setRightChild(tree.getRightChild().getRightChild())
            checks[0]=1
        if tree.getRightChild().getNodeValue() not in names and "-" not in tree.getRightChild().getNodeValue() :
            tree.setRightChild(None)
            checks[1]=1
            if tree.getLeftChild().getRightChild()!=None:#the other part of the tree is not a leaf    
                tree.setrChildLength(tree.getLeftChild().getrChildLength()+tree.getrChildLength())
                tree.setlChildLength(tree.getLeftChild().getlChildLength()+tree.getlChildLength())
                tree.setRightChild(tree.getLeftChild().getRightChild())
                tree.setLeftChild(tree.getLeftChild().getLeftChild())
        if checks[0]==1:
            if checks[1]==1:
                return
            else:
                return pruneTree(tree.getRightChild(), names)
        if checks[1]==1:
            return pruneTree(tree.getLeftChild(), names)
        else:
            return pruneTree(tree.getLeftChild(), names),pruneTree(tree.getRightChild(), names)        
    else:
        return
def likelihoods(seqs, u_t, l, names2, max_len):
    #check for valid sequences
    order_seqs=names2
    probs_matrix=[]
    bp_order=["A","T","G","C"]
    pi=[0.25, 0.25, 0.25, 0.25] #background probabilities vector
    for j in range(0,int(len(l)/2)):
        #for each part of the tree, iterate through performing the algorithm
        if l[2*j].getNodeValue() in order_seqs:
            #leaf cases. for each indice u, determine if the probability is 0 or 1 and store that 4*n matrix in probs_matrix
            rows=[]
            for p in range(0,len(seqs[order_seqs.index(l[2*j].getNodeValue())])):
                columns=[]
                if seqs[order_seqs.index(l[2*j].getNodeValue())][p] not in bp_order: #if missing data, give this segment no impact on probability
                    columns=[1,1,1,1]
                else:
                    for b in range(0,4):
                        if seqs[order_seqs.index(l[2*j].getNodeValue())][p]==bp_order[b]:
                            columns.append(1)
                        else:
                            columns.append(0)
                rows.append(columns)
            probs_matrix.append(rows)
        elif "-" in l[2*j].getNodeValue() or l[2*j].getNodeValue()=="root":  #nonbase nodes
            rows=[]
            leftchildindx=l.index(l[2*j].getLeftChild())
            rightchildindx=l.index(l[2*j].getRightChild())
            #predetermine the possible values for probabilities of the array including the branch length
            cond_prob_array=[1+3*math.exp(-4*u_t*l[2*j].getrChildLength()), 1-math.exp(-4*u_t*l[2*j].getrChildLength()),
                1+3*math.exp(-4*u_t*l[2*j].getlChildLength()), 1-math.exp(-4*u_t*l[2*j].getlChildLength())]
            gap_data=[1,1,1,1]
            for p in range(0,max_len):     #for each indice
                columns=[]
                if p<len(probs_matrix[int(rightchildindx/2)]) and p<len(probs_matrix[int(leftchildindx/2)]):
                    for b in range(0,4):
                        #cycle through the different bases. for each base, take the summation of the child sequence's likelihood 
                        #multiplied by the mutation expression vector. multiply this value from both child sequences
                        #append that from each base to form a 4*n matrix for this node
                        b1=0
                        b2=0
                        for c in range(0,4):
                            if b==c:
                            #check through child sequences for bp matches
                                b1=b1+(0.25*cond_prob_array[0])*probs_matrix[int(rightchildindx/2)][p][c]       #equation for match case 
                            if b!=c:
                                b1=b1+(0.25*cond_prob_array[1])*probs_matrix[int(rightchildindx/2)][p][c]        #not match case
                            if b==c:    
                                b2= b2+(0.25*cond_prob_array[2])*probs_matrix[int(leftchildindx/2)][p][c]       #equation for match case 
                            if b!=c:
                                b2=b2+(0.25*cond_prob_array[3])*probs_matrix[int(leftchildindx/2)][p][c]        #not match case
                        columns.append(b1*b2)
                else:
                    #if there is not data for both sequences at this point, use the availabile datad, else disregard this 
                        if p>=len(probs_matrix[int(rightchildindx/2)]) and p>=len(probs_matrix[j-1]):
                            columns=gap_data    
                        elif p>=len(probs_matrix[int(rightchildindx/2)]):
                           columns= probs_matrix[int(leftchildindx/2)][p]
                        else:
                           columns= probs_matrix[int(rightchildindx/2)][p] 

                rows.append(columns)
            probs_matrix.append(rows)                
    #terminating of sequence, sum up the probability for all bases, add log likelihoods for each sequence
    seq_logs=[]
    indexProb=[]
    for j in range(0,int(len(l)/2)):
        curr_seq=[]
        tot_prob=0                
        for p in range(0,len(probs_matrix[j])):
            curr_prob=0
            for b in range(0,4):
                curr_prob+=0.25*probs_matrix[j][p][b]
            curr_seq.append(log(curr_prob))
            tot_prob+=log(curr_prob)
        indexProb.append(curr_seq)
        seq_logs.append(tot_prob)
    return seq_logs, probs_matrix, indexProb

class BinaryTree():
    
    def __init__(self, name, left=None, right=None, lChildLength=0, rChildLength=0):
      self.root=name
      self.left = left
      self.right = right
      self.rChildLength= rChildLength
      self.lChildLength= lChildLength
    def getLeftChild(self):
        return self.left
    def getRightChild(self):
        return self.right
    def setNodeValue(self,value):
        self.root = value
    def getrChildLength(self):
        return self.rChildLength
    def getlChildLength(self):
        return self.lChildLength    
    def getNodeValue(self):
        return self.root
    def setLeftChild(self,value):
        self.left = value
    def setRightChild(self,value):
        self.right = value
    def setlChildLength(self,value):
        self.lChildLength = value
    def setrChildLength(self,value):
        self.rChildLength = value
def Forward(non_matrix, con_matrix, mu, vu):    #states is a 4x4 matrix of substitution rates
    bp_order=["A","T","G","C"]
    pi=[0.25, 0.25, 0.25, 0.25]
    forward_matrix=[]
    forward_matrix=[[log(0.5),log(0.5)]]
    if len(non_matrix)!=len(con_matrix):
        print ("ERROR, ERROR. MATRICES ARE NOT SAME LENGTH")
    else:
        for i in range(1,len(non_matrix)):
            column=[]
            column.append(non_matrix[i]+sumLogProb(log(1-vu)+forward_matrix[i-1][0],
                log(mu)+forward_matrix[i-1][1])) #non
            column.append(con_matrix[i]+sumLogProb(log(1-mu)+forward_matrix[i-1][1],
                log(vu)+forward_matrix[i-1][0])) #con
            forward_matrix.append(column)
    return forward_matrix
def Backward(non_matrix, con_matrix, mu, vu):    #states is a 4x4 matrix of substitution rates
    bp_order=["A","T","G","C"]
    pi=[0.25, 0.25, 0.25, 0.25]
    
    if len(non_matrix)!=len(con_matrix):
        print ("ERROR, ERROR. MATRICES ARE NOT SAME LENGTH")
    else:
        backward_matrix=np.zeros((len(con_matrix), 2))
        backward_matrix[len(con_matrix)-1]=[1,1]
        for i in range(1,len(non_matrix)):
            column=[]
            column.append(sumLogProb(log(1-vu)+backward_matrix[len(con_matrix)-i][0]
                +non_matrix[len(con_matrix)-i],log(mu)+backward_matrix[len(con_matrix)-i][1]
                +con_matrix[len(con_matrix)-i])) #non
            column.append(sumLogProb(log(1-mu)+backward_matrix[len(con_matrix)-i][1]
                +con_matrix[len(con_matrix)-i],log(vu)+backward_matrix[len(con_matrix)-i][0]
                +non_matrix[len(con_matrix)-i])) #con
            backward_matrix[len(con_matrix)-1-i]=column
    return backward_matrix

def ForBack(forward_matrix, backward_matrix):
    if len(forward_matrix)!=len(backward_matrix):
        print ("ERROR, ERROR. MATRICES ARE NOT SAME LENGTH")
    else:
        ForBackMat=[]
        for i in range(0,len(backward_matrix)):
            denom=sumLogProb(forward_matrix[i][0]+backward_matrix[i][0],forward_matrix[i][1]
                +backward_matrix[i][1])
            ForBackMat.append([forward_matrix[i][0]+backward_matrix[i][0]-denom,forward_matrix[i][1]
                +backward_matrix[i][1]-denom])
        return ForBackMat
def Viterbi(non_matrix, con_matrix, mu, vu):
    viterbiMat=[]
    backPoint=[[2,2]]
    viterbiMat.append([log(0.5)+non_matrix[0],log(0.5)+con_matrix[0]]) #non,con
    sums=0
    for i in range(1,len(non_matrix)):
        column=[]
        nonSt=[log(1-vu)+viterbiMat[i-1][0],log(mu)+viterbiMat[i-1][1]]
        conSt=[log(vu)+viterbiMat[i-1][0],log(1-mu)+viterbiMat[i-1][1]]
        column.append(non_matrix[i]+max(nonSt))
        column.append(con_matrix[i]+max(conSt))
        viterbiMat.append(column)
        backPoint.append([nonSt.index(max(nonSt)),conSt.index(max(conSt))])       #0 or 1 to indicate the state that caused the value
    #create a path from the states
    stPath=[]
    curr=viterbiMat[len(non_matrix)-1].index(max(viterbiMat[len(non_matrix)-1]))
    for i in range(0,len(non_matrix)):
        if curr==0:
            stPath.insert(0, False)
        else:
            stPath.insert(0,True)
        indicer=len(non_matrix)-i-1
        #print pointers[999-i]
        curr=backPoint[indicer][int(curr)]
    return stPath                
        
    #use the backpointers to create a path
    #print viterbiMat
def phastconners(genes_loc, src_gene_files, mu):
    #    genes_loc='C:\Users\Zachary_Roga\Documents\cs4775\\apoe.fa'
        trees_loc= src_gene_files+"newick.txt"
        seqs=[]
        checker=0
        sequenced=[]
        names=[]
        names2=[]
        #split open the file and separate the names from the actual sequences
        with open(genes_loc, "r") as infile:
            for line in infile:
                if line[0]=='>':
                    names.append(line.strip())
                    names2.append(line.strip()[1:])
                    if checker==1:
                        seqs.append(''.join(sequenced))
                    sequenced=[]
                    checker=1
                elif checker==1:
                    sequenced.append(line.strip())
            seqs.append(''.join(sequenced))
        infile.close()
        #open the newick tree and get the data
        origin_tree_data=[]
        with open(trees_loc, "r") as infile2:
            for line2 in infile2:
                origin_tree_data.append(line2.strip())
        origin_tree_data=list("".join(origin_tree_data))
        #create an array wtih the lengths of sequences
        Q=[]    #substition
        bp_order=["A","T","G","C"]
        pi=[0.25, 0.25, 0.25, 0.25] #background probabilities vector
        p=0.325       #scaling parameter
        vu=mu
        lens=[]
        for i in seqs:
            lens.append(len(i))
        #make the info into a parseable tree
        #pico=0  #len(origin_tree_data)
        #create the binary tree, remove non-sequenced nodes, get a traversal order.
        mu1=2
        vu1=mu1
        while abs(mu-mu1)>0.001 or abs(vu-vu1)>0.001:
            vu1=vu
            mu1=mu
            B=[0.001010, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001,0.000001, 
            0.000001,0.000001, 0.000001, 0.001907, 0.000001, 0.000953, 0.000001,0.001010, 0.000001, 0.000001]
            B_order=['FostersO,','FostersB','CLIB382','FostersO-FostersB', ]
            BaseTree=BinaryTree("root", BinaryTree("CLIB382-FostersO-FostersB-VL3-CBS7960-AWRI1631-Kyokai7-M22", 
                BinaryTree("CLIB382-FostersO-FostersB-VL3-CBS7960-AWRI1631-Kyokai7",
                BinaryTree('CLIB382-FostersO-FostersB-VL3-CBS7960-AWRI1631', 
                BinaryTree('CLIB382-FostersO-FostersB-VL3-CBS7960',BinaryTree('CLIB382-FostersO-FostersB-VL3',
                BinaryTree('CLIB382-FostersO-FostersB', BinaryTree('CLIB382'),BinaryTree('FostersO-FostersB',
                BinaryTree('FostersO'), BinaryTree('FostersB'), 0.001010, 0.000001), 0.000001, 0.000001), 
                BinaryTree("VL3"), 0.000001, 0.000001), BinaryTree("CBS7960"),0.000001,0.000001), BinaryTree("AWRI1631"),
                0.000001, 0.000001), BinaryTree("Kyokai7"),0.000001, 0.001907), BinaryTree("M22"), 0.000001, 0.000953),
                BinaryTree("EC1118-BC187-AWRI796", BinaryTree("EC1118"),BinaryTree("BC187-AWRI796",BinaryTree("BC187"),
                BinaryTree("AWRI796"),0.000001,0.001010), 0.000001, 0.000001), 0.000953, 0.000001)   
            
            BaseTree=genBaseTree(1, B) if gene_parser_type == "Y" else genTestTree(p, B)

            names2=["H","M","R","D"] if gene_parser_type != "Y" else names2
            pruneTree(BaseTree, names2)
            newick=[]
            Transverse_tree(BaseTree,newick)
            newick.append(0.0)
            nonliker_seq, nonliker_matrix, nonStates=likelihoods(seqs, vu, newick, names2, max(lens))
            #conserve_start_low=-999999.9
            #ender=-999999.0
            #while ender>conserve_start_low:#attempt at the second part of algorithm, however, their technique isnt in rightup and this didnt converge
            #conserve_start_low=ender

            conserveTree=genBaseTree(p, B) if gene_parser_type == "Y"  else genTestTree(p, B)

            pruneTree(conserveTree, names2)
            newickCon=[]
            Transverse_tree(conserveTree,newickCon)
            newickCon.append(0.0)        
            coliker_seq, coliker_matrix, conStates=likelihoods(seqs, mu, newickCon, names2, max(lens))
            #ender=coliker_seq[len(coliker_seq)-1]
        #run viterbi algorithm to potimize tranistion parameters by EM
            paths=Viterbi(nonStates[len(nonStates)-1], conStates[len(nonStates)-1], mu, vu)
            cons_path=[]
            nums_cons=0
            
            for i in range(0,len(paths)):
                if paths:
                    nums_cons+=1
            nums_nons=len(paths)-nums_cons
            switch_nc=0
            switch_cn=0
            if paths[0]:
                cons_path.append('(')
                cons_path.append(str(i+1))
                cons_path.append(',')
            prev_stat=paths[0]
            for i in range(1,len(paths)):
                if prev_stat!= paths[i]:
                    if paths[i]:
                        switch_nc+=1
                        cons_path.append('(')
                        cons_path.append(str(i+1))
                        cons_path.append(',')
                    else:
                        switch_cn+=1
                        cons_path.append(str(i))
                        cons_path.append(')')
                prev_stat=paths[i]
            if nums_cons<0.05*len(paths):
                nums_cons=0.05*len(paths)
            if nums_nons<0.05*len(paths):
                nums_nons=0.05*len(paths)
            if switch_cn<nums_cons*0.01:
                switch_cn=nums_cons*0.01
            if switch_nc<nums_nons*0.01:
                switch_nc=nums_nons*0.01
            mu=switch_cn/nums_cons
            vu=switch_nc/nums_nons
            
            if len(cons_path)>=1:
                if cons_path[len(cons_path)-1]!=')':
                    cons_path.append('1000')
                    cons_path.append(')')                
        print (''.join(cons_path))
        #now we have to run the forward backward algorithm to parse ideal path between the matrices
        fors=[]
        backs=[]
        test_data=[]
        for i in range(0,len(nonStates)):
            fors.append(Forward(nonStates[i], conStates[i], mu, vu))    
            backs.append(Backward(nonStates[i], conStates[i], mu, vu))    
            test_data.append(ForBack(fors[i],backs[i]))
        model_check=[]
        tester_data=[]
        for i in range(0,len(test_data[len(test_data)-1])):
            model_check.append(math.exp(test_data[len(test_data)-1][i][1]))
            tester_data.append(math.exp(test_data[len(test_data)-1][i][0]))
        print ("final values of mu and vu: "+str(mu)+ ", "+ str(vu))
        print (nonliker_seq[len(nonliker_seq)-1], "non-conserved sequence log likelihood")#, genes
        print (coliker_seq[len(coliker_seq)-1], "conserved sequence log likelihood")#, genes    
        return model_check, lens

    
def main():
    src_gene_files = "C:/Users/Zachary_Roga/Documents/past semesters/cs4775/Final_project_data/"
    dashes = [20, 10, 10, 5]
    if not os.path.exists(src_gene_files):
        src_gene_files = str(input("please link the Final_project_data folder"))
    gene_results_path = src_gene_files + "Final_project_datatest"
    if gene_parser_type == 'Y':
        gene_list = ["YMR307W_", "YOL030W_","YFR030W_", "YGR155W_",
        "YBR213W_", "YKR069W_", "YAL022C_", "YLR289W_"]
        for genes in gene_list: #
            genes_loc=gene_results_path+ genes+".txt"
            model_check, lens = phastconners(genes_loc, src_gene_files, 0.05)
            print ("testing model on gene "+str(genes))
    else:
        genes_loc = src_gene_files + 'apoe.fa'
        model_check, lens = phastconners(genes_loc, src_gene_files, 0.25)
    
    plt.figure()
    plt.plot(model_check)
    for i in range(0,len(lens)):
        line1=plt.axvline(x=lens[i], color="r", lw=0.5, )
        line1.set_dashes(dashes)
    plt.xlabel( "length of bp") #str(genes)+
    plt.ylabel("probability of conservation model")
    plt.show()
        
#plots looked very similar so checked that they aren't actually symmetric
#    equality=[]
#    for i in range(0,len(test_data[len(test_data)-1])/2):
#        equality.append(test_data[len(test_data)-1][i][1]==test_data[len(test_data)-1][-i][1])
#    plt.plot(equality)
if __name__ == "__main__":
    gene_parser_type = ""
    while gene_parser_type != "Y" and gene_parser_type!="N":
        gene_parser_type =  str(input("Run full Program? Y/N"))
    main()
#(((((((YLR293C_CLIB382:0.000001,(YLR293C_FostersO:0.001010,YLR293C_FostersB:0.000001):0.000001):0.000001,YLR293C_VL3:0.000001):0.000001,YLR293C_CBS7960:0.000001):0.000001,YLR293C_AWRI1631:0.000001):0.000001,YLR293C_Kyokai7:0.001907):0.000001,YLR293C_M22:0.000953):0.000953,YLR293C_LalvinQA23:0.000001,(YLR293C_EC1118:0.000001,(YLR297W_BC187:0.000001,YLR293C_AWRI796:0.001010):0.000001):0.000001);
