This file explains how to use the Project contained within to detect conservation within genes across the strains contained in this folder.

This project uses python 3.8.1 and requires numpy and matplotlib

First, use cs_final_print_req_seqs.py . This file is used to organize the genetic data from the genome files. To use this on your machine, you must edit the path names by using the input to the file directory or src_gene_files to match the  local machine.
This file also creates multiple output files.
The file works by detecting the copy of the ORF in each strain and concatenating it to a file for each gene in the gene_src. This is important for letting the phastcons parse the gene data easily. 

Second, use phastcons.py. Specify the genes you are interested testing over in the first line of main file ("for genes in []"). fill in the array with the systemic names of desired genes with a "_" appended.
The program operates in 2 main steps. Firstly, the program optimizes values for vu, and mu. It does this by first creating a binary tree and removing leaf nodes for which the gene sequence is unavailable.
Next a transversal order is generated from the tree. Then a likelihood is generated for the tree, including a matrix of the probabilities iterating along the tree length.
These steps regarding the generation of the tree and processing of it are repeated again for a tree with different branch lengths, the conserved model.
Using data from both of these trees, a viterbi algorithm is used to generate a path between the conserved and non-conserved models, using mu and vu as transition probabilities.
The number of transitions in the viterbi path is calculated. These numbers are the number of transitions from the conservative to nonconservative model, from the nonconservative to conservative model, to the conservative model, and to the nonconservative model. 
Pseudocounts are added to these number of transitions such that a floor for mu and vu is set.
By dividing the number of transitions from the conservative to nonconservative model by the number of transiitons to the nonconservative model and from the nonconservative to conservative model by the number of transiitons to the conservative model, we can calculate new values for mu and vu.
These steps are repeated until the difference between old and new values for mu and vu drops below 0.001.
Once vu and mu have been maximized, a forward backward algorithm is run on the data, using the data from the tree transition probabilites.
After the relative probabilties of the models is calculated, the probability of the conserved model is plotted.
Verticle lines are also added to this plot where each sequence ends.


phastcons_check.py represents the same implementation as phastcons but for a different genetic dataset and tree.
