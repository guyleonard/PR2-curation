<img align="right" width="80" alt="logo (1)" src="https://user-images.githubusercontent.com/71433342/188953119-c830769a-4c40-47c5-9ca6-f62dd45804cb.png">

# Final Master Thesis - PR2 curation: a shiny app version of the EukRef data curation pipeline

EukRef is an initiative that allows users interested in a particular group of microbial eukaryotes to retrieve all sequences belonging to that group, place those sequences in a phylogenetic tree, and curate taxonomic and environmental information for the group. EukRef is currently part of the Protist Ribosomal Reference (PR2) database and its curation workflow has been incorporated to the database. Here we present a web application developed in the R environment that adapts the above pipeline to make the user experience friendlier and simpler. PR2 curation allows:

(i) downloading data from PR2,

(ii) building the initial tree,

(iii) building the tree again applying the necessary modifications,

(iv) filtering the initial fasta file based on the edited tree.

## 1. Installation

### 1.1 Software needed

Tools necessary to install to run the pipeline are listed below. You should install each of these programs on your computer.

-   [Vsearch](https://github.com/torognes/vsearch): a 64-bit, open source, free, multithreaded tool for processing and preparing metagenomics, genomics, and population genomics nucleotide sequence data. VSEARCH uses a global optimal aligner (Needleman-Wunsch full dynamic programming).

-   [MAFFT](https://mafft.cbrc.jp/alignment/server/): is a high-speed multiple sequence alignment program for unix-like operating systems.

-   [trimAl](http://trimal.cgenomics.org/): tool for automated removal of spurious sequences or misaligned regions from a multiple sequence alignment. It implements a series of automated algorithms that trim the alignment by finding optimal thresholds based on inherent features of the input alignment, which will be used to increase the signal-to-noise ratio after the alignment trimming phase.
-   [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/): is a program for sequential and parallel Maximum Likelihood based inference of large phylogenetic trees. It can also be used for post analyses of sets of phylogenetic trees, analyses of alignments, and evolutionary placement of short reads.

-   [R](https://cran.r-project.org/bin/windows/base/): is a language and environment for statistical computing and graphics.
-   [RStudio](https://www.rstudio.com/): an integrated development environment for R, with a console, syntax-highlighting editor that supports direct code execution, and tools for plotting, history, debugging and workspace management.

You can install the latest development version of the code using the _shiny_ R package.

```{r, eval=FALSE}
## Install devtools, if you haven't already.
# install.packages("shiny")
library(shiny)

runGitHub(repo = "PR2-curation", username = "delCampoLab", ref = "main")
```

## 2. Pipeline improvements
We present the standardized guidelines and open source operational tools developed by this application. The final results are a phylogenetic reference and alignment tree and a selected reference database with selected accession numbers and taxonomic string. To enable efficiency and consistency, the pipeline was developed to select and annotate diverse eukaryotic lineages to comprehensively capture their existing sequence diversity (Figure 1b). Curation begins with a broadly sampled alignment and corresponding 18S rRNA phylogeny. The initial set of sequences is downloaded from PR2 and becomes the input for the PR2 curation workflow. The set of input sequences and relevant outgroups are aligned using MAFFT. It is then automatically trimmed using trimAl and used for phylogenetic inference with RAxML, which can easily handle typical datasets of hundreds to several thousand sequences. The resulting tree is the starting point for sequence curation and classification. Next, curators manually examine the tree to identify discrepancies, such as long branches, which may be potential artifacts or chimeras that escaped the initial filtering, which must be removed. 
PR2 curation has the interface of an application and the user can do all the manual curation with clicks. After the removal of this problematic data, a new alignment and a tree with the remaining sequences are constructed. This information, together with the tree, is the starting point for the classification of the group and for each sequence. These results are combined with previous taxonomic knowledge and improved metadata are manually incorporated throughout this process.

The reference trees serve us in addressing the goal of improving the eukaryotic taxonomic framework and creating a repository of correctly annotated high-throughput environmental sequencing (HTES) reads.

![pipeline](https://user-images.githubusercontent.com/71433342/189580353-437fd447-f4f1-4bf3-af55-dbbf44381c78.png)

## 3. Features

**a) Default features**

- When clustering sequences with VSEARCH, a similarity threshold of 97% is set. The code chooses the longest sequence as the representative sequence for each cluster (these will be used later for alignment). Sequences shorter than 500 bp are also removed.
- In trimAl all positions in the alignment with gaps in 70% or more of the sequences are removed. A threshold of 0.001 is set for the minimum allowed average similarity. 
- The algorithm selected in RAxML was a fast bootstrap analysis and search for the best scoring ML tree in a single program run.


**b) Functions modifiable by the user**

The user can make modifications of a visual nature, modifications that do not change the properties of the tree, or of an analytical nature, modifications that change the properties of the tree. All of them are defined in Table 1. 

| Function | Description |
|:-------------------------------|:---------------------------------------| 
| <em>Visual<em> | | 
| <b>Rotate<b> | To facilitate the exploration of the tree structure, there is the option to rotate the selected clade by 180 degrees indicating the node in question. In other words, for a given node, rotate swaps the position of two clades descending from this node. | 
| <b>Flip<b> | The position of the immediate descendant clades of the internal node can be exchanged with this option. Therefore, exchange positions of 2 clades that share a parent node. | 
| <em>Analytical<em> | |
| <b>Reroot<b> | Reroots a phylogenetic tree with respect to the specified outgroup or at the node specified in node. In case the tree is not rooted with the outgroup, the user can employ the use of this option. |
| <b>Remove in tree<b> | Remove branches from a phylogenetic tree due to low sequence quality, errors in sequence assembly, alignment error, phylogenetic inference error, among others. This modification will be shown in the Cladogram display. For this edition, the user must attach three files: an .xlsx file with the names of the branches to be deleted, the .cluster file where the seed sequences are identified and the FA file containing the DNA sequences. Internally, the application looks for the name of the branches provided by the user in the cluster file. These branches are deleted and, if they are seed, the branches within them are also deleted. Finally, a new phylogenetic tree is built. | 
| <b>Remove in input dataset<b> | Remove sequences from the initial FA file. The process is the same as the one described in the "Remove in tree" section but a new FA file is obtained with the resulting sequences. | 
| <b>Rename in tree<b> | Rename the branches in the tree display. For example, if the accession numbers are listed. The user must provide a txt file where in the first column the current names are specified and in a second column (separated by tabulation) the desired taxonomic studies. |
| <b>Rename in input dataset<b> | Rename the branches in the initial FA file with the DNA sequences. For example, if the accession numbers are listed. The user must provide a txt file where in the first column the current names are specified and in a second column (separated by tabulation) the desired taxonomic studies. Internally, the application searches for the name of the branches provided in the first column and replaces them with the names in the second column, providing a new FA file. |
| | | 
| <b>Save<b> | Option to download the tree creation in PDF format or, also, a new tree file to save the modifications made in the interface. |

## 4. Suessiales of a case study

To illustrate the usefulness and functioning of PR2 curation, we use Suessiales, an order of unicellular organisms of the class Dinophyceae, as a case study. 
The app consists of a left panel with all the options for modifying and editing the tree and a right panel where the phylogenetic trees are displayed, always showing the last tree in the first position (Figure 2A). 

### 1. Creation of phylogenetic tree

Two options:

#### 1.a - Contribute my tree file

The user already has a tree file. In this case only the first option has to be selected and a browser is displayed to provide the file. The phylogenetic tree is automatically displayed in the right panel.

#### 1.b - Search lineage in PR2 database

This option is used in case the user does not yet have a tree file and wants to search for a lineage in the PR2 reference database. In this case the second option is selected and three inputs are displayed. For our example we select _Order_ under '**Taxonomic category**' and write _Suessiales_ under '**Lineage group of interest**'. Internally the application searches for the DNA sequences in the PR2 reference database and generates a fasta file with the information in its working directory.

The last step is to enter the browser and input the FA file to get the tree file. This step may require some time because internally the application runs the pipeline where clustering, alignment, trimming and inference are performed. It is very important to remember that the alignment requires the presence of the outgroup sequences. It is important that the outgroup is formed by sequences that we know for sure which elements compose it (more information in the article). This file should be FA and in the same directory as the file with the PR2 sequences. Finally, the phylogenetic tree is automatically generated in the right panel.

![doble](https://user-images.githubusercontent.com/71433342/188953966-7832416a-4a3a-485d-88cb-c45e0cd3fdd8.jpg)

Next figutr shows the generated tree. The type of plot generated by the application is a Cladogram (scaled rectangular diagram) with the names of each branch on the right. Bootstrap value information is always displayed in a white box. The node number, the blue box, informs us of the value assigned to each node of the tree to refer to in case we want to make modifications. 
  
![phylo_tree (1)_page-0001](https://user-images.githubusercontent.com/71433342/189581933-a5c58a0f-dcae-40fa-b6a5-d18f6e66c50d.jpg)


### 2. Modification of the phylogenetic tree

These modifications are visual, they do not change the properties of the tree. You may want to *rotate* a node to make it easier to explore the tree structure information. To do this we must indicate the number of the node in question. In this case we rotate node 49 to have the separation of branches from the earliest to the latest.

There is also the possibility to *flip* two nodes for tree exploration. For this example, two nodes have been randomly flipped (node 87 with node 72). As in the rotation, when flipping two branches, the application highlights in two boxes the branches that are affected by this action.

### 3. Editing the phylogenetic tree

In this section the modifications are of analytical type, that is to say, they change the properties of the tree. *Reroot* allows us to reorient the tree with respect to the specified node in case it is not rooted in the outgroup. As an example, we randomly reroot at node 70. 

The user might want to *remove* branches from the tree for various reasons, for example an error in alignment, phylogenetic inference or sequence assembly. In this case we want to remove the long branch "AB693194" (Figure 6) because we suspect that it is not correct or that it is a chimera. We attach the xlsx file with the name of the mentioned branch, the cluster file and the FA file with the sequences. The process of running the pipeline internally requires a certain execution time. The output is the new reference tree and the FA file will be the new reference database. 

In some cases the accession numbers may appear as branch names. The application incorporates the option to rename. To do this, we provide a text file with a first column containing the current names, and a second column with the new ones. Remember that these columns must be separated by tabulation.

### 4. Save results

The first option is to download the Cladogram graph in PDF format. The other option is to download a tree file containing the phylogenetic tree resulting from the modifications made. Both options will always save the last printed tree in the right panel of the application.
