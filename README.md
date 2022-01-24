# NIDM: network impulsive dynamics on multiplex biological network for disease-gene prediction.
The prediction of genes related to diseases is important to the study of the diseases due to high cost and time consumption of biological experiments. Network propagation is a popular strategy for disease-gene prediction. However, existing methods focus on the stable solution of dynamics while ignoring the useful information hidden in the dynamical process, and it is still a challenge to make use of multiple types of physical/functional relationships between proteins/genes to effectively predict disease-related genes. Therefore, we proposed a framework of network impulsive dynamics on multiplex biological network (NIDM) to predict disease-related genes, along with four variants of NIDM models and four kinds of impulsive dynamical signatures (IDSs). NIDM is to identify disease-related genes by mining the dynamical responses of nodes to impulsive signals being exerted at specific nodes. By a series of experimental evaluations in various types of biological networks, we confirmed the advantage of multiplex network and the important roles of functional associations in disease-gene prediction, demonstrated superior performance of NIDM compared with four types of network-based algorithms and then gave the effective recommendations of NIDM models and IDS signatures. To facilitate the prioritization and analysis of (candidate) genes associated to specific diseases, we developed a user-friendly web server, which provides three kinds of filtering patterns for genes, network visualization, enrichment analysis and a wealth of external links (http://bioinformatics.csu.edu.cn/DGP/NID.jsp). NIDM is a protocol for disease-gene prediction integrating different types of biological networks, which may become a very useful computational tool for the study of disease-related genes.


## Requirements
Matlab 2016 or above.   


## Codes 
#runCV_NIDM.m: cross-validation code.  <br>
This code allows parallel execution. You can change "parfor" to "for" to cancel parallel execution  <br>
 
#A_NIDM.m: a version of NIDM in the study. <br>   
[TableScores ] = A_NIDM(AdjSet, P0, MDL )<br>
% Input:  <br>
% AdjSet: set of gene-gene matrices <br> 
% P0: a vector marking excitation nodes of signals  <br>
% MDL: M (default), C, R or L <br>
% Ouput: <br>
% TableScores: a table whos variable record the scores of genes.  <br> 


## Dataset
Data is located in the directory: ./data <br>
./demoDataset.mat, including disease-gene associations (DG), and multiple types of gene-gene associations (GG);  <br> 


## Results 
The results will be automatically saved into the directory: results.  

## cite
If you use NIDM in your research, please cite: <br> 
Ju Xiang, Jiashuai Zhang, Ruiqing Zheng, Xingyi Li, Min Li. NIDM: network impulsive dynamics on multiplex biological network for disease-gene prediction, Briefings in Bioinformatics, Volume 22, Issue 5, September 2021, bbab080, https://doi.org/10.1093/bib/bbab080
.


## contact<br>
Email: xiang.ju@foxmail.com 
