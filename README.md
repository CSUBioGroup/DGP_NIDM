# NIDM
NIDM: network impulsive dynamics on multiplex biological network for disease-gene prediction.


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