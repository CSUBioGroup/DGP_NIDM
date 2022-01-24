fdata = 'demoDataset.mat '
if ~exist( fdata,'file')
    clear all
    %     dataname = 'Dataset_HumanPPI_MencheScience2015_OMIM_HPOSim_meshDis_SymptomSim_DisGene' ;
    %PPItype   = 'SinglePPI'
    PPItype    = 'multiple' 
    fdatestrForfile = datestr(now,'yyyy.mmm.dd-HH.MM.SS')  
    netname    = 'PPI';   
    DisSetName = 'Mesh';
    % %  

    load('updateID2019_newDataset_HumanPPI_DisGene_DisSim_AllData_noGTExCoExp_weight.mat');  % if nargin<1     
    if strcmpi(netname,'PPI') 
        AdjGfG00 = sparse( data_all.DataGeneNet_1Net.matrix  ); 
        MatrixSet_gg = [];  
        if strcmpi('SinglePPI', PPItype)  
            MatrixSet_gg.Science2015PPI = AdjGfG00;    
        else
            MatrixSet_gg.matrix_binary    = data_all.DataSetPPI.matrix_binary;     
            MatrixSet_gg.matrix_complexes = data_all.DataSetPPI.matrix_complexes;     
            MatrixSet_gg.matrix_kinase    = data_all.DataSetPPI.matrix_kinase;     
            MatrixSet_gg.matrix_literature= data_all.DataSetPPI.matrix_literature;     
            MatrixSet_gg.matrix_metabolic = data_all.DataSetPPI.matrix_metabolic;     
            MatrixSet_gg.matrix_regulatory= data_all.DataSetPPI.matrix_regulatory;     
            MatrixSet_gg.matrix_signaling = data_all.DataSetPPI.matrix_signaling;            
        end 
    end
    %  
    AdjGfD            = logical(data_all.DataDiseaseNet_MeSH70DiseaseGeneSet.Matrix_Disease_Gene)';  
    Matrix_gene_dis00 = logical(data_all.DataDiseaseNet_MeSH70DiseaseGeneSet.Matrix_Disease_Gene)';
    % TableDiseases     = data_all.DataDiseaseNet_MeSH70DiseaseGeneSet.Table70MeshDiseases; 
    data_all       = [] 
    % %      
    CoEXPname = ''
    % CoEXPname  = 'CoEXP';  
    if ~isempty( CoEXPname )
        dataCoEXP  = load('Net_GETxGeneCoExp2019filteredByHumanPPI_MencheScience2015.mat'); 
        MatrixSet_gg.(['CoExp'])                    = dataCoEXP.(['GeneCoExp_pearson',num2str(20)]);                     
        dataCoEXP =[];   
    end 
    save  data\demoDataset.mat 
       
end 
% % % %
datadir    = 'data';    
ResPerfDir = 'results'; if ~exist(ResPerfDir,'dir');mkdir(ResPerfDir); end
load([datadir,filesep,'demoDataset.mat'])
 
SubNetworkNames =  fieldnames(MatrixSet_gg )';  
disp(['Net:',netname,'      DisSet:',DisSetName,'   CoEXPname:',CoEXPname,'      Date:',fdatestrForfile]);
disp( SubNetworkNames'  ) 


% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
rng(1111)    
nCVTimes  = 5;  nFold  = 5;  CVtype = [num2str(nFold),'FCV'];  MinSizeDisGeSet = nFold;     
%
n_disgenes_eachdisease = sum(Matrix_gene_dis00,1)';    
dis_IDset = find(n_disgenes_eachdisease>=MinSizeDisGeSet); 
% dis_IDset = dis_IDset(1:2)   %%%for test only %%%%%%%%%%%    
n_neighbors_artificial_linkage_interval = 99; 
% 
n_disease_in_Table   = length( dis_IDset   ); 
nCV_list             = zeros( n_disease_in_Table, 1 );   
matAUROC_nCVTimes    = cell(nCVTimes,1);
matAUPRC_nCVTimes    = cell(nCVTimes,1);
matAURecall_nCVTimes = cell(nCVTimes,1);
matAUPrec_nCVTimes   = cell(nCVTimes,1);
matRec_nCVTimes      = cell(nCVTimes,1);
matPrec_nCVTimes     = cell(nCVTimes,1);
methodset_nCVTimes   = cell(nCVTimes,1);   
date_start           = datestr(now,'yyyy.mmm.dd-HH.MM.SS');
for  i_cv = 1:nCVTimes    
    disp(['i_cv-',num2str(i_cv) ]) 
    MatrixSet_gg_copy= MatrixSet_gg; 
    %
    matAUROC    = [];
    matAUPRC    = [];
    matAURecall = []; 
    matAUPrec   = []; 
    matRec      = []; 
    matPrec     = [];    
	methodset   = [] ;
    idx_res     = 0;     
    for ii_dis = 1:n_disease_in_Table
        tic
        Matrix_gene_dis_copy = Matrix_gene_dis00 ;  
        ID_dis               = dis_IDset(ii_dis);  
        disp(['i_cv-',num2str(i_cv),'; ii_dis-',num2str(ii_dis),'; ID_dis-',num2str(ID_dis)]) 
        ac_gene_dis00 = Matrix_gene_dis_copy(:,ID_dis ); 
        idx_pos       = find( ac_gene_dis00 );  n_pos = length( idx_pos); 
        idx_neg       = find( ~ac_gene_dis00 ); n_neg = length( idx_neg); 
        n_fold_real   = min(nFold, n_pos) ;   
        ind_fold_pos  = crossvalind('Kfold', n_pos, n_fold_real ) ; 
        %ind_fold_neg  = crossvalind('Kfold', n_neg, n_fold_real ) ; 
        for i_fold = 1:n_fold_real  
            % idx_pos_train = idx_pos(ind_fold_pos~=i_fold);
            idx_pos_test    = idx_pos(ind_fold_pos==i_fold);  n_pos_test =length(idx_pos_test); 
            %   
            idx_neg_test_X100_RC      = idx_neg(   randperm(n_neg,  min(n_pos_test*n_neighbors_artificial_linkage_interval,n_neg)  ) ); n_neg_test_x100 = length(  idx_neg_test_X100_RC  ) ; % 每个测试基因选择100个对应控制基因，总共n_pos_test*100
            idx_test_pos_neg_X100_RC  = [idx_neg_test_X100_RC; idx_pos_test ] ;   
            %
            idx_neg_test_WG        = idx_neg ;  n_neg_test_all = length(  idx_neg_test_WG  ) ;               
            idx_test_pos_neg_WG    = [idx_neg_test_WG; idx_pos_test ] ;  
            %
            AdjGfD                       = Matrix_gene_dis_copy; 
            AdjGfD(idx_pos_test,ID_dis ) = 0 ;   
            P0                           = AdjGfD(:,ID_dis );
            % 
            
            %
            TableScores = table ;  
            %  M    
            [TableScores_NIDM ] = A_NIDM(MatrixSet_gg_copy, P0 ) ; 
            TableScores_NIDM.Properties.VariableNames =  strcat( TableScores_NIDM.Properties.VariableNames, [CoEXPname]  ); 
            TableScores = [TableScores,TableScores_NIDM];   
           
            % 
            [TableScores_RWRM ] = A_RWRMcp1(MatrixSet_gg_copy, P0, [],[], 'col') ; 
            TableScores = [TableScores, TableScores_RWRM];   
            [TableScores_RWRM ] = A_RWRMcp2(MatrixSet_gg_copy, P0, [],[], 'col') ; 
            TableScores = [TableScores, TableScores_RWRM];   
            
            % % % % % % % % % % % % % % % % % % % % %    
            methodset = TableScores.Properties.VariableNames ;   
            test_real = ac_gene_dis00(idx_test_pos_neg_WG);  
            [AUROCset, AUPRCset , Rec5set  , Prec5set, AURecallset, AUPrecset   ] = getPerf(test_real,TableScores(idx_test_pos_neg_WG,:) ) ; 
            %
            test_real_RC = ac_gene_dis00(idx_test_pos_neg_X100_RC);  
            [AUROCset_RC, AUPRCset_RC , Recset_RC  , Precset_RC, AURecallset_RC, AUPrecset_RC  ] = getPerf(test_real_RC,TableScores(idx_test_pos_neg_X100_RC,:) ) ; 
            idx_res = idx_res +1 ; 
            matAUROC(idx_res,:) = [AUROCset_RC, AUROCset   ];
            matAUPRC(idx_res,:) = [AUPRCset_RC, AUPRCset   ];
            matAURecall(idx_res,:)  = [AURecallset_RC,  AURecallset    ]; 
            matAUPrec(idx_res,:)    = [AUPrecset_RC,  AUPrecset  ];         
            matRec(idx_res,:)   = [Recset_RC,  Rec5set    ]; 
            matPrec(idx_res,:)  = [Precset_RC, Prec5set   ];         
            % % % % %  
        end 
    end 
    %  
	matAUROC_nCVTimes{i_cv}   = mean(matAUROC,1);
	matAUPRC_nCVTimes{i_cv}   = mean(matAUPRC,1) ;
	matAURecall_nCVTimes{i_cv}= mean(matAURecall,1); 
	matAUPrec_nCVTimes{i_cv}  = mean(matAUPrec,1); 
	matRec_nCVTimes{i_cv}     = mean(matRec,1); 
	matPrec_nCVTimes{i_cv}    = mean(matPrec,1); 
    methodset_nCVTimes{i_cv}  = methodset;    
	toc 
    disp('  ') 
    mean(matAUPRC,1)
    
end 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
%  
matAUROC_nCVTimes = cat(1,matAUROC_nCVTimes{:});    
matAUPRC_nCVTimes = cat(1,matAUPRC_nCVTimes{:});    
matAURecall_nCVTimes  = cat(1,matAURecall_nCVTimes{:}); 
matAUPrec_nCVTimes    = cat(1,matAUPrec_nCVTimes{:}); 
matRec_nCVTimes  = cat(1,matRec_nCVTimes{:}); 
matPrec_nCVTimes = cat(1,matPrec_nCVTimes{:}); 
% 
vars = [strcat( methodset_nCVTimes{1}, '_RC'),  methodset_nCVTimes{1} ] ; 
matRESmean        = [   mean(matAUROC_nCVTimes,1);... 
                        mean(matAUPRC_nCVTimes,1);...  
                        mean(matAURecall_nCVTimes,1);... 
                        mean(matAUPrec_nCVTimes,1);...  
                        mean(matRec_nCVTimes,1); ... 
                        mean(matPrec_nCVTimes,1) ] ;  
tbRESmean         = array2table(matRESmean, 'VariableNames', vars, 'RowNames',{'AUROC','AUPRC','AURecall','AUPrec','Rec','Prec'})   
% 
%  save  
date_cmplt  = datestr(now,'yyyy.mmm.dd-HH.MM.SS');
parastr     = sprintf('CVtype=%s_CVtime=%d_MSDGS%d', CVtype  ,  nCVTimes, MinSizeDisGeSet );   
outfile     = [ResPerfDir,filesep,'ResPerf_NIDM_',netname,'_',DisSetName,'_',parastr,'_',date_cmplt,'.mat'] 
save([outfile],  'tbRESmean',  'date_start',  'date_cmplt' ,  '-v7.3' )   ;   
 

 
    
    
    
    
    
    
