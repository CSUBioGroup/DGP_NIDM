% % function [TableScores ] = A_NIDM(AdjSet, P0, aaa, aaa, aaa, aaa, aaa, aaa, IndividualNetRes )
function [TableScores ] = A_NIDM(AdjSet, P0, MDL, IndividualNetRes )
    % % % % % % % % % % % % % % % % % % % %  
    % NormalizationType ='LaplacianNormalizationMeanDegree', 'ProbabilityNormalizationColumn', 'LaplacianNormalization'     
    % %     if length( ID_dis )>1; error('The number of disease ids in ID_dis is larger than 1'); end 
	if ~exist('MDL','var') || isempty(MDL); MDL = 'M' ; end  
	if ~exist('IndividualNetRes','var') || isempty(IndividualNetRes);  IndividualNetRes = false; end  
    tau_vec                     = [];
    isdebug                     = false;   
    ImpulseRestart_Times        = [ 1000 ,  1    ];   % period and times     
	weight_inter_net            = 1;    %             
    DiffusionConstant           = 1 ; 
	delta_t_input               = 0.005;   
    IndividualNormalizationType = 'None' ;  %   Normalize each sub-network?
    %
    switch MDL
        case 'M'; NormalizationType='LaplacianNormalizationMeanDegree';  
        case 'C'; NormalizationType='ProbabilityNormalizationColumn';   
        case 'R'; NormalizationType='ProbabilityNormalizationRow';      
        case 'L'; NormalizationType='LaplacianNormalization';         
        otherwise; error('There is no definition.');
    end 
	% 控制脉冲时间周期；大于1 表示多久出现一次脉冲；若小于1 表示 脉冲出现的时间比，如0.1表示1/10每10个时间步出现一次脉冲；若等于1，表示持续出现  
    if ImpulseRestart_Times(1)>1
        ImpulseRestart  = ceil( ImpulseRestart_Times(1) ); 
    elseif 0<ImpulseRestart_Times(1) && ImpulseRestart_Times(1)<1 
        ImpulseRestart  = max(1, round( 1/ImpulseRestart_Times(1)   )  );     
    else
        error('restart is wrong');        
    end  
    ImpulseTimes   = max(1, ceil( ImpulseRestart_Times(2) ) );  % 保证至少计算一个周期 
    N_max_time_seq = ImpulseRestart*ImpulseTimes ;
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
    % %     disp(['Impulse type:',NormalizationType])       
    % IsEachAdjColNormalized = false;    % % 文章中原始方法 未对每个网络标准化   
    % IsEachAdjColNormalized = true ;    % %  对每个网络标准化 ，从概率上更好理解  
 
    InterLayerMappingType   = 'FullMapping'; 
    SetIsolatedNodeSelfLoop = true;   
    [A_nLxnL,N_node, L_net] = getMultiplexMatrixFromAdjSet_IN(AdjSet, weight_inter_net,  IndividualNormalizationType, InterLayerMappingType,  SetIsolatedNodeSelfLoop, DiffusionConstant) ; 
    A_nLxnL = getNormalizedMatrix_IN( A_nLxnL , NormalizationType, true  );
    A_nLxnL = diag(sum(A_nLxnL,2))-A_nLxnL ; 
 
    %     if ~exist('tau_vec','var') || isempty(tau_vec)
    %         tau_vec = repmat(1/L_net,L_net,1); 
    %     end         
    if ~issparse(  A_nLxnL )
        A_nLxnL = sparse( A_nLxnL );
    end
    % %  
    if isempty(delta_t_input)
        delta_t = 0.005  ; %    0.01
    else
        delta_t = delta_t_input; %  
    end
    % 
	% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   
    beta =  1  ; 
    %     delta_t =0.005; 
    %     ImpulseRestart = 1000;   
    Ut     = [ 1; 0 ]*1; 
    P0     = reshape(P0,[],1)   ; 
    B      = [ P0, ~P0 ] ; 
    Bm     = repmat(B, L_net, 1 );  
    BmUt   = ( (1/delta_t)*Bm*Ut )  ;  
    Xt_seq = zeros(N_node*L_net,N_max_time_seq ) ;     
    X0     =  ( ones( size(BmUt) ) *0 ); 
    Xt1    = X0 ; 
    for T = 1: N_max_time_seq
        if mod( T, ImpulseRestart ) ==1 
            Xt2 = Xt1 + (  - A_nLxnL*Xt1  - beta*( Xt1-X0    ) + (  BmUt   ) ).*delta_t  ; % impulsive signal 
        else 
            Xt2 = Xt1 + ( - A_nLxnL*Xt1 - beta*( Xt1-X0  ) ).*delta_t ;
        end
        Xt1           = Xt2 ;  
        Xt_seq(:, T ) = Xt2;  
        % check
        if isdebug && ( any(   isnan( Xt2(:)  )   ) )  || any(   isinf( Xt2(:)  )   )
            disp( [ 'There is nan or inf in Step: ', num2str( T ) ] ); 
            disp( full( sum(A_nLxnL(:)    )     )   ); 
            disp( full(    max(Xt2)     )   ); 
            disp( full(    min(Xt2)     )   ); 
            if any(   isinf( Xt2(:)  )   ) || any(   isnan( Xt2(:)  )   )
                find(   Xt2( isinf( Xt2(:)  ) )     )
                find(   Xt2( isnan( Xt2(:)  ) )     )
                error('There is nan or inf : Impulse dynamics. ' );  
            end
        end
    end
 
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
    i_seq = 1  ; 
    Xt_seq_final        = Xt_seq(:, [N_max_time_seq-(i_seq)*ImpulseRestart+1: N_max_time_seq-(i_seq-1)*ImpulseRestart ]     );  % exract final circle   
    Xt_seq_final_NxIRxL = zeros(N_node,ImpulseRestart,  L_net );  % put into 3D array  
    for ii = 1: L_net
        ids                         = (ii-1)*N_node  + [1:N_node];   %
        Xt_seq_final_NxIRxL(:,:,ii) = Xt_seq_final(ids,:);
    end
    % 
    % %     EffSignatures = {'Xmax','Xmean','Vmt', 'cX5'   } ;  
	[Res_allnets]           = getTimeSeqFeaturesInOneCycle_IN( Xt_seq_final_NxIRxL  ) ;  
	
    TableScores             = table; 
    AdjSetfieldnameset      = fieldnames( AdjSet); 
    Resfieldnameset_allnets = fieldnames( Res_allnets);  
    % RankGeoMean
    for jj = 1: length( Resfieldnameset_allnets )  
        Scores = Res_allnets.(Resfieldnameset_allnets{jj})(:,1,:) ; 
        Ranks  = getRankingOfScoreList_IN(   Scores, 'descend' ) ;    
        % TableScores.([Resfieldnameset_allnets{jj},'_','rgm'] ) = -geomean( Ranks , 3  ) ;    
        TableScores.([Resfieldnameset_allnets{jj}] ) = -geomean( Ranks , 3  ) ;    
    end   
% % %     %
% % %     if IndividualNetRes
% % %         for ii_L_net = 1: L_net
% % %             % ids = (ii_L_net-1)*N_node  + [1:N_node];   %
% % %             for jj = 1: length( Resfieldnameset_allnets )
% % %                 TableScores.([AdjSetfieldnameset{ii_L_net},'_',Resfieldnameset_allnets{jj}] ) = Res_allnets.(Resfieldnameset_allnets{jj})(:,1,ii_L_net) ;   
% % %             end 
% % %         end 
% % %     end  
    
    TableScores.Properties.VariableNames =  strcat(['NIDM_',MDL,'_'] , TableScores.Properties.VariableNames ); 
    
    if isempty(TableScores)
        error('TableScores is empty') 	
    end	
	    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [Res] = getTimeSeqFeaturesInOneCycle_IN( Xt_seq , minusFirst  )
% % Xt_seq 每行对应一个节点的时序值；每一列对应一个时间点上所有节点的状态值; 第三维对应不同网络，每个网络的结果对应一层
    if ~exist('minusFirst','var') || isempty( minusFirst )
        % %         minusFirst = true;  % 是否计算增量  
        minusFirst = false;  % 是否计算增量  
    end
    % %     minusFirst = 0;
    if minusFirst
        Xt_seq = Xt_seq-repmat(Xt_seq(:,1,:), 1, size(Xt_seq,2)  ); 
    end
    [Xt_max, T_max ]= max(Xt_seq, [], 2);  
    Xt_mean = mean(Xt_seq, 2);  
    %
    Res.Xmx = Xt_max;    %峰值
    Res.Xmn = Xt_mean;                
    Res.Vmt = Xt_max./(T_max +eps );  
    Res.Xmt = sum(Xt_seq(:,1:5,:),2 )   ;   % default cX5  
    
    % %     Res.Tinv = 1./(T_max +eps );   %  峰值 时间 的倒数 
    % % %     Res.Xmin = Xt_min;     %% Xt_min(1:5,:,:)            % 

%     j_dX = 0; 
%     Res.(['X',num2str(j_dX +1)]) = Xt_seq(:,j_dX+1,:)   ;  j_dX = j_dX +1 ; 
%     Res.(['X',num2str(j_dX +1)]) = Xt_seq(:,j_dX+1,:)   ;  j_dX = j_dX +1 ; 
%     Res.(['X',num2str(j_dX +1)]) = Xt_seq(:,j_dX+1,:)   ;  j_dX = j_dX +1 ; 
%     Res.(['X',num2str(j_dX +1)]) = Xt_seq(:,j_dX+1,:)   ;  j_dX = j_dX +1 ; 
%     Res.(['X',num2str(j_dX +1)]) = Xt_seq(:,j_dX+1,:)   ;  j_dX = j_dX +1 ; 
%     Res.(['X',num2str(j_dX +1)]) = Xt_seq(:,j_dX+1,:)   ;  j_dX = j_dX +1 ; 
%     Res.(['X',num2str(j_dX +1)]) = Xt_seq(:,j_dX+1,:)   ;  j_dX = j_dX +1 ;  
% % % % % % %     %
%     j_cumdX = 0;
%     Res.(['cX',num2str(j_cumdX +1)]) = sum(Xt_seq(:,1:j_cumdX+1,:),2 )   ;  j_cumdX = j_cumdX +1 ; 
%     Res.(['cX',num2str(j_cumdX +1)]) = sum(Xt_seq(:,1:j_cumdX+1,:),2 )   ;  j_cumdX = j_cumdX +1 ; 
%     Res.(['cX',num2str(j_cumdX +1)]) = sum(Xt_seq(:,1:j_cumdX+1,:),2 )   ;  j_cumdX = j_cumdX +1 ; 
%     Res.(['cX',num2str(j_cumdX +1)]) = sum(Xt_seq(:,1:j_cumdX+1,:),2 )   ;  j_cumdX = j_cumdX +1 ; 
%     Res.(['cX',num2str(j_cumdX +1)]) = sum(Xt_seq(:,1:j_cumdX+1,:),2 )   ;  j_cumdX = j_cumdX +1 ; 
%     Res.(['cX',num2str(j_cumdX +1)]) = sum(Xt_seq(:,1:j_cumdX+1,:),2 )   ;  j_cumdX = j_cumdX +1 ; 
%     Res.(['cX',num2str(j_cumdX +1)]) = sum(Xt_seq(:,1:j_cumdX+1,:),2 )   ;  j_cumdX = j_cumdX +1 ;   
    
% %     Res.base = rand( size(Xt_max) );  
%     Res.base = ones( size(Xt_max) );   %  性能评价函数将对其完全随机化 得到随机化结果 
% % % %     Res.base2 = ones( size(Xt_max) );   
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [Ranks, ord ]= getRankingOfScoreList_IN(   ScoreList, sorttype, OP_IdenticalValues )
% by Xiang 
% 2019-2 
    if ~exist('OP_IdenticalValues','var') || isempty(OP_IdenticalValues)
       OP_IdenticalValues = 'MeanRank';  % % for elements with the same values 
    end
    % 
    [~, ord ] = sort(ScoreList , 1,  sorttype );  
    
    IsSparse_ScoreList = issparse( ScoreList );
    if IsSparse_ScoreList   
        ScoreList = full( ScoreList ); 
    end
    
    IDlist = [1: size( ScoreList, 1 ) ]'; 
    Ranks = zeros( size( ScoreList ) ); 
    rank_t = zeros( size(IDlist) ); 
    for d2=1:size( ScoreList ,2 )
        for d3 = 1:size( ScoreList ,3)
            for d4 = 1:size( ScoreList ,4)
                score_t    = ScoreList(:,d2,d3,d4); 
                rank_t(ord(:,d2,d3,d4)) = IDlist; 
                % Ranks( ord(:,d2,d3,d4),   d2,d3,d4) = IDlist ; 
                % 
                if ~strcmpi(OP_IdenticalValues, 'None')   % 对于相同元素值，采取操作 
% %                 if strcmpi(OP_IdenticalValues, 'MeanRank') || strcmpi(OP_IdenticalValues, 'RandPermutation')  % 对于相同元素值，采用平均 rank 
                    [uniqueScores, ~,ic] = unique( score_t ) ;
                    if length( ic  ) ~= length( uniqueScores  )
                        for ii_uniqueScorese = 1: length( uniqueScores  ) 
                            idx = ( ic== ii_uniqueScorese ) ; 
                            n_thisscore = nnz( idx )   ;
                            if n_thisscore>1
                                if strcmpi(OP_IdenticalValues, 'MeanRank') 
                                    % Ranks( idx,   d2,d3,d4 ) = mean(Ranks(idx,   d2,d3,d4), 1);  %对于相同分值的元素赋予相同的所有ranks的均值
                                    rank_t( idx ) = mean(rank_t(idx), 1);  %对于相同分值的元素赋予相同的所有ranks的均值
                                    % % sum( labels_ord( idx ) )/nnz( idx )
                                elseif strcmpi(OP_IdenticalValues, 'RandPermutation') 
                                    ind = find(idx); 
                                    ind_randperm = ind(  randperm( n_thisscore  )  ); 
                                    rank_t( ind ) = rank_t( ind_randperm ); 
                                else
                                    error('There is no definition of OP_IdenticalValues');
                                end
                            end  
                            % %sum(labels_ord) 
                        end
                    end
                end
                Ranks( : ,   d2,d3,d4) = rank_t;  
                % % 
            end
        end
    end 
    
    if IsSparse_ScoreList   % 
        Ranks = sparse( Ranks );  
    end

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [A_nLxnL,N_node, L_net] = getMultiplexMatrixFromAdjSet_IN(AdjSet, pro_jump, NormalizationType, InterLayerMappingType, SetIsolatedNodeSelfLoop, DiffusionConstant)
    if ~exist('NormalizationType','var') 
        NormalizationType = 'col';
    elseif isempty(NormalizationType) || strcmpi(NormalizationType,'none') 
        NormalizationType = 'None';
    end
    %
    if ~exist('InterLayerMappingType','var') 
        InterLayerMappingType ='FullMapping';
    % % else
    % %     InterLayerMappingType ='PartialMapping'; % 
    end
    %
    if ~exist('SetIsolatedNodeSelfLoop','var') 
        SetIsolatedNodeSelfLoop = true;  
    end
    %
    if ~exist('DiffusionConstant','var') || isempty(DiffusionConstant)
        DiffusionConstant = 1-pro_jump; 
    end

    delta =pro_jump;
    %
    if isa(AdjSet,'struct') || isa(AdjSet,'table')
        if isa(AdjSet,'struct')
            fieldnameset = fieldnames( AdjSet); 
        elseif isa(AdjSet,'table')
            fieldnameset = AdjSet.Properties.VariableNames; 
        end
        L_net  =  length( fieldnameset ); 
        N_node =  length( AdjSet.(fieldnameset{1}));
        NxL = N_node*L_net ; 
        %
        if L_net==1;  delta = 0; end %% no jumping for only single layer. 
        %         
        A_nLxnL = repmat( (delta/(L_net-1+eps)).*speye( N_node, N_node), L_net,L_net);            
        for ii_net = 1: L_net
            idx = N_node*(ii_net-1)+[1: N_node ] ; 
            if strcmpi(NormalizationType,'None') 
                A_nLxnL(idx,idx)=  (DiffusionConstant).*AdjSet.(fieldnameset{ii_net}); 
            else
                A_nLxnL(idx,idx)=  (DiffusionConstant).*getNormalizedMatrix_IN( AdjSet.(fieldnameset{ii_net}) , NormalizationType, SetIsolatedNodeSelfLoop ); 
            end                   
        end  
        
    else
        error(['AdjSet is wrong. It should be a cell matrix or struct.' ]);
        
    end 
        
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function WAdj = getNormalizedMatrix_IN(Adj, NormalizationType, SetIsolatedNodeSelfLoop  ,IsMeanDegreeForExistNodeSet )
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% Adj  adjecent matrix
% % NormalizationType: 
% % 'probability normalization'  
% % 'laplacian normalization' 
% SetIsolatedNodeSelfLoop    set isolated node
% >= Matlab 2016
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
%     if ~issparse(Adj)
%         Adj = sparse( Adj );
%     end   
    if ~exist('IsMeanDegreeForExistNodeSet','var') || isempty(IsMeanDegreeForExistNodeSet)
        IsMeanDegreeForExistNodeSet = false;
    end
 

    if ischar(NormalizationType)
    %         NormalizationType =  (NormalizationType);
        switch  lower( NormalizationType )
            case lower( { 'column','col',  ...
                    'ProbabilityNormalizationColumn','ProbabilityNormalizationCol',...
                    'ProbabilityColumnNormalization','ProbabilityColNormalization',...
                    'NormalizationColumn','NormalizationCol' , ...
                    'ColumnNormalization','ColNormalization'   })
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =1;
            case lower({ 'row' ,'ProbabilityNormalizationRow' ,'NormalizationRow' ,'ProbabilityRowNormalization' ,'RowNormalization'   })
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =2;
            case lower('LaplacianNormalization')
                NormalizationName = NormalizationType; 
            case lower('LaplacianNormalizationMeanDegree')
                NormalizationName = NormalizationType; 
            case lower('ColNorm2')
                NormalizationName = NormalizationType; 
            case lower('RowNorm2')
                NormalizationName = NormalizationType; 
            case lower({'none', 'None', 'NONE'})
                % NormalizationName = 'None'; 
                WAdj = Adj; 
                return; 
            otherwise
                error(['There is no type of normalization: ',char( string(NormalizationType) )] );
        end
        
    elseif isnumeric(  NormalizationType   ) 
        NormalizationName =  ( 'ProbabilityNormalization' ) ;  %  'Random Walk'  
        dim = NormalizationType; 
        
    elseif isempty( NormalizationType )
        WAdj = Adj; 
        return;  
        
    else; error('There is no defintion of NormalizationType')
    end 
    % NormalizationName = lower( NormalizationName );
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     %
    matlabveryear = textscan(version('-release'),'%d') ; matlabveryear=matlabveryear{1};
    if matlabveryear>=2016 
        switch lower( NormalizationName )
            case lower( 'ProbabilityNormalization' )
                degrees = sum(Adj,dim);
                if any( degrees~=1)
                    WAdj = Adj./ ( degrees+eps  );           
                    % % WAdj = Adj./ repmat( degrees +eps,[size(Adj,1),1]); 
                else
                    WAdj = Adj; 
                end
                % 
                if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2) 
                    ii = find( ~degrees ); 
                    idx = sub2ind( size(Adj), ii,ii ); 
                    WAdj(idx) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( 'LaplacianNormalization')
                deg_rowvec = ( sum(Adj,1) ).^0.5;  
                deg_colvec = ( sum(Adj,2) ).^0.5;   
                WAdj = (Adj./(deg_colvec+eps))./(deg_rowvec+eps) ;    
                % 
                if SetIsolatedNodeSelfLoop && size(Adj,1)==size(Adj,2)
                    ii = find( ~sum(Adj,2) ) ; 
                    % size(  WAdj )
                    % size(  Adj )
                    WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( 'LaplacianNormalizationMeanDegree')
    %             n_node = length( Adj ); 
                k_col = sum( Adj, 2 ); 
                k_row = sum( Adj, 1 );             
                if IsMeanDegreeForExistNodeSet 
                    km1 = sum(k_col)./ (nnz(k_col)+eps);  
                    km2 = sum(k_row)./ (nnz(k_row)+eps);                  
                else
                    km1 = sum(k_col)./ length(k_col);  
                    km2 = sum(k_row)./ length(k_row);                  
                end
                WAdj = Adj./( (km1.^0.5)*(km2.^0.5)  +eps) ;    
                % 
                if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2)
                    ii = find( ~sum(Adj,2) ); 
                    WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( {'ColNorm2'} )   
                WAdj = Adj./ ( sqrt(sum( Adj.^2 ,1 )) +eps ); 

            case lower( {'RowNorm2'} )    
                WAdj = Adj./ ( sqrt(sum( Adj.^2 ,2 )) +eps ); 

            case lower( {'None','none'} )
                WAdj = Adj;   % 不做任何处理  
            otherwise
                error(['NormalizationName is wrong: ',char(string(NormalizationName) )   ]);
        end

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
    else
    %     error('errorerrorerrorerrorerrorerrorerrorerror');
        switch lower( NormalizationName )
            case lower( 'ProbabilityNormalization' )
                degrees = sum(Adj,dim);
                if any( degrees~=1)
                    % WAdj = Adj./ ( degrees+eps  ); 
                    WAdj = getMatrixOperation(Adj, degrees, './') ;
                     % % WAdj = Adj./ repmat( degrees +eps,[size(Adj,1),1]); 
                else
                    WAdj = Adj; 
                end
                % 
                if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2) 
                    ii = find( ~degrees ); 
                    idx = sub2ind( size(Adj), ii,ii ); 
                    WAdj(idx) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( 'LaplacianNormalization')
                deg_rowvec = ( sum(Adj,1) ).^0.5;  
                deg_colvec = ( sum(Adj,2) ).^0.5;   
                % %             WAdj = (Adj./(deg_colvec+eps))./(deg_rowvec+eps) ; 
                WAdj = getMatrixOperation(Adj,  deg_colvec, './') ;
                WAdj = getMatrixOperation(WAdj, deg_rowvec, './') ;
                % 
                if SetIsolatedNodeSelfLoop && size(Adj,1)==size(Adj,2)
                    ii = find( ~sum(Adj,2) ) ; 
                    % size(  WAdj )
                    % size(  Adj )
                    WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( 'LaplacianNormalizationMeanDegree')
                n_node = length( Adj );
                km = sum( Adj(:) )./ n_node;  
                WAdj = Adj./( (km.^0.5)*(km.^0.5)  +eps) ;    
                % 
                if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2)
                    ii = find( ~sum(Adj,2) ); 
                    WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( {'ColNorm2'} )   
                % % WAdj = Adj./ ( sqrt(sum( Adj.^2 ,1 )) +eps ); 
                WAdj = getMatrixOperation(Adj,  sqrt(sum( Adj.^2 ,1 )), './') ;

            case lower( {'RowNorm2'} )    
                % % WAdj = Adj./ ( sqrt(sum( Adj.^2 ,2 )) +eps ); 
                WAdj = getMatrixOperation(Adj,  sqrt(sum( Adj.^2 ,2 )), './') ;

            case lower( {'None','none'} )
                WAdj = Adj;   % 不做任何处理  
            otherwise
                error(['NormalizationName is wrong: ',char(string(NormalizationName) )   ]);
        end 

    end
 
end



