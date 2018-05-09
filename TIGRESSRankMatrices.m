function [AUROC_score_orig,AUROC_score_area, TPR_array_area, FPR_array_area,TPR_array_orig, FPR_array_orig]=...
    TIGRESSRankMatrices(A_app_ind, A, L,R)
% TIGRESSRankMatrices computes the AUROC of the list of R infered matrices, A_app_ind (G*G*R)
% compared the A matrix (G*G) with a threshold L_thr. This threshold is
% used instead of L in the canonical TIGRESS. TIGRESSRankMatrices outputs
% the two AUROCs and the list of couples (TPR,FPR) in case one wants to plot
% the ROC for the two methods of TIGRESS (original and area).

% Compute the TIGRESS socres
Score_mat_area=sum(A_app_ind.*(L+1-A_app_ind),4)/R/L;   
Score_mat_orig=sum(logical(A_app_ind>0),4)/R;        

% Sort the scores and store in I their ranks with respect to the whole 
% matrix
[~,I_area]=sort(Score_mat_area(:),'descend');
[~,I_area]=sort(I_area);
I_area=reshape(I_area,size(Score_mat_area));

[~,I_orig]=sort(Score_mat_orig(:),'descend');
[~,I_orig]=sort(I_orig);
I_orig=reshape(I_orig,size(Score_mat_orig));

Rnk_array_TIGRESS_area=zeros(size(A_app_ind,1),size(A_app_ind,2),numel(A));
Rnk_array_TIGRESS_orig=Rnk_array_TIGRESS_area;

TPR_array_orig=zeros(1,length(numel(A)+1));
FPR_array_orig=zeros(1,length(numel(A)+1));

TPR_array_area=zeros(1,length(numel(A)+1));
FPR_array_area=zeros(1,length(numel(A)+1));

listnnz_A=1:size(A,1);
A_wo_diag = A - diag(diag(A));
listnz_A=intersect(find(~logical(sum(A_wo_diag,1))),find(~logical(sum(A_wo_diag,2))));
listnnz_A = setdiff(listnnz_A,listnz_A);
A_wo_diag_red=A_wo_diag(listnnz_A,listnnz_A);

% Threshold the ranks based on k, compute the true and flase positive rates
% through CompRoc
for k=1:(numel(A)+1)
   A_app_Rnk=(I_area<k);
   Rnk_array_TIGRESS_area(:,:,k)=A_app_Rnk;
   [TPR_array_area(k), FPR_array_area(k)]=CompROC(A_app_Rnk,A_wo_diag_red,listnnz_A);   
 
   A_app_Rnk=(I_orig<k);
   Rnk_array_TIGRESS_orig(:,:,k)=A_app_Rnk;
   [TPR_array_orig(k), FPR_array_orig(k)]=CompROC(A_app_Rnk,A_wo_diag_red,listnnz_A);
end

AUROC_score_orig=trapz(FPR_array_orig(1:end),TPR_array_orig(1:end));
AUROC_score_area=trapz(FPR_array_area(1:end),TPR_array_area(1:end));
end