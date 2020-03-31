function [ROC_score_area_mat, PR_score_area_mat, TPR_array_area, FPR_array_area,...
    PPV_array_area, A_app_array_Rnk0, ROC_score_orig_mat, PR_score_orig_mat, TPR_array_orig, FPR_array_orig,...
    PPV_array_orig]=ScoreMat_TIGRESS(A_app_array_ind, A,L_array,R, bool_orig)
% ScoreMat_TIGRESS stores the AUROC values over the parameter arrays L_array
% in ROC_score_**_mat (L_array*1).
% It also exports the true positive rate (TPR) and false positive rates 
% if we wish to plot the ROC. It takes as input the array of inferred
% matrices A_app_array_ind(G*G*L_array*R), the true binary matrix A (G*G),
% the parameter arrays, L_array and L_thr_array, and R.
% ScoreMat_TIGRESS is meant to explore the parameter array. TIGRESSRankMatrices
% suffices if there is only one L to investigate.
% Pierre-Cyril Aubin-Frankowski, 2018
ROC_score_orig_mat=zeros(length(L_array),1);  
ROC_score_area_mat=zeros(length(L_array),1);  
PR_score_orig_mat=zeros(length(L_array),1);  
PR_score_area_mat=zeros(length(L_array),1);  
A_app_array_Rnk0=repmat(A,1,1,length(L_array));

if nargin<5
    bool_orig=false;
end


for i=1:length(L_array)
    A_app_ind=double(A_app_array_ind(:,:,i,:));
%     [AUROC_score_orig,AUROC_score_area, TPR_array_area, FPR_array_area,TPR_array_orig, FPR_array_orig ]=...
%     TIGRESSRankMatrices(A_app_ind, A, L_array(i),R); % 
    [AUROC_score_area, AUPR_score_area, TPR_array_area, FPR_array_area,...
    PPV_array_area, A_app_Rnk0, AUROC_score_orig, AUPR_score_orig, TPR_array_orig, FPR_array_orig,...
    PPV_array_orig]= TIGRESSRankMatrices(A_app_ind, A, L_array(i),R, bool_orig);
    A_app_array_Rnk0(:,:,i)=A_app_Rnk0;
    ROC_score_orig_mat(i)=AUROC_score_orig;
    ROC_score_area_mat(i)=AUROC_score_area;
    PR_score_orig_mat(i)=AUPR_score_orig;
    PR_score_area_mat(i)=AUPR_score_area;
end
end