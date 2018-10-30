function [ROC_score_orig_mat,ROC_score_area_mat, TPR_array_area, FPR_array_area,TPR_array_orig, FPR_array_orig]=...
    ScoreMat_TIGRESS(A_app_array_ind, A,L_array,R)
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

for i=1:length(L_array)
    A_app_ind=double(A_app_array_ind(:,:,i,:));
    [AUROC_score_orig,AUROC_score_area, TPR_array_area, FPR_array_area,TPR_array_orig, FPR_array_orig ]=...
    TIGRESSRankMatrices(A_app_ind, A, L_array(i),R); %       
    ROC_score_orig_mat(i)=AUROC_score_orig;
    ROC_score_area_mat(i)=AUROC_score_area;
end
end