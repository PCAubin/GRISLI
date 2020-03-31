function [TPR, FPR, PPV]=CompROC(A_app_Rnk,A_wo_diag_red,listnnz_A)
% CompROC computes the true positive ratio (TPR) and false positive ratio 
% (FPR) of A_app_Rnk (G*G) with respect to A_wo_diag_red which is 
% the submatrix of the binary matrix A (G*G) where only
% the non-diagonal elements were kept and such that for every i between 1
% and G there is at least a 1 in either row i or column i of A_wo_diag_red 
% (meaning that the transcription factor i has to have an edge in A).
% It takes as inputs a rank matrix A_app_Rnk(GxG), a submatrix A_wo_diag_red of A
% and the list listnnz_A of the active transcription factors of A (the ones with at least
% one nonzero value).

% Pierre-Cyril Aubin-Frankowski, 2018


% Compute the submatrix of A_app_Rnk based on the submatrix A_wo_diag_red
% of A and the rules (non-diagonal, non empty row and column in A)
A_app_wo_diag = A_app_Rnk - diag(diag(A_app_Rnk));
A_app_wo_diag_red=A_app_wo_diag(listnnz_A,listnnz_A);

P=nnz(A_wo_diag_red);
N=size(A_wo_diag_red,1)^2-size(A_wo_diag_red,1)-P;

TN=sum(sum((A_app_wo_diag_red==0)&(A_wo_diag_red==0)))-size(A_app_wo_diag_red,1);
TP=sum(sum((A_app_wo_diag_red~=0)&(A_wo_diag_red~=0)));

FPR=1-TN/N;
if N==0
    FPR=1;
end
TPR=TP/P;
if P==0
    TPR=0;
end
FP=N-TN;
PPV=TP/(TP+FP);
if TP+FP==0
    PPV=1;
end
end