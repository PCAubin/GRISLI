function [AUROC_GRISLI, AUPR_GRISLI, elapsedTime, TPR_array_area, FPR_array_area,...
    PPV_array_area, A_app_array_Rnk0]=Test_GRISLI_realdata(A,X,L_array,Alpha,R,alpha_min,saveResultsBool, saveFileName,bool_orig,V,nbtry)
% Test_GRISLI_realdata computes the AUROC scores of the matrices predicted
% through the GRISLI Algorithm (Aubin-Frankowski and Vert, 2018) as compared to the 
% reference matrix A (G*G) for the given dataset found at data*filenumber*/.
% If the boolean saveResultsBool is true it can 
% export the results in the txt file saveFileName under the form:
% ['#cells','R','L','alpha','AUROC','time']
% For instance GRISLI parameters could be:
% R=1500;
% L_array=20:10:90;
% alpha_min=.3;
% saveResultsBool=false;
% Pierre-Cyril Aubin-Frankowski, 2018
if nargin<11
    nbtry=1;
end
if nargin<10
    K = SpacetimeKernel(X,Alpha);
    Knorm = KernelNormalization(K);
    V=VelocityInference(X,Knorm);
%     V0=V;
%     save('V_new_G_mat.mat','V0')
end
if nargin<9
    bool_orig=false;
end
if nargin<8
    saveFileName='GRISLI_results.txt';
end

if nargin<7
    saveResultsBool=false;
end

numcells=size(X,1);
lambda_array=[0.0001];%In case we wish to use an array of lambda values over which we compute the lasso, rather than using L


for count=1:nbtry
tic
%%4-D array to store the ranks of the inferred A of the R tests for each L in L_array
A_app_array_ind=uint32(zeros(size(A,1),size(A,2),length(L_array),R));
size_batch=floor(numcells/2);
%%Compute frequencies
for j=1:floor(R/2)
index1=sort(randperm(size(X,1),size_batch));    
index2=setdiff(1:size(X,1),index1);
A_app_array_ind(:,:,:,j)=A_array_ind_TIGRESS_Lasso(X,V,index1,Alpha,alpha_min,lambda_array,L_array);
A_app_array_ind(:,:,:,R+1-j)=A_array_ind_TIGRESS_Lasso(X,V,index2,Alpha,alpha_min,lambda_array,L_array);
end
elapsedTime = toc;

%%Compute scores
[ROC_score_area_mat, PR_score_area_mat, TPR_array_area, FPR_array_area,...
    PPV_array_area, A_app_array_Rnk0]=ScoreMat_TIGRESS(A_app_array_ind, A,L_array,R,bool_orig);

% Max_orig=max(max(ROC_score_orig_mat));
AUROC_GRISLI=max(max(ROC_score_area_mat));
AUPR_GRISLI=max(max(PR_score_area_mat));
%%Export
if saveResultsBool
%     fileID = fopen(saveFileName,'a');
%     fmt = '%6s %4s %4s %5s %5s %6s\r\n';
%     fprintf(fileID,fmt,'#cells','R','L','alpha','AUROC','time');
%     fmt = '%6d \t %4d \t %4d \t %1.3f \t %5.5f \t %6.2f\n';
%     for i=1:length(L_array)
%     fprintf(fileID,fmt, [numcells R L_array(i) alpha_min ROC_score_orig_mat(i) elapsedTime]);
%     end
%     fclose(fileID);
    fileID = fopen(saveFileName,'a');
    % fmt = '%6s %4s %4s %5s %5s %6s\r\n';
    % fprintf(fileID,fmt,'#cells','R','L','alpha','AUROC','time');
    fmt = '%6d \t %4d \t %4d \t %1.3f \t %5.5f \t %6.2f\n';
    for i=1:length(L_array)
    fprintf(fileID,fmt, [numcells R L_array(i) alpha_min ROC_score_area_mat(i) elapsedTime]);
    end
    fclose(fileID);
end
end
end
