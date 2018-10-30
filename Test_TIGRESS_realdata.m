function [AUROC_score_area]=Test_TIGRESS_realdata(A,X,L_array,L_max,R,alpha_min,saveResultsBool, saveFileName)
% Test_TIGRESS_realdata computes the AUROC scores of the matrices predicted
% through area TIGRESS (Haury et al. 2012) compared to the reference matrix A (G*G). 
% If the boolean saveResultsBool is true it can export the results in the txt file
% saveFileName under the form:
% [#cells R L alpha_min AUROC_score_area elapsedTime]
% L_array contains the L values to investigate. L_max has to be larger than
% largest value in L_array.
% Test_TIGRESS_realdata calls the Matlab functions provided by (Haury et al. 2012).
% For instance TIGRESS parameters could be:
% saveResults=false;
% L_max=90;
% L_array=1:5:L_max;
% alpha_min=.2
% R=1500;
% Pierre-Cyril Aubin-Frankowski, 2018

if nargin<8
    saveFileName='TIGRESS_results.txt';
end

if nargin<7
    saveResultsBool=false;
end


data.expdata=X(:,2:end);
genenames = strsplit(int2str(1:100));
data.genenames=genenames';

tic
freq=tigress(data,'R',R,'alpha',alpha_min,'L',L_max,'LarsAlgo', 'spams', 'verbose',false);
for L=L_array
Score_mat_area=score_edges(freq,'method','area','L',L);
Score_mat_area=Score_mat_area';
% edges = predict_network(Score_mat_area,1:100,'genenames',data.genenames);

% Score_mat_orig=score_edges(freq,'method','original','L',L);


A_app_ind=A;
% Sort the scores and store in I their ranks with respect to the whole 
% matrix
[~,I_area]=sort(Score_mat_area(:),'descend');
[~,I_area]=sort(I_area);
I_area=reshape(I_area,size(Score_mat_area));

% [~,I_orig]=sort(Score_mat_orig(:),'descend');
% [~,I_orig]=sort(I_orig);
% I_orig=reshape(I_orig,size(Score_mat_orig));

Rnk_array_TIGRESS_area=zeros(size(A_app_ind,1),size(A_app_ind,2),numel(A));
% Rnk_array_TIGRESS_orig=Rnk_array_TIGRESS_area;

% TPR_array_orig=zeros(1,length(numel(A)+1));
% FPR_array_orig=zeros(1,length(numel(A)+1));

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
 
%    A_app_Rnk=(I_orig<k);
%    Rnk_array_TIGRESS_orig(:,:,k)=A_app_Rnk;
%    [TPR_array_orig(k), FPR_array_orig(k)]=CompROC(A_app_Rnk,A_wo_diag_red,listnnz_A);
end

% AUROC_score_orig=trapz(FPR_array_orig(1:end),TPR_array_orig(1:end));
AUROC_score_area=trapz(FPR_array_area(1:end),TPR_array_area(1:end));
elapsedTime = toc;
% ROCPlot([FPR_array_area;FPR_array_orig],[TPR_array_area;TPR_array_orig],["TIGRESS area","TIGRESS orig"]);
if saveResultsBool
    fileID = fopen(saveFileName,'a');
    % fmt = '%6s %4s %4s %5s %6s\r\n';
    % fprintf(fileID,fmt,'#cells','L','L_thr','AUROC','time');
    fmt = '%6d \t %4d \t %4d \t %1.3f \t %5.5f \t %6.2f\n';
    fprintf(fileID,fmt,[size(X,1) R L alpha_min AUROC_score_area elapsedTime]);
    fclose(fileID);
end
end
end
