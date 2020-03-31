function [AUROC_scode, AUPR_scode, elapsedTime,...
    TPR_array_SCODE, FPR_array_SCODE, PPV_array_SCODE,I_scode]=Test_SCODE_realdata(filenumber,A,D,X,...
    number_tries,number_average,saveResultsBool,saveFileName)
% Test_SCODE_realdata computes the AUROC scores of the matrices predicted
% through the SCODE Algorithm (Matsumoto et al., 2017) as compared to the 
% reference matrix A (G*G) for the given dataset found at data*filenumber*/.
% If the boolean saveResultsBool is true it can 
% export the results in the txt file saveFileName under the form:
% [#cells number_tries number_average AUROC_scode elapsedTime]
% Test_SCODE_realdata calls the R function provided by (Matsumoto et al.,
% 2017), it thus requires installing R for it to work.
% For instance SCODE parameters could be:
% D=4;
% number_tries=100;
% number_average=50;
% saveResultsBool=false;
% Pierre-Cyril Aubin-Frankowski, 2018

% IMPORTANT: this function requires Rscript to run as it calls R from
% MATLAB (as SCODE was written in R).

if nargin<8
    saveFileName='SCODE_results.txt';
end

if nargin<7
    saveResultsBool=false;
end

tic
cd SCODE-master
command = ['Rscript SCODE.R data' num2str(filenumber) '/datamatrix.txt data'...
    num2str(filenumber) '/pseudotime.txt out ' num2str(size(A,1)) ' ' num2str(D)...
    ' ' num2str(size(X,1)) ' ' num2str(number_tries) ' ' num2str(number_average)];
system(command);
cd ..
numcells=size(X,1);
A_SCODE = dlmread('SCODE-master/A_SCODE.txt','\t');
elapsedTime = toc;

listnnz_A=1:size(A,1);
A_wo_diag = A - diag(diag(A));
listnz_A=intersect(find(~logical(sum(A_wo_diag,1))),find(~logical(sum(A_wo_diag,2))));
listnnz_A = setdiff(listnnz_A,listnz_A);
A_wo_diag_red=A_wo_diag(listnnz_A,listnnz_A);

[~,I_scode]=sort(abs(A_SCODE(:)),'descend');
[~,I_scode]=sort(I_scode);
I_scode=reshape(I_scode,size(A_SCODE));

% Rnk_array_SCODE=zeros(size(A_SCODE,1),size(A_SCODE,2),numel(A));
TPR_array_SCODE=zeros(1,length(numel(A)+1));
FPR_array_SCODE=zeros(1,length(numel(A)+1));
PPV_array_SCODE=zeros(1,length(numel(A)+1));
for k=1:(numel(A)+1)
   A_app_Rnk=(I_scode<k);
%    Rnk_array_SCODE(:,:,k)=A_app_Rnk;
   [TPR_array_SCODE(k), FPR_array_SCODE(k), PPV_array_SCODE(k)]=CompROC(A_app_Rnk,A_wo_diag_red,listnnz_A);   
end
% [M,I] = max(abs(A_SCODE(:)));
% [I_row, I_col] = ind2sub(size(A_SCODE),I);

AUROC_scode=trapz(FPR_array_SCODE(1:end),TPR_array_SCODE(1:end));
AUPR_scode=trapz(TPR_array_SCODE(1:end),PPV_array_SCODE(1:end));

if saveResultsBool
    fileID = fopen(saveFileName,'a');
    % fmt = '%6s %4s %4s %5s %6s\r\n';
    % fprintf(fileID,fmt,'#cells','L','L_thr','AUROC','time');
    fmt = '%6d \t %4d \t %4d \t %5.5f \t %6.2f\n';
    fprintf(fileID,fmt, [size(X,1) number_tries number_average AUROC_scode elapsedTime]);
    fclose(fileID);
end
end
