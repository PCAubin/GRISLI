function [AUROC_scode, elapsedTime]=Score_SCODE_realdata(filenumber,A,D,X,number_tries,number_average,saveResults)
% Score_SCODE_realdata evaluates the performance of SCODE for the given
% dataset found at data*filenumber*/ based on the pseudotime of the
% measurements. If saveResults=true, then it saves the runtime
% (elapsedTime) and the AUROC score in Score_SCODE_realdata.txt. 
% IMPORTANT: this function requires Rscript to run as it calls R from
% MATLAB (SCODE was written in R).
tic
cd SCODE-master
command = ['Rscript SCODE.R data' num2str(filenumber) '/datamatrix.txt data'...
    num2str(filenumber) '/pseudotime.txt out ' num2str(size(A,1)) ' ' num2str(D)...
    ' ' num2str(size(X,1)) ' ' num2str(number_tries) ' ' num2str(number_average)];
system(command);
cd ..
numcells=size(X,1);
A_SCODE = dlmread('SCODE-master/A.txt','\t');
elapsedTime = toc;

listnnz_A=1:size(A,1);
A_wo_diag = A - diag(diag(A));
listnz_A=intersect(find(~logical(sum(A_wo_diag,1))),find(~logical(sum(A_wo_diag,2))));
listnnz_A = setdiff(listnnz_A,listnz_A);
A_wo_diag_red=A_wo_diag(listnnz_A,listnnz_A);

[~,I_scode]=sort(abs(A_SCODE(:)),'descend');
[~,I_scode]=sort(I_scode);
I_scode=reshape(I_scode,size(A_SCODE));

Rnk_array_SCODE=zeros(size(A_SCODE,1),size(A_SCODE,2),numel(A));
TPR_array_SCODE=zeros(1,length(numel(A)+1));
FPR_array_SCODE=zeros(1,length(numel(A)+1));
for k=1:(numel(A)+1)
   A_app_Rnk=(I_scode<k);
   Rnk_array_SCODE(:,:,k)=A_app_Rnk;
   [TPR_array_SCODE(k), FPR_array_SCODE(k)]=CompROC(A_app_Rnk,A_wo_diag_red,listnnz_A);   
end

AUROC_scode=trapz(FPR_array_SCODE(1:end),TPR_array_SCODE(1:end));

if saveResults
    fileID = fopen('Score_SCODE_realdata.txt','a');
    fmt = '%6d \t %4d \t %4d \t %5.5f \t %6.2f\n';
    fprintf(fileID,fmt, [size(X,1) number_tries number_average AUROC_scode elapsedTime]);
    fclose(fileID);
end
end
