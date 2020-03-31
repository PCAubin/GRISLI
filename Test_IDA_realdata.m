function [AUROC_ida, AUPR_ida, elapsedTime,...
    TPR_array_IDA, FPR_array_IDA, PPV_array_IDA,I_ida]=Test_IDA_realdata(filenumber,A,IDA_selection_crit)
% Test_IDA_realdata computes the AUROC scores of the matrices predicted
% through the IDA algorithm (pcalg package) as compared to the 
% reference matrix A (G*G) for the given dataset found at data*filenumber*/.
% Unlike the SCODE testing procedure where R was called within Matlab, the
% R files to run IDA are to be found in SCODE-master. The follwing code
% only reads the pre-computed matrices without doing any computation. 
% As IDA sometimes outputs multiples values for a given edge, one has to
% select a value among them.
% There are two choices implemented to select the edges based on the outputs
% of IDA. IDA_selection_crit = 'min' or 'mean'. They give similar results.
% Pierre-Cyril Aubin-Frankowski, 2020

tic

A_IDA = dlmread(['SCODE-master/data' num2str(filenumber) '/A_IDA_' IDA_selection_crit '.txt'],'\t');
elapsedTime = toc;

listnnz_A=1:size(A,1);
A_wo_diag = A - diag(diag(A));
listnz_A=intersect(find(~logical(sum(A_wo_diag,1))),find(~logical(sum(A_wo_diag,2))));
listnnz_A = setdiff(listnnz_A,listnz_A);
A_wo_diag_red=A_wo_diag(listnnz_A,listnnz_A);

A_IDA = A_IDA - diag(diag(A_IDA));
[~,I_ida]=sort(abs(A_IDA(:)),'descend');
[~,I_ida]=sort(I_ida);
I_ida=reshape(I_ida,size(A_IDA));

TPR_array_IDA=zeros(1,length(numel(A)+1));
FPR_array_IDA=zeros(1,length(numel(A)+1));
PPV_array_IDA=zeros(1,length(numel(A)+1));
for k=1:(numel(A)+1)
   A_app_Rnk=(I_ida<k);
   [TPR_array_IDA(k), FPR_array_IDA(k), PPV_array_IDA(k)]=CompROC(A_app_Rnk,A_wo_diag_red,listnnz_A);   
end

AUROC_ida=trapz(FPR_array_IDA(1:end),TPR_array_IDA(1:end));
AUPR_ida=trapz(TPR_array_IDA(1:end),PPV_array_IDA(1:end));

% if saveResultsBool
%     fileID = fopen(saveFileName,'a');
%     % fmt = '%6s %4s %4s %5s %6s\r\n';
%     % fprintf(fileID,fmt,'#cells','L','L_thr','AUROC','time');
%     fmt = '%6d \t %4d \t %4d \t %5.5f \t %6.2f\n';
%     fprintf(fileID,fmt, [size(X,1) number_tries number_average AUROC_ida elapsedTime]);
%     fclose(fileID);
% end
end
