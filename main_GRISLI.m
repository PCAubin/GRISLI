warning('off')
addpath(genpath('spams-matlab-v2.6'))
addpath(genpath('SCODE-master'))
warning('on')
%% Retrieve the SCODE datasets
filenumber=2;%Choose the murine data (3-SCODE Data2) or the human data (2-SCODE Data3)
path_name=['SCODE-master/data' num2str(filenumber) '/']; 
data_matrix = dlmread([path_name 'datamatrix.txt'],'\t');
pseudotime_array = dlmread([path_name 'pseudotime.txt'],'\t');
realtime_array = dlmread([path_name 'realtime.txt'],'\t');
A = dlmread([path_name 'A.txt'],'\t');%The real binary graph from the literature
X=[pseudotime_array,data_matrix'];
[~,I]=sort(X(:,1));
X=X(I,:);%X is the spacetime matrix of the data (of size Cx(1+G)), with the
%chosen (real or pseudo) time in the first column and the gene expression in the others 

% dlmwrite('myFile.txt',M,'delimiter','\t','precision',3)
% listnnz_A=1:size(A,1);
% A_wo_diag = A - diag(diag(A));
% listnz_A=intersect(find(~logical(sum(A_wo_diag,1))),find(~logical(sum(A_wo_diag,2))));
% listnnz_A = setdiff(listnnz_A,listnz_A);
% A_wo_diag_red=A_wo_diag(listnnz_A,listnnz_A);
%% TESTING SCODE

% D=4;
% number_tries=100;
% number_average=50;
% number_test=30;
% saveResults=true;
% for i=1:number_test
% [numcells, AUROC_scode, elapsedTime]=Test_Scode_realdata(filenumber,A,D,X,number_tries,number_average,saveResults);
% end

% Results_matrix = dlmread('Compared_perf_GRISLI_vs_SCODE.txt','\t', 1, 0);
  %% Kernel and plot

Alpha=@(Kx,Dt,sigx,sigt)exp(-(Kx.^2)/(2*sigx^2)).*exp(-(Dt.^2)/(2*sigt^2)).*(Dt.^2);%The Peanian kernel we use

% The other GRISLI parameters are :ç_
% R=1000;
L_array=70;%1:15;%1:15; 20:5:90
alpha_min=.4;

numcells=size(X,1);
lambda_array=[0.0001];%In case we wish to use an array of lambda values over which we compute the lasso, rather than using L



K = SpacetimeKernel(X,Alpha);
Knorm = KernelNormalization(K);
V=VelocityInference(X,Knorm);

for count=1:20 %alpha_min=0.05:0.05:0.95
for R=[10, 50,80,100, 200,500, 800,1000,1500,2000,3000, 4000,5000]

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
[ROC_score_orig_mat,ROC_score_area_mat]=...
    ScoreMat_TIGRESS(A_app_array_ind, A,L_array,R);%, TPR_array_area, FPR_array_area,TPR_array_orig, FPR_array_orig

Max_orig=max(max(ROC_score_orig_mat));
Max_area=max(max(ROC_score_area_mat));

% ROCPlot([FPR_array_area;FPR_array_orig],[TPR_array_area;TPR_array_orig],["TIGRESS area","TIGRESS orig"]);
%%Export
fileID = fopen('RGRISLI_orig_perf.txt','a');
% fmt = '%6s %4s %4s %5s %5s %6s\r\n';
% fprintf(fileID,fmt,'#cells','R','L','alpha','AUROC','time');
fmt = '%6d \t %4d \t %4d \t %1.3f \t %5.5f \t %6.2f\n';
for i=1:length(L_array)
fprintf(fileID,fmt, [numcells R L_array(i) alpha_min ROC_score_orig_mat(i) elapsedTime]);
end
fclose(fileID);
fileID = fopen('RGRISLI_area_perf.txt','a');
% fmt = '%6s %4s %4s %5s %5s %6s\r\n';
% fprintf(fileID,fmt,'#cells','R','L','alpha','AUROC','time');
fmt = '%6d \t %4d \t %4d \t %1.3f \t %5.5f \t %6.2f\n';
for i=1:length(L_array)
fprintf(fileID,fmt, [numcells R L_array(i) alpha_min ROC_score_area_mat(i) elapsedTime]);
end
fclose(fileID);
end
end
  %% Export
% tic
% Score_export=[0 L_thr_array; [L_array' , ROC_score_orig_mat]];
% datafilename=['Data2_ROC_orig_alph_' num2str(10*alpha_min) '_R_' num2str(R) '.mat'];
% save(datafilename, 'Score_export');
% Score_export=[0 L_thr_array; [L_array' , ROC_score_area_mat]];
% datafilename=['Data2_ROC_area_alph_' num2str(10*alpha_min) '_R_' num2str(R) '.mat'];
% save(datafilename, 'Score_export');
% toc
% figure
% plot(L_array,ROC_score_area_mat)
% hold on
% plot(L_array,ROC_score_orig_mat)