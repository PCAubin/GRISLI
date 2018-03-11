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

D=4;
number_tries=100;
number_average=50;
number_test=30;
saveResults=true;
for i=1:number_test
[numcells, AUROC_scode, elapsedTime]=Test_Scode_realdata(filenumber,A,D,X,number_tries,number_average,saveResults);
end

% Results_matrix = dlmread('Compared_perf_GRISLI_vs_SCODE.txt','\t', 1, 0);
  %% Kernel and plot

Alpha=@(Kx,Dt,sigx,sigt)exp(-(Kx.^2)/(2*sigx^2)).*exp(-(Dt.^2)/(2*sigt^2)).*(Dt.^2);%The Peanian kernel we use
tic
% The other GRISLI parameters are :
R=500;
L_array=75; 
L_thr_array=60;
alpha_min=.4;

numcells=size(X,1);
lambda_array=[0.0001];%In case we wish to use an array of lambda values over which we compute the lasso, rather than using L

%4-D array to store the ranks of the inferred A of the R tests for each L in L_array
A_app_array_ind=uint32(zeros(size(A,1),size(A,2),length(L_array),R));
size_batch=floor(numcells/2);

for j=1:floor(R/2)
index1=sort(randperm(size(X,1),size_batch));    
index2=setdiff(1:size(X,1),index1);
A_app_array_ind(:,:,:,j)=A_array_ind_TIGRESS_Lasso(X,index1,Alpha,alpha_min,lambda_array,L_array);
A_app_array_ind(:,:,:,R+1-j)=A_array_ind_TIGRESS_Lasso(X,index2,Alpha,alpha_min,lambda_array,L_array);
end
toc

%%Compute scores
[ROC_score_orig_mat,ROC_score_area_mat, TPR_array_area, FPR_array_area,TPR_array_orig, FPR_array_orig]=...
    ScoreMat_TIGRESS(A_app_array_ind, A, L_array,L_thr_array,R);

Max_orig=max(max(ROC_score_orig_mat));
Max_area=max(max(ROC_score_area_mat));

ROCPlot([FPR_array_area;FPR_array_orig],[TPR_array_area;TPR_array_orig],["TIGRESS area","TIGRESS orig"]);
%% Export
tic
Score_export=[0 L_thr_array; [L_array' , ROC_score_orig_mat]];
datafilename=['Data2_ROC_orig_alph_' num2str(10*alpha_min) '_R_' num2str(R) '.mat'];
save(datafilename, 'Score_export');
Score_export=[0 L_thr_array; [L_array' , ROC_score_area_mat]];
datafilename=['Data2_ROC_area_alph_' num2str(10*alpha_min) '_R_' num2str(R) '.mat'];
save(datafilename, 'Score_export');
toc
