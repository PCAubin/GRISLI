function Rnk_array_TIGRESS_area_L=Compute_A_app_wo_ref(X,L_array,Alpha,...
    R,alpha_min)
% Test_GRISLI_realdata computes the AUROC scores of the matrices predicted
% through the GRISLI Algorithm (Aubin-Frankowski and Vert, 2018) as compared to the 
% reference matrix A (G*G) for the given dataset found at data*filenumber*/.
% For instance GRISLI parameters could be:
% R=1500;
% L_array=20:10:90;
% alpha_min=.3;
% saveResultsBool=false;
% Pierre-Cyril Aubin-Frankowski, 2018

numcells=size(X,1);
lambda_array=[0.0001];%In case we wish to use an array of lambda values over
%which we compute the lasso, rather than using L

K = SpacetimeKernel(X,Alpha);
Knorm = KernelNormalization(K);
V=VelocityInference(X,Knorm);

tic
%%4-D array to store the ranks of the inferred A of the R tests for each L in L_array
A_app_array_ind=uint32(zeros(size(X,2)-1,size(X,2)-1,length(L_array),R));
size_batch=floor(numcells/2);
%%Compute frequencies
for j=1:floor(R/2)
index1=sort(randperm(size(X,1),size_batch));    
index2=setdiff(1:size(X,1),index1);
A_app_array_ind(:,:,:,j)=A_array_ind_TIGRESS_Lasso(X,V,index1,Alpha,alpha_min,lambda_array,L_array);
A_app_array_ind(:,:,:,R+1-j)=A_array_ind_TIGRESS_Lasso(X,V,index2,Alpha,alpha_min,lambda_array,L_array);
end


%%Compute ranks
Rnk_array_TIGRESS_area_L=zeros(size(X,2)-1,size(X,2)-1,length(L_array));
for i=1:length(L_array)
    A_app_ind=double(A_app_array_ind(:,:,i,:));
    % Compute the TIGRESS scores
    L=L_array(i);
    Score_mat_area=sum(A_app_ind.*(L+1-A_app_ind),4)/R/L;   
    % Sort the scores and store in I their ranks with respect to the whole 
    % matrix
    [~,I_area]=sort(Score_mat_area(:),'descend');
    [~,I_area]=sort(I_area);
    Rnk_array_TIGRESS_area_L(:,:,i)=reshape(I_area,size(Score_mat_area));    
end
elapsedTime = toc;
disp(['It took ', num2str(elapsedTime), 's to generate the ',...
            num2str(length(L_array)), ' experimental matrices, with R=',num2str(R)])

end