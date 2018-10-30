function A_app_rank_appearance=A_array_ind_TIGRESS_Lasso(X,V,index,Alpha,alpha_min,lambda_array,L_array)
% A_array_ind_TIGRESS_Lasso computes the rank matrix, A_app_ind (G*G) based on the
% lasso-approximated A matrix, A_app (G*G). It takes as inputs the
% spacetime matrix, X (G*(C+1)), an index list, index, which is a subset of [1:C],
% the kernel, Alpha, the lower threshold of the randomizing factor,
% alpha_min, the array of lambda, lambda_array, over which the lasso is
% performed, or the array of maximum steps L, L_array, which supersedes the
% lambda approach.
% Pierre-Cyril Aubin-Frankowski, 2018


% Extract the submatrix X_red out of X
X_red=X(index,:);
V_red=V(index,:);
% Compute the miflux matrix V_red of X_red based on the kernel Alpha
% K = SpacetimeKernel(X_red,Alpha);
% Knorm = KernelNormalization(K);
% V_red=VelocityInference(X_red,Knorm);

% Randomize the lines of X_red
beta=alpha_min+(1-alpha_min)*rand(1,size(X_red,2)-1);
X_red(:,2:end)=beta.*X_red(:,2:end);

% Perform the lasso on (X_red, V_red) with lambda or L as parameters, then
% sort the absolute values of the inferred matrix A_app to get the ranks of
% each edge on its row, which area stored in A_app_ind.
A_app_array_ind=lasso_lambda_L(V_red,X_red,lambda_array,L_array);
A_app_rank_appearance=zeros(size(A_app_array_ind));
for i=1:length(L_array)
[~,I]=sort(A_app_array_ind(:,:,i),2,'descend');
[~,A_app_rank] = sort(uint32(I),2);
A_app_rank(A_app_rank>L_array(i))=0;
A_app_rank_appearance(:,:,i)=A_app_rank;
end
end