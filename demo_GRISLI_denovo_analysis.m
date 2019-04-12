% Demo for GRISLI
% Pierre-Cyril Aubin-Frankowski, 2018
%% Please make sure to be located in GRISLI/.

clear all
close all
addpath(genpath('./'))%Add the necessary path
fprintf('This is a demo to run GRISLI on a dataset and produce experimental regulatory networks.\n')

%% Retrieve the SCODE datasets || OR plug in your own data !

filenumber=2;%Choose the murine data (2 for SCODE Data2) or the human data (3 for SCODE Data3)
path_name=['SCODE-master/data' num2str(filenumber) '/']; 
data_matrix = dlmread([path_name 'datamatrix.txt'],'\t');
pseudotime_array = dlmread([path_name 'pseudotime.txt'],'\t');
realtime_array = dlmread([path_name 'realtime.txt'],'\t');
A = dlmread([path_name 'A.txt'],'\t');%The real binary graph from the literature

t_array=pseudotime_array;%Choose the time label (real or pseudo) 
%(Advice: pseudo for SCODE Data2, real for SCODE Data3)

X=[pseudotime_array,data_matrix'];
[~,I]=sort(X(:,1));
X=X(I,:);%X is the spacetime matrix of the data (of size Cx(1+G)), with the
%chosen (real or pseudo) time in the first column and the gene expression in the others 

%% TESTING GRISLI

Alpha=@(Kx,Dt,sigx,sigt)exp(-(Kx.^2)/(2*sigx^2)).*exp(-(Dt.^2)/(2*sigt^2)).*(Dt.^2);%The kernel we use
R=100;
L_array=20:30:90;
alpha_min=.3;

%Rnk_array_TIGRESS_area_L (G*G*length(L_array)) contains the predicted edges
%of the network. Its last index is related to L=L_array(i).
Rnk_array_TIGRESS_area_L=Compute_A_app_wo_ref(X,L_array,Alpha,...
    R,alpha_min);