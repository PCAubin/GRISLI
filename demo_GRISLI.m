% Demo for GRISLI
% Pierre-Cyril Aubin-Frankowski, 2018
%% Please make sure to be located in GRISLI/.

clear all
close all
addpath(genpath('./'))%Add the necessary path
fprintf('This is a demo to run GRISLI.\n')

%% Retrieve the SCODE datasets

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
L_array=70;%20:10:90;
alpha_min=.3;
saveResultsBool=false;
saveFileName='AUROC_files/GRISLI_results.txt';

[AUROC_GRISLI, elapsedTime, TPR_array_area, FPR_array_area]=Test_GRISLI_realdata(A,X,L_array,Alpha,...
    R,alpha_min,saveResultsBool, saveFileName);
ROCPlot(FPR_array_area,TPR_array_area,"GRISLI"); %ROC curves

%% TESTING SCODE

D=4;
number_tries=100;
number_average=50;
number_test=1;
saveResultsBool=false;
saveFileName='AUROC_files/SCODE_results.txt';

for i=1:number_test
[AUROC_SCODE, elapsedTime]=Test_SCODE_realdata(filenumber,A,D,X,number_tries,number_average,saveResultsBool,saveFileName);
end

%% TESTING TIGRESS

alpha_min=.3;
saveResults=true;
L_max=90;
L_array=90;%20:10:L_max;
R=100;
number_test=1;
saveResultsBool=false;
saveFileName='AUROC_files/TIGRESS_results.txt';

for count=1:number_test
    [AUROC_TIGRESS, elapsedTime]=Test_TIGRESS_realdata(A,X,L_array,L_max,R,alpha_min,saveResultsBool, saveFileName);
end

%Plotting GRISLI results as in (Aubin-Frankowski and Vert, 2018)
%% Boxplot of compared performances
addpath(genpath('AUROC_boxplots'))
Data2_SCODE_matrix=dlmread('AUROC_boxplots/Data2_SCODE_perf.txt','\t', 0, 0);
Data2_TIGRESS_matrix=dlmread('AUROC_boxplots/Data2_TIGRESS_perf.txt','\t', 0, 0);
Data2_GRISLI_matrix=dlmread('AUROC_boxplots/Data2_GRISLI_area_perf.txt','\t', 0, 0);

Data3_SCODE_matrix=dlmread('AUROC_boxplots/Data3_SCODE_perf.txt','\t', 0, 0);
Data3_TIGRESS_matrix=dlmread('AUROC_boxplots/Data3_TIGRESS_perf.txt','\t', 0, 0);
Data3_GRISLI_matrix=dlmread('AUROC_boxplots/Data3_GRISLI_area_perf.txt','\t', 0, 0);

maxsize=max([size(Data2_GRISLI_matrix,1),size(Data2_TIGRESS_matrix,1),size(Data2_SCODE_matrix,1),...
    size(Data3_GRISLI_matrix,1),size(Data3_TIGRESS_matrix,1),size(Data3_SCODE_matrix,1)]);

z2S=nan(maxsize,2);
z2S(1:size(Data2_SCODE_matrix,1),:)=Data2_SCODE_matrix(:,end-1:end);
z2A=nan(maxsize,2);
z2A(1:size(Data2_GRISLI_matrix,1),:)=Data2_GRISLI_matrix(:,end-1:end);
z2T=nan(maxsize,2);
z2T(1:size(Data2_TIGRESS_matrix,1),:)=Data2_TIGRESS_matrix(:,end-1:end);

z3S=nan(maxsize,2);
z3S(1:size(Data3_SCODE_matrix,1),:)=Data3_SCODE_matrix(:,end-1:end);
z3A=nan(maxsize,2);
z3A(1:size(Data3_GRISLI_matrix,1),:)=Data3_GRISLI_matrix(:,end-1:end);
z3T=nan(maxsize,2);
z3T(1:size(Data3_TIGRESS_matrix,1),:)=Data3_TIGRESS_matrix(:,end-1:end);

figure
subplot(1,2,1)
boxplot([z2S(:,1),z2T(:,1),z2A(:,1)],'Labels',{'SCODE','TIGRESS','GRISLI'})
 h = findobj(gca,'Tag','Box');
 set(gca,'fontsize',20)
 set(findobj(gca,'type','line'),'linew',2)
 color_list='y';
 ylim([0.5 0.58])
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'), color_list,'FaceAlpha',.2);
 end

 hold on
 title('Murine dataset')
 subplot(1,2,2)

boxplot([z3S(:,1),z3T(:,1),z3A(:,1)],'Labels',{'SCODE','TIGRESS','GRISLI'})
 title('Human dataset')
  set(gca,'fontsize',20)
 h = findobj(gca,'Tag','Box');
 set(findobj(gca,'type','line'),'linew',2)
 color_list='y';
 ylim([0.5 0.58])
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'), color_list,'FaceAlpha',.2);
 end

%% Effect of R on the AUROC perfomance
RArea_Results_matrix = dlmread('AUROC_files/GRISLI_Data2_varyingR.txt','\t');
R_array=sort(unique(RArea_Results_matrix(:,2)))';
%[10, 50,80,100, 200,500, 800,1000,1500,2000,3000,4000,5000];

mat_to_plot=nan(size(RArea_Results_matrix,1),length(R_array));
mat_to_plot(any(isnan(mat_to_plot),2),:) = [];
for i=1:length(R_array)
    data=RArea_Results_matrix(RArea_Results_matrix(:,2)==R_array(i),5);
    mat_to_plot(1:size(data,1),i)=data;    
end    
mat_to_plot(mat_to_plot==0)=NaN;
figure
hold on
line(1:(length(R_array)),nanmean(mat_to_plot))
boxplot(mat_to_plot,'Labels',...
    {num2str(R_array')})
xlabel('R')
ylabel('AUROC')
% ylim([0.5 0.58])

%% Heatmap of values for varying GRISLI/TIGRESS parameters alpha_min and L

numcells=758;
Results_matrix = dlmread('AUROC_files/GRISLI_Data2-3_varyingL_alpha.txt','\t', 2,0);
Results_matrix_d = Results_matrix(Results_matrix(:,1)==numcells,:);
L_array=sort(unique(Results_matrix_d(:,3)))';
alpha_array=sort(unique(Results_matrix_d(:,4)))';

R_alpha=reshape(Results_matrix_d(1:285,5),length(L_array),length(alpha_array));

figure
imagesc(R_alpha)
colorbar
colormap(jet(256))
caxis([0.5 0.58])
axis xy
% axis([0 1 1 70])
xlabel('\alpha')
ylabel('L')
xticks([1 5 10 15 20])
yticks([1 5 10 15])
xticklabels({'0','0.25','0.5','0.75','1'})
% yticklabels({'20','40','65','90'})


