% Demo for GRISLI
% Pierre-Cyril Aubin-Frankowski, 2018
%% Please make sure to be located in GRISLI/.
% cd 'C:\Users\pierr\OneDrive\Documents\Biologie_Vert\Code_bio\Code GRISLI v2'
clear all
close all
addpath(genpath('./'))%Add the necessary path
fprintf('This is a demo run of GRISLI intented to reproduce the results of the article.\n')

%% Retrieve the SCODE datasets
filenumber=4;%Choose the murine data (2 for SCODE Data2) or the human data 
%(3 for SCODE Data3) or the pancreatic data (4 for scvelo dataset)
path_name=['SCODE-master/data' num2str(filenumber) '/']; 
data_matrix = dlmread([path_name 'datamatrix.txt'],'\t');
% data_matrix = dlmread([path_name 'data' num2str(filenumber) 'zinb.txt'],'\t');
pseudotime_array = dlmread([path_name 'pseudotime.txt'],'\t');
if filenumber ~=4
	realtime_array = dlmread([path_name 'realtime.txt'],'\t');
    % no real time is available for scvelo data4
end
if filenumber ==4
	Vvelocity_matrix= dlmread([path_name 'velocity.txt'],' ');
    % scevlo velocity is only available for scvelo data4
end

A = dlmread([path_name 'A.txt'],'\t');%The real binary graph from the literature

t_array=pseudotime_array;%Choose the time label (real or pseudo) 
%(Advice: pseudo for SCODE Data2, real for SCODE Data3, pseudo for scvelo Data4)
X=[t_array,data_matrix'];

[~,I]=sort(X(:,1));
X=X(I,:);%X is the spacetime matrix of the data (of size Cx(1+G)), with the
%chosen (real or pseudo) time in the first column and the gene expression in the others 
%% TESTING GRISLI
Alpha=@(Kx,Dt,sigx,sigt)exp(-(Kx.^2)/(2*sigx^2)).*exp(-(Dt.^2)/(2*sigt^2)).*(Dt.^2);%The kernel we use
R=100;
L_array=100;
alpha_min=.4;
saveResultsBool=false;
nbtry=1;
%%
saveFileName='AUROC_files/GRISLI_Data4_Gvelo.txt';
[AUROC_GRISLI, AUPR_GRISLI, elapsedTime, TPR_array_area, FPR_array_area,...
    PPV_array_area,A_app_array_Rnk0]=Test_GRISLI_realdata(A,X,L_array,Alpha,R,alpha_min,saveResultsBool, saveFileName,false);%,V0,nbtry
ROCPlot(FPR_array_area,TPR_array_area,"GRISLI");
% save('Data4_GRISLI_GV','AUROC_GRISLI', 'AUPR_GRISLI', 'TPR_array_area',...
%     'FPR_array_area','PPV_array_area','A_app_array_Rnk0')
%%
saveFileName='AUROC_files/GRISLI_Data4_scvelo.txt';
[AUROC_GRISLI_scV, AUPR_GRISLI_scV, elapsedTime_scV, TPR_array_area_scV, FPR_array_area_scV,...
    PPV_array_area_scV,A_app_array_Rnk0_scV]=Test_GRISLI_realdata(A,X,L_array,Alpha,R,alpha_min,saveResultsBool, saveFileName,false,Vvelocity_matrix,nbtry);
% save('Data4_GRISLI_scV','AUROC_GRISLI_scV', 'AUPR_GRISLI_scV', 'TPR_array_area_scV',...
%     'FPR_array_area_scV','PPV_array_area_scV','A_app_array_Rnk0_scV')
ROCPlot(FPR_array_area_scV,TPR_array_area_scV,"GRISLI + scvelo");
%% TESTING the algorithms IDA (pcalg) and MI (mutual inference)
[AUROC_ida, AUPR_ida, elapsedTime,...
    TPR_array_IDA, FPR_array_IDA, PPV_array_IDA,I_ida]=Test_IDA_realdata(filenumber,A,'min');
[AUROC_MI, AUPR_MI, elapsedTime,...
    TPR_array_MI, FPR_array_MI, PPV_array_MI,I_MI]=Test_MI_realdata(filenumber,A)

%% TESTING SCODE
D=4;
number_tries=20;%100;
number_average=10;%50;
number_test=30;
saveResultsBool=true;
saveFileName='AUROC_files/Data4new_SCODE.txt';

for i=1:number_test
[AUROC_scode, AUPR_scode, elapsedTime, TPR_array_SCODE, FPR_array_SCODE,...
    PPV_array_SCODE,I_scode]=Test_SCODE_realdata(filenumber,A,D,X,number_tries,number_average,saveResultsBool,saveFileName);
end
save('Data4_new_scode','AUROC_scode', 'AUPR_scode', 'TPR_array_SCODE',...
    'FPR_array_SCODE','PPV_array_SCODE','I_scode')
%% TESTING TIGRESS
alpha_min=.4;
saveResults=true;
L_max=100;
L_array=100;%20:10:L_max;
R=200;
number_test=29;
saveResultsBool=true;
saveFileName='AUROC_files/Data4new_TIGRESS.txt';

for count=1:number_test
    [AUROC_TIGRESS,AUPR_TIGRESS, elapsedTime,TPR_array_TIGRESS,...
    FPR_array_TIGRESS, PPV_array_TIGRESS]=Test_TIGRESS_realdata(A,X,L_array,L_max,R,alpha_min,saveResultsBool, saveFileName);
end
% save('Data4_new_TIGRESS','AUROC_TIGRESS', 'AUPR_TIGRESS', 'TPR_array_TIGRESS',...
%     'FPR_array_TIGRESS','PPV_array_TIGRESS')
%% RETRIEVING THE BEST INTERACTIONS FOR GRISLI
fileID = fopen([path_name 'tf.txt'],'r');
TF_names = textscan(fileID,'%s \t %s'); TF_names=TF_names{1};
fclose(fileID);
% load('Data4_GRISLI_scV')
% Replace A_app_array_Rnk0 by I_ida or I_MI for non-GRISLI procedures
% A_app_small=I_ida(:,:);
k0=20; A_app_small=A_app_array_Rnk0(:,:); A_app_small(A_app_small>k0)=0;%_scV
TF_names_A_Rnk0=cell(k0,2);

saveFileName=['TFs_Data2_A.txt'];%num2str(filenumber)
[row,col] = find(A_app_small);
for k=1:k0
    TF_names_A_Rnk0{A_app_small(row(k),col(k)),1}=TF_names{row(k)};
    TF_names_A_Rnk0{A_app_small(row(k),col(k)),2}=TF_names{col(k)};
end
[nrows,ncols] = size(TF_names_A_Rnk0);
fileID = fopen(saveFileName,'a');
for rowIdx = 1:nrows
    fprintf(fileID,'%s \t & %s \\\\ \n',TF_names_A_Rnk0{rowIdx,:});
end
fclose(fileID);

%% ROC and PRC plots of comapred performance
load('Data4_GRISLI_new_GV');load('Data4_GRISLI_new_scV');load('Data4_new_scode');load('Data4_new_TIGRESS');

FPR_array_multiple=[FPR_array_area;FPR_array_SCODE;FPR_array_TIGRESS;FPR_array_area_scV];
TPR_array_multiple=[TPR_array_area;TPR_array_SCODE;TPR_array_TIGRESS;TPR_array_area_scV];
PPV_array_multiple=[PPV_array_area;PPV_array_SCODE;PPV_array_TIGRESS;PPV_array_area_scV];
ROCPlot(FPR_array_multiple,TPR_array_multiple,["GRISLI","SCODE","TIGRESS","scvelo+GRISLI"]); %ROC curves
PRPlot(TPR_array_multiple,PPV_array_multiple,["GRISLI","SCODE","TIGRESS","scvelo+GRISLI"],0.15);
PRPlot(TPR_array_multiple,PPV_array_multiple,["GRISLI","SCODE","TIGRESS","scvelo+GRISLI"],1);
% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'PR_Data4_compa_small','-dpdf')


%Plotting GRISLI results as in (Aubin-Frankowski and Vert, 2020)
%% Boxplot of compared performances
addpath(genpath('AUROC_boxplots'))
% Data2_SCODE_matrix=dlmread('AUROC_boxplots/Data2_SCODE_perf.txt','\t', 0, 0);
% Data2_TIGRESS_matrix=dlmread('AUROC_boxplots/Data2_TIGRESS_perf.txt','\t', 0, 0);
% Data2_GRISLI_matrix=dlmread('AUROC_boxplots/Data2_GRISLI_area_perf.txt','\t', 0, 0);
% 
% Data3_SCODE_matrix=dlmread('AUROC_boxplots/Data3_SCODE_perf.txt','\t', 0, 0);
% Data3_TIGRESS_matrix=dlmread('AUROC_boxplots/Data3_TIGRESS_perf.txt','\t', 0, 0);
% Data3_GRISLI_matrix=dlmread('AUROC_boxplots/Data3_GRISLI_area_perf.txt','\t', 0, 0);

Data4_GRISLI_scvelo_matrix=dlmread('AUROC_files/Data4new_scvelo.txt','\t', 0, 0);
Data4_GRISLI_Gvelo_matrix=dlmread('AUROC_files/Data4new_Gvelo.txt','\t', 0, 0);
Data4_SCODE_matrix=dlmread('AUROC_files/Data4new_SCODE.txt','\t', 0, 0);
Data4_TIGRESS_matrix=dlmread('AUROC_files/Data4new_TIGRESS.txt','\t', 0, 0);

% maxsize=max([size(Data2_GRISLI_matrix,1),size(Data2_TIGRESS_matrix,1),size(Data2_SCODE_matrix,1),...
%     size(Data3_GRISLI_matrix,1),size(Data3_TIGRESS_matrix,1),size(Data3_SCODE_matrix,1)]);
% 
% z2S=nan(maxsize,2);
% z2S(1:size(Data2_SCODE_matrix,1),:)=Data2_SCODE_matrix(:,end-1:end);
% z2A=nan(maxsize,2);
% z2A(1:size(Data2_GRISLI_matrix,1),:)=Data2_GRISLI_matrix(:,end-1:end);
% z2T=nan(maxsize,2);
% z2T(1:size(Data2_TIGRESS_matrix,1),:)=Data2_TIGRESS_matrix(:,end-1:end);
% 
% z3S=nan(maxsize,2);
% z3S(1:size(Data3_SCODE_matrix,1),:)=Data3_SCODE_matrix(:,end-1:end);
% z3A=nan(maxsize,2);
% z3A(1:size(Data3_GRISLI_matrix,1),:)=Data3_GRISLI_matrix(:,end-1:end);
% z3T=nan(maxsize,2);
% z3T(1:size(Data3_TIGRESS_matrix,1),:)=Data3_TIGRESS_matrix(:,end-1:end);
maxsize=30;
z4S=nan(maxsize,2);
z4S(1:size(Data4_SCODE_matrix,1),:)=Data4_SCODE_matrix(:,end-1:end);
z4A=nan(maxsize,2);
z4A(1:size(Data4_GRISLI_scvelo_matrix,1),:)=Data4_GRISLI_scvelo_matrix(:,end-1:end);
z4B=nan(maxsize,2);
z4B(1:size(Data4_GRISLI_Gvelo_matrix,1),:)=Data4_GRISLI_Gvelo_matrix(:,end-1:end);
z4T=nan(maxsize,2);
z4T(1:size(Data4_TIGRESS_matrix,1),:)=Data4_TIGRESS_matrix(:,end-1:end);

figure
% subplot(1,2,1)
% boxplot([z2S(:,1),z2T(:,1),z2A(:,1)],'Labels',{'SCODE','TIGRESS','GRISLI'})
%  h = findobj(gca,'Tag','Box');
%  set(gca,'fontsize',20)
%  set(findobj(gca,'type','line'),'linew',2)
%  color_list='y';
%  ylim([0.5 0.58])
%  for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'), color_list,'FaceAlpha',.2);
%  end
% 
%  hold on
%  title('Murine dataset')
%  subplot(1,2,2)

% boxplot([z3S(:,1),z3T(:,1),z3A(:,1)],'Labels',{'SCODE','TIGRESS','GRISLI'})
boxplot([z4S(:,1),z4T(:,1),z4A(:,1),z4B(:,1)],'Labels',{'SCODE','TIGRESS','scvelo+GRISLI','GRISLI'})
%  title('Human dataset')
title('Pancreatic (murine) dataset')
  set(gca,'fontsize',10)
 h = findobj(gca,'Tag','Box');
 set(findobj(gca,'type','line'),'linew',1)
 color_list='y';
 ylim([0.49 0.55])
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'), color_list,'FaceAlpha',.2);
 end
fig = gcf; fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'Data4_boxplot_compa','-dpdf')

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

numcells=373;
Results_matrix = dlmread('AUROC_files/GRISLI_Data2-3_varyingL_alpha.txt','\t', 2,0);
Results_matrix_d = Results_matrix(Results_matrix(:,1)==numcells,:);
Results_matrix_d=Results_matrix_d(1:285,:);
L_array=sort(unique(Results_matrix_d(:,3)))';
alpha_array=sort(unique(Results_matrix_d(:,4)))';

R_alpha=reshape(Results_matrix_d(:,5),length(L_array),length(alpha_array));

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
yticklabels({'20','40','65','90'})