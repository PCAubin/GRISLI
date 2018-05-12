%% alpha vs L
list_alpha=0.05:0.05:0.95;
L_array=1:15;%20:5:90;
Results_matrix = dlmread('RGRISLI_orig_perf.txt','\t', 2,0);
Results_matrix_d = Results_matrix(Results_matrix(:,1)==758,:);
R_alpha=reshape(Results_matrix_d(:,5),length(L_array),length(list_alpha));

% Alph_Score_mat_area=zeros(length(L_array),length(list_alpha));
% count=1;
% R=1000;
% for alpha_min=list_alpha
%     datafilename=['Alph_Data3_ROC_area_alph_' num2str(10*alpha_min) '_R_' num2str(R) '.mat'];
%     load(datafilename);
%     Alph_Score_mat_area(:,count)=Score_export(1:71,2);
%     count=count+1;
% end
% Alph_Score_mat_area=[Score_export(1:71,1),Alph_Score_mat_area];
% Alph_Score_mat_area(1,2:end)=list_alpha;
% [M,I] = max(Alph_Score_mat_area(:));
% [I_row, I_col] = ind2sub(size(Alph_Score_mat_area),I);
% datafilename='Data3_ROC_area_alph_vs_lthr.mat';
% load(datafilename,'Alph_Score_mat_area');
% Alph_Score_mat_area=Score_export;

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
%% 
Score_export(Score_export==0)=NaN;
% L_thr_list = cellstr(num2str(Score_export(1,2:end)));
L_thr_list=sprintfc('%u',Score_export(1,2:end));
L_list=sprintfc('%u',Score_export(2:end,1));
figure
pcolor(Score_export(2:end,2:end))
colorbar
% yticks(1:40)
% yticklabels(L_list)
xticks(1:14)
xticklabels(L_thr_list(1:1:end))
xlabel('L_{thr}')
ylabel('L')
%%
list_R=[10,50,100,200,400, 500, 800,1000,1500,2000];
R_Score_list_area=zeros(1,length(R));
R_Score_list_orig=zeros(1,length(R));
count=1;
for R=list_R
    datafilename=['Data2_ROC_area_alph_4_R_' num2str(R) '.mat'];
    load(datafilename);
    R_Score_list_area(:,count)=max(max(Score_export(2:end,2:end)));
    datafilename=['Data2_ROC_orig_alph_4_R_' num2str(R) '.mat'];
    load(datafilename);
    R_Score_list_orig(:,count)=max(max(Score_export(2:end,2:end)));
    count=count+1;
end
list_R=[list_R 5000];
R_Score_list_orig=[R_Score_list_orig, 0.5461];
R_Score_list_area=[R_Score_list_area, 0.5524];
figure
semilogx(list_R,R_Score_list_orig)
hold on
semilogx(list_R,R_Score_list_area)
grid on
ylim([.5, inf])
xlabel('R')
ylabel('AUROC')
legend('original Tigress','area Tigress')

%%

SCODE_Results_matrix = dlmread('Compared_perf_GRISLI_vs_SCODE.txt','\t', 1, 0);
Area_Results_matrix = dlmread('GRISLI_area_perf.txt','\t', 0, 0);
Orig_Results_matrix = dlmread('GRISLI_orig_perf.txt','\t', 0, 0);

Data2_SCODE_matrix=SCODE_Results_matrix(find(SCODE_Results_matrix(:,1)==373),:);
Data2_SCODE_matrix=Data2_SCODE_matrix(find(Data2_SCODE_matrix(:,3)==50),:);

Data2_Area_matrix=Area_Results_matrix(find(Area_Results_matrix(:,1)==373),:);

Data2_Orig_matrix=Orig_Results_matrix(find(Orig_Results_matrix(:,1)==373),:);

Data3_SCODE_matrix=SCODE_Results_matrix(find(SCODE_Results_matrix(:,1)==758),:);
Data3_SCODE_matrix=Data3_SCODE_matrix(find(Data3_SCODE_matrix(:,3)==50),:);

Data3_Area_matrix=Area_Results_matrix(find(Area_Results_matrix(:,1)==758),:);

Data3_Orig_matrix=Orig_Results_matrix(find(Orig_Results_matrix(:,1)==758),:);


maxsize=max([size(Data2_Area_matrix,1),size(Data2_Orig_matrix,1),size(Data2_SCODE_matrix,1),...
    size(Data3_Area_matrix,1),size(Data3_Orig_matrix,1),size(Data3_SCODE_matrix,1)]);
%%
z2S=nan(maxsize,2);
z2S(1:size(Data2_SCODE_matrix,1),:)=Data2_SCODE_matrix(:,end-1:end);

z2A=nan(maxsize,2);
z2A(1:size(Data2_Area_matrix,1),:)=Data2_Area_matrix(:,end-1:end);

z2O=nan(maxsize,2);
z2O(1:size(Data2_Orig_matrix,1),:)=Data2_Orig_matrix(:,end-1:end);

z3S=nan(maxsize,2);
z3S(1:size(Data3_SCODE_matrix,1),:)=Data3_SCODE_matrix(:,end-1:end);

z3A=nan(maxsize,2);
z3A(1:size(Data3_Area_matrix,1),:)=Data3_Area_matrix(:,end-1:end);

z3O=nan(maxsize,2);
z3O(1:size(Data3_Orig_matrix,1),:)=Data3_Orig_matrix(:,end-1:end);
% figure
% boxplot(zS2(:,4),'Labels',{['SCODE M=50 T=' num2str(nanmean(zS2(:,5)),'%10.1f') 's']})
%  h = findobj(gca,'Tag','Box');
% patch(get(h(1),'XData'),get(h(1),'YData'),'y','FaceAlpha',.2);
% hLegend = legend(findall(gca,'Tag','Box'), {['Group A ' num2str(nanmean(zS2(:,5)),'%10.1f') 's']});

%%
figure
subplot(1,2,1)
boxplot([z2S(:,1),z2O(:,1),z2A(:,1)],'Labels',...
    {['SCODE'],...
    ['GRISLI with Original TIGRESS'],...% num2str(nanmean(z2O(:,2)),'%10.1f') 's'],...
    ['GRISLI with Area TIGRESS']})
 h = findobj(gca,'Tag','Box');
 color_list=['y','b','g','y','b','g'];
 ylim([0.5 0.58])
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'), color_list(j),'FaceAlpha',.2);
 end
 
 hold on
 subplot(1,2,2)

boxplot([z3S(:,1),z3O(:,1),z3A(:,1)],'Labels',...
    {['SCODE'],...
    ['GRISLI with Original TIGRESS'],...% num2str(nanmean(z2O(:,2)),'%10.1f') 's'],...
    ['GRISLI with Area TIGRESS']})
 h = findobj(gca,'Tag','Box');
 color_list=['y','b','g','y','b','g'];
 ylim([0.5 0.58])
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'), color_list(j),'FaceAlpha',.2);
 end
 
 %% Boxplot with varying R
 

RArea_Results_matrix = dlmread('RGRISLI_area_perf.txt','\t', 631,0);
R_array=[10, 50,80,100, 200,500, 800,1000,1500,2000,3000, 4000,5000];

% maxsize=max([size(Data2_Area_matrix_100,1),size(Data2_Area_matrix_500,1),size(RArea_Results_matrix ,1)]);
mat_to_plot=nan(20,length(R_array));

for i=1:length(R_array)
    data=RArea_Results_matrix(RArea_Results_matrix(:,2)==R_array(i),5);
    mat_to_plot(1:size(data,1),i)=data;    
end    
% mat_to_plot(1:size(Data2_Area_matrix_100,1),length(R_array)+1)=Data2_Area_matrix_100(487:746,4);
% mat_to_plot(1:size(Data2_Area_matrix_100,1),length(R_array)+1)=Data2_Area_matrix_100(:,4);
% mat_to_plot(1:size(Data2_Area_matrix_500,1),length(R_array)+2)=Data2_Area_matrix_500(:,4);

%  idx = [1     2     3     10     4     11  5     6     7     8     9   ];
% mat_to_plot = mat_to_plot(:,idx);
mat_to_plot(mat_to_plot==0)=NaN;
figure
hold on
line(1:(length(R_array)),nanmean(mat_to_plot))
boxplot(mat_to_plot,'Labels',...
    {num2str(R_array')})
xlabel('R')
ylabel('AUROC')
ylim([0.5 0.58])


 h = findobj(gca,'Tag','Box');
 color_list=['y','b','g','y','b','g'];
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'), 'y','FaceAlpha',.2);
 end
%% 
figure
box on
plot([373,758],[nanmean(z2A5(:,5)),nanmean(zA5(:,5))],'oy')
hold on
plot([373,758],[nanmean(z2S50(:,5)),nanmean(zS50(:,5))],'og')
xlabel('nNumber of cells')
ylabel('Runtime (s)')