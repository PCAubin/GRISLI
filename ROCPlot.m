function ROCPlot(FPR_array_multiple,TPR_array_multiple,array_labels)
% ROCPlot plot the ROC curve based on the couples of true positive and false
% positive rates(TPR,FPR) FPR_array_multiple (numlabels*G^2) and TPR_array_multiple
% (numlabels*G^2) based on the array_label (1*numlabels) such as 'TIGRESS
% orig', 'TIGRESS area' or 'SCODE'
% Pierre-Cyril Aubin-Frankowski, 2018

checksizes=(size(FPR_array_multiple,1)==size(TPR_array_multiple,1))&&(length(array_labels)==size(TPR_array_multiple,1));
if ~checksizes
   warning('Dimension mismatch of the arguments of ROCPlot')
end

numlabels=length(array_labels);
legendInfo=cell(1,numlabels);

figure
hold on
for i=1:numlabels
plot([0 FPR_array_multiple(i,1:(end-1))],[0 TPR_array_multiple(i,1:(end-1))])
AUROC=trapz(FPR_array_multiple(i,:),TPR_array_multiple(i,:));
legendInfo{i} = strcat(string(array_labels(i))," AUROC= ", num2str(AUROC));
end
box on
plot([0 1], [0 1])
axis tight
legend(legendInfo)
xlabel('FPR (fall-out)')
ylabel('TPR (sensitivity)')
end