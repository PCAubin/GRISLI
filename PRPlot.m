function PRPlot(TPR_array_multiple,PPV_array_multiple,array_labels,max_rec)
% PRPlot plot the PR curve based on the couples of true positive rate and precision (TPR,PPV)
% TPR_array_multiple (numlabels*G^2) and PPR_array_multiple
% (numlabels*G^2) based on the array_label (1*numlabels) such as 'TIGRESS
% orig', 'TIGRESS area' or 'SCODE'
% Pierre-Cyril Aubin-Frankowski, 2018
if nargin<4
    max_rec=1;
end
checksizes=(size(TPR_array_multiple,1)==size(PPV_array_multiple,1))&&(length(array_labels)==size(PPV_array_multiple,1));
if ~checksizes
   warning('Dimension mismatch of the arguments of PRPlot')
end

numlabels=length(array_labels);
legendInfo=cell(1,numlabels);

figure
hold on
for i=1:numlabels
plot([0 TPR_array_multiple(i,1:(end-1))],[0 PPV_array_multiple(i,1:(end-1))])
AUPR=trapz(TPR_array_multiple(i,:),PPV_array_multiple(i,:));
legendInfo{i} = strcat(string(array_labels(i))," AUPR= ", num2str(AUPR));
end
box on
% plot([0 max_rec], [0 1])
axis tight
axis([0 max_rec 0 1])
legend(legendInfo)
xlabel('TPR (recall)')
ylabel('PPV (precision)')
end