function Knorm = KernelNormalization(K)
% KernelNormalization normalizes the rows of the weight matrix K (GxG) so
% that the subrow on the left of the diagonal (past) and the rows on its right
% (future) sum to one, to make it a convex combination. Knorm is of size
% (GxG).
% Pierre-Cyril Aubin-Frankowski, 2018

Kf = triu(K);
normrow=sum(Kf,2).^-1;
normrow(isinf(normrow)) = 0;
Kfnorm=Kf.*normrow;
Kp = tril(K);
normrow=sum(Kp,2).^-1;
normrow(isinf(normrow)) = 0;
Kpnorm=Kp.*normrow;
Knorm=Kfnorm+Kpnorm;
end