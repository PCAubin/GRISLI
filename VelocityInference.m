function [V,Vf,Vp]=VelocityInference(X,Knorm)
% VelocityInference computes the midflux matrix (V), the influx matrix (Vp)
% and the outflux matrix (Vf)
% It takes as inputs the spacetime (time in first column) matrix X (Cx(1+G))
% and the matrix of the weights summing to one in past and future Knorm (GxG)
Dt = bsxfun(@minus,X(:,1),X(:,1)');
invDt=Dt.^-1;
invDt(isinf(invDt)) = 0;
Dx = bsxfun(@minus,reshape(X(:,2:end),size(X,1),1,size(X,2)-1),reshape(X(:,2:end),1,size(X,1),size(X,2)-1));
Vf=squeeze(sum(invDt.*triu(Knorm).*Dx,2));
Vp=squeeze(sum(invDt.*tril(Knorm).*Dx,2));
V=(Vf+Vp)/2;
end