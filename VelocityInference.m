function [V,Vf,Vp]=VelocityInference(X,Knorm)
% VelocityInference computes the midflux matrix (V), the influx matrix (Vp)
% and the outflux matrix (Vf)
% It takes as inputs the spacetime (time in first column) matrix X (Cx(1+G))
% and the matrix Knorm (GxG) of the weights summing to one in past and future 
% Pierre-Cyril Aubin-Frankowski, 2018
Dt = bsxfun(@minus,X(:,1),X(:,1)');
invDt=Dt.^-1;
invDt(isinf(invDt)) = 0;
if size(X,1)<2000
Dx = bsxfun(@minus,reshape(X(:,2:end),size(X,1),1,size(X,2)-1),reshape(X(:,2:end),1,size(X,1),size(X,2)-1));
Vf=squeeze(sum(invDt.*triu(Knorm).*Dx,2));
Vp=squeeze(sum(invDt.*tril(Knorm).*Dx,2));
V=(Vf+Vp)/2;
else
    tic
    V=zeros(size(X,1),size(X,2)-1);
    suKnorm=-triu(invDt.*Knorm);
    slKnorm=-tril(invDt.*Knorm);
    for i=1:size(X,1)
        Vf=suKnorm(i,:)*(X(:,2:end)-X(i,2:end));
        Vp=slKnorm(i,:)*(X(:,2:end)-X(i,2:end));
        V(i,:)=(Vf+Vp)/2;
    end
    elapsedTime=toc;
    disp(['Time to compute velocities for large matrices =' num2str(elapsedTime) 's'])
end

end