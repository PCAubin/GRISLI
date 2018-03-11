function [K,sigx,sigt] = SpacetimeKernel(X,Alpha,sigx,sigt)
% SpacetimeKernel computes the weight matrix K (GxG).
% It takes as inputs the spacetime (time in first column) matrix X
% (Cx(1+G)), the Peanian kernel Alpha, and its parameters sigx and sigt (if
% those are not provided they are taken as the sqaure root of the tenth
% percentile of the distribution of distances in space and time).
Kx = squareform(pdist(X(:,2:end),'euclidean'));
Dt = bsxfun(@minus,X(:,1),X(:,1)');

if nargin<4
    sigt=sqrt(prctile(reshape(Dt.^2,1,[]), 10));
    if sigt==0
        sigt=mean2(Dt.^2);
    end
end
if nargin<3
    sigx=sqrt(prctile(reshape(Kx,1,[]), 10));
    if sigx==0
        sigx=mean2(Kx);
    end
end
if nargin<2
    Alpha=@(Kx,Dt,sigx,sigt)exp(-(Kx.^2)/(2*sigx^2)).*exp(-(Dt.^2)/(2*sigt^2)).*logical(Dt);
end
K=Alpha(Kx,Dt,sigx,sigt);

%check if every point gives a nonnull weight to at least one another point.
if(sum(any(K))~=size(K,1))
    disp('Min real value of Matlab reached, change parameters');
end
end
