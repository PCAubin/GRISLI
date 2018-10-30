function VelocityPlot(X,V,Vf,Vp,scale)
% VelocityPlot computes (x,v) plots showing the measured points x with
% their midflux (and influx Vf and outflux Vp) as an arrow.
% It takes as inputs the spacetime (time in first column) matrix X (Cx(1+G))
% and flux matrices V (CxG)
% It plots the first component of the fluxes in the (t,x_1) plane (thus getting 2-D
% plots) as well a the point cloud projected to the (t,x_1) plane.
% Pierre-Cyril Aubin-Frankowski, 2018
N=size(X,1);
figure
if nargin<5
    scale=.2;
end
hold on
plot(X(:,1),X(:,2),'o');
hold on
quiver(X(:,1),X(:,2),ones(N,1),V,scale);

if nargin>3
    hold on
    quiver(X(:,1),X(:,2),-ones(N,1),-Vp,scale);
end
if nargin>2
    hold on
    quiver(X(:,1),X(:,2),ones(N,1),Vf,scale);
end
xlabel('t')
ylabel('x')
legend('measured point','midflux','influx','outflux')
ax = gca;
ax.FontSize=30;
end