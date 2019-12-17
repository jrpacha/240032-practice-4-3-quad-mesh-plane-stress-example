clearvars 
close all
eval('meshPlacaForatQuad')
vertexs=[0,0;5,0;5,2;0,2;0,0]
plot(vertexs(:,1),vertexs(:,2),'-','lineWidth',1,'color','black')
%axis off
axis tight
axis([0,7.5,-1, 3])
set(gca,'xtick',[0,1,2,3,4,5]);
set(gca,'ytick',[-1,0,1,2,3]);
set(gca,'box','off');
hold on
theta=linspace(0,2*pi,361);
radius=0.5;
xx=1+radius*cos(theta);
yy=1+radius*sin(theta);
plot(xx,yy,'k-','lineWidth',2)
forceLoad=[5000;0];
ndiv=10;
scale=1;
signe=1;
plotEdgeConstantBC(vertexs(2,:),vertexs(3,:),forceLoad,ndiv,scale,signe)
text(5.75,1,'$\tau = 5000\,\mathrm{N/mm}$','interpreter','LaTeX','fontSize',8)
saveas(gcf,'placaAmbForat.png');
hold off

