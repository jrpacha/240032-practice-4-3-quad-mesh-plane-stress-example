%% Plane Elasticity: Quadrangular Meshes
% *Example: Constant rigth traction*
% 
% The $5\times 2\times 0.05,\mathrm{mm}$ plate shown in the 
% Figure has hole of $1\, \mathrm{mm}$ dimeter, centered 
% $1\, \mathrm{mm}$ right and $1\, \mathrm{mm}$ up from its 
% lower left corner. The plate is clamped along the hole, 
% whereas its right edge is under a constant traction force 
% of $\tau = 5000\, \mathrm{N/mm}.$ 
% 
% <<placaAmbForatCrop.png>>
%
% Apply the FEM, using the mesh in |meshPlacaForatQuad.m| to 
% compute the displacements of the nodes, the stress-strain and
% the von Mises stress. Material properties: Young modulus
% $E = 10^{7}\, \mathrm{N/mm}^{2}$, Poisson ratio
% $\nu = 0.3$.
%
% See Toni Susin's Numerical Factory,
% <https://numfactory.upc.edu/web/FiniteElements/Pract/P10-PlaneElasticity/html/PlaneElasticityQuadMesh.html FEM course, practice 4.3>
% 

clearvars
close all

excelFileDispl='displacements.xlsx';
excelFileStress = 'stress.xlsx';

%Data
E=1.0e+7;                  %Young Modulus (N/mm^2)
mu=0.3;                    %Poisson's ratio (adimensional)
th=0.05;                   %thickness (in mm)
forceLoad=[5.0e+3; 0.0e0]; %[Fx,Fy] traction force (in N/mm)

modelProblem='stress';

eval('meshPlacaForatQuad');
[numNod,ndim]=size(nodes);
numElem=size(elem,1);
numbering=0;
plotElements(nodes, elem, numbering)
hold on

%%%
% Select boundary points
indRight=find(nodes(:,1)>4.99);
indCirc=find(sqrt((nodes(:,1)-1.0e0).^2 + (nodes(:,2)-1.0e0).^2) < 0.51);
 
plot(nodes(indRight,1),nodes(indRight,2),'ok','lineWidth',1,...
    'markerFaceColor','green','markerSize',4,...
    'markerEdgeColor','green')
plot(nodes(indCirc,1),nodes(indCirc,2),'ok','lineWidth',1,...
    'markerFaceColor','red','markerSize',4,...
    'markerEdgeColor','red')
hold off

%% Material properties
% * |modelProblem == 'srtess'|: plane stress problem (thickness must be specified)
% * |modelProblem == 'strain'|: plane strain problem (thickness set to 1)

switch lower(modelProblem)
    case 'stress'
        c11=E/(1-mu^2);
        c22=c11;
        c12=mu*c11;
        c21=c12;
        c33=E/(2*(1+mu));   
    case 'strain'
        th=1;
        c11=E*(1-mu)/((1+mu)*(1-2*mu)); 
        c22=c11;
        c12=c11*mu/(1-mu);
        c21=c12;
        c33=E/(2*(1+mu));
    otherwise
        error('Please, specify whether it''s a ''stress'' or a ''strain'' problem!')
end

C = [c11, c12, 0; c21, c22, 0; 0, 0, c33];

%Computation of the stiffness matrix
K=zeros(ndim*numNod);
F=zeros(ndim*numNod,1);
Q=zeros(ndim*numNod,1);
for e=1:numElem
    Ke=planeElastQuadStiffMatrix(nodes,elem,e,C,th);
    %
    % Assemble the stiffness matrices
    %
    row=[2*elem(e,1)-1; 2*elem(e,1); 2*elem(e,2)-1; 2*elem(e,2); 
         2*elem(e,3)-1; 2*elem(e,3); 2*elem(e,4)-1; 2*elem(e,4)];
    col=row;
    K(row,col)=K(row,col)+Ke;
end

%% Natural B.C.
% constant traction  on the right edge
nodLoads=indRight'; %nodes the traction is applied at
Q=applyLoadsQuad(nodes,elem,nodLoads,Q,forceLoad);

%% Essential B.C.
% set displacements along the hole to zero
fixedNodes=[ndim*indCirc'-1; ndim*indCirc'];
freeNodes=setdiff(1:ndim*numNod,fixedNodes);
u=zeros(ndim*numNod,1); %initialize the solution to u=0
u(fixedNodes)=0.0;

%Reduced system
%Remark: the linear system is not modified. This is only valid if the BC=0
Km=K(freeNodes,freeNodes);
Qm=Q(freeNodes);

%solve the reduced system
um=Km\Qm;
u(freeNodes)=um;

%% Post process 
% Plot displacements stress-strain and von Mises stress
scale=7;
plotPlaneNodElemDespl(nodes, elem, u, scale);

%Contour Solution: Desp.X
valueToShow=u(1:2:end,:);
titol='Desp. X';
colorTable='jet';
plotContourSolution(nodes,elem,valueToShow,titol,colorTable);

%Compute Strainand StressElement
[stress,vonMises]=computeQuadStrainStressVM(nodes,elem,u,C);

%Plot VM stress
valueToShow=vonMises;
titol='VonMises Stress';
colorTable='jet';
plotContourSolution(nodes,elem,valueToShow,titol,colorTable);

%% Post process
% Write tables with displacements and stress
format short e
displTable=table(int64(1:numNod)',nodes(:,1),nodes(:,2),...
                  u(1:2:end),u(2:2:end),...
                  'VariableNames',{'NumNod','X','Y','UX','UY'});

%Write table to an Excel file                      
writetable(displTable,excelFileDispl);
%And show the displacements just for the last 10 nodes
displTable(end-9:end,:)

% Table with the stress tensor SX, SY, SXY and von Misses
format short e
stressTable=table(int64(1:numNod)',stress(:,1),stress(:,2),stress(:,3),...
                  vonMises',...
                  'VariableNames',{'NumElem','SX','SY','SXY','VM'});
%Write table to an Excel file                      
writetable(stressTable,excelFileStress);                            
%And show the stress just for the last 10 results
stressTable(end-9:end,:)
