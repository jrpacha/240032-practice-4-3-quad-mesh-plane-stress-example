clearvars
close all

%Data
E=1.0e+7;                  %Young Modulus (N/mm^2)
nu=0.3;                    %Poisson's ratio (adimensional)
th=0.05;                   %thickness (mm)
forceLoad=[5.0e+3; 0.0e0]; %(Fx,Fy) traction force (in n/mm)

%Define the plane elasticity problem: 
%modelProblem=1; %plane stress
%modelProblem=2; %plane strain (2)

modelProblem=1; %Plane stress

eval('meshPlacaForatQuad');
[numNod,ndim]=size(nodes);
numElem=size(elem,1);
numbering=0;
plotElements(nodes, elem, numbering)
hold on

%%Select Boundary points
indRight=find(nodes(:,1)>4.99);
indCirc=find(sqrt((nodes(:,1)-1.0e0).^2 + (nodes(:,2)-1.0e0).^2) < 0.51);
 
plot(nodes(indRight,1),nodes(indRight,2),'ok','lineWidth',1,...
    'markerFaceColor','green','markerSize',4)
plot(nodes(indCirc,1),nodes(indCirc,2),'ok','lineWidth',1,...
    'markerFaceColor','red','markerSize',4)
hold off

switch modelProblem
    case 1
        c11=E/(1-nu^2);
        c22=c11;
        c12=nu*c11;
        c21=c12;
        c33=E/(2*(1+nu));
        fprintf('Plane stress problem\n')
    case 2
        th=1.0;
        c11=E*(1-nu)/((1+nu)*(1-2*nu));
        c22=c11;
        c12=c11*nu/(1-nu);
        c21=c12;
        c33=E/(2*(1+nu));
        fprintf('Plane strain problem\n')
    otherwise
        error('modelProblem should be 1 (stress) or 2 (strain)');
end
C=[c11, c12, 0; c21, c22, 0; 0, 0, c33];

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

%%Boundary Conditions
% Natural B.C.: constant traction  on the right edge
nodLoads=indRight'; %nodes the traction is applied at
Q=applyLoadsQuad(nodes,elem,nodLoads,Q,forceLoad);

%% Essential B.C.: 
% set displacements along the hole to zero
fixedNodes=[ndim*indCirc'-1; ndim*indCirc'];
freeNodes=setdiff(1:ndim*numNod,fixedNodes);
u=zeros(ndim*numNod,1); %initialize the solution to u=0
u(fixedNodes)=0.0;

%Reduced system
%Remark: the linear system is not modified. This is only valid if the BC=0
Km=K(freeNodes,freeNodes);
Qm=Q(freeNodes);

%solve the system
um=Km\Qm;
u(freeNodes)=um;

%%Post Process
%Compute Strain and StressElement
[stress,vonMisses]=computeQuadStrainStressVM(nodes, elem, u, C);

%Output (fancy output: don't try this at the exams!)
%Displacements only the last 10 nodes:
displacements=[u(1:2:end),u(2:2:end)];
tableDispl=[(1:numNod)',nodes(:,1),nodes(:,2),displacements];
fprintf('\n%51s\n\n','Displacements (only for the last 10 nodes)')
fprintf('%7s%8s%12s%12s%11s\n',...
    'Num.Nod.','X','Y','U','V')
fprintf('%4d%16.4e%12.4e%12.4e%12.4e\n',tableDispl(end-9:end,:)');    
%Stress:
tableStress=[(1:numNod)',nodes(:,1),nodes(:,2),stress,vonMisses']; 
fprintf('\n%48s\n\n','Stress (only for the last 10 nodes)')
fprintf('%7s%8s%12s%13s%12s%12s%15s\n','Num.Nod.','X','Y',...
        'SXX','SYY','SXY','vonMisses')
fprintf('%4d%16.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n',...
        tableStress(end-9:end,:)')

%Graphical output (displacements and stress)
esc=10; %scale for the displacements
plotPlaneNodElemDespl(nodes, elem, u, esc);
valueToShow=abs(u(1:2:end)); %x displacement in each node
titol='Desp. X';
plotContourSolution(nodes,elem,valueToShow,titol,'Jet');
titol='VonMisses Stress';
plotContourSolution(nodes,elem,vonMisses,titol,'Jet');