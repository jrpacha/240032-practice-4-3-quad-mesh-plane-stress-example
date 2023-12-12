clearvars
close all

%Data
E=1.0e+8;                  %Young Modulus (N/m^2)
nu=0.25;                   %Poisson's ratio (adimensional)
th=0.01;                   %thickness (m)
forceLoad=[530.0; 0.0e0];  %(Fx,Fy) traction force (in N/mm)

%Define the plane elasticity problem
modelProblem='stress'; %Plane stress

nodes = [
    -2, -2;
    2, -2;
    2, 2;
    -2, 2;
    -1, -1;
    1, -1;
    1, 1;
    -1, 1
    ];

elem = [
    1, 2, 6, 5;
    2, 3, 7, 6;
    3, 4, 8, 7;
    4, 1, 5, 8
    ];

[numNodes,ndim]=size(nodes);
numElem=size(elem,1);

Fx = forceLoad(1);
Fy = forceLoad(2);
L23 = norm(nodes(2,:)-nodes(3,:));

%Plot elements and constant traction
numbering=1;
plotElementsOld(nodes, elem, numbering)

v2 = nodes(2,:);
v3 = nodes(3,:);
h = norm(nodes(2,:)-nodes(3,:));
Fx = forceLoad(1);
scale = 0.5;
ndiv = 10;
signe = 1;
plotEdgeConstantBC(v2,v3,Fx,ndiv,scale,signe);

fprintf('Prob 4. Part A:\n')

switch modelProblem
    case 'stress'
        c11=E/(1-nu^2);
        c22=c11;
        c12=nu*c11;
        c21=c12;
        c33=E/(2*(1+nu));
        fprintf('Plane stress problem\n')
    case 'strain'
        th=1.0;
        c11=E*(1-nu)/((1+nu)*(1-2*nu));
        c22=c11;
        c12=c11*nu/(1-nu);
        c21=c12;
        c33=E/(2*(1+nu));
        fprintf('Plane strain problem\n')
    otherwise
        error('modelProblem should be stress or strain');
end
C=[c11, c12, 0; c21, c22, 0; 0, 0, c33];

%Computation of the stiffness matrix
K=zeros(ndim*numNodes);
F=zeros(ndim*numNodes,1);
Q=zeros(ndim*numNodes,1);
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

%% Boundary Conditions
indLeft = [1;4];
fixedNods = [ndim*indLeft'-1, ndim*indLeft'];
freeNods = setdiff(1:ndim*numNodes,fixedNods);
 
%Natural B.C.: horizontal constant traction [580;0] N/m is applied on the
%right external boundary

%Q=applyLoadsQuad(nodes,elem,nodLoads,Q,forceLoad);

%nod 2;
inod = 2;
Q(ndim*inod-1) = L23*Fx/2;
Q(ndim*inod) = L23*Fy/2;
%nod 3
inod = 3;
Q(ndim*inod-1) = L23*Fx/2;
Q(ndim*inod) = L23*Fy/2;

% Essential B.C.: 
% set displacements along the hole to zero
u=zeros(ndim*numNodes,1); %initialize the solution to u=0
u(fixedNods)=0.0;

%Reduced system
%Remark: the linear system is not modified. This is only valid if the BC=0
Km=K(freeNods,freeNods);
Qm=Q(freeNods);

%solve the system
um=Km\Qm;
u(freeNods)=um;

fprintf('max UY = %.4e\n',max(abs(u(2:2:end))))
fprintf('Hint. UX = %.4e\n',max(abs(u(1:2:end))))

%Graphical output (displacements and stress)
esc=30; %scale for the displacements
plotPlaneNodElemDespl(nodes, elem, u, esc);