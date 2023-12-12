clearvars
close all

% Now, suppose that no loads are applied, and just the own weight of the 
% piece is taken into account, Furthermore, the density of the material is
% such that the vector of local force loads can be thought of constant with 
% the same value, Fe = [0;-390;0;-390;0;-390;0;-390;0;-390] for all the 
% elements

%Data
E=1.0e+8;                  %Young Modulus (N/m^2)
nu=0.25;                   %Poisson's ratio (adimensional)
th=0.01;                   %thickness (m)
fe=-340;                   %constant local force

%Define the plane elasticity problem
modelProblem='stress'; %Plane stress

%eval('meshPlacaForatQuad');

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

%Plot elements
numbering=1;
plotElementsOld(nodes, elem, numbering)

fprintf('Prob 4. Part C:\n')

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

Fe = [0;fe;0;fe;0;fe;0;fe];

for e=1:numElem
    Ke=planeElastQuadStiffMatrix(nodes,elem,e,C,th);
    %
    % Assemble the stiffness matrices
    %
    row=[2*elem(e,1)-1; 2*elem(e,1); 2*elem(e,2)-1; 2*elem(e,2); 
         2*elem(e,3)-1; 2*elem(e,3); 2*elem(e,4)-1; 2*elem(e,4)];
    col=row;
    K(row,col)=K(row,col)+Ke;
    F(row) = F(row) + Fe;
end

%% Boundary Conditions
indLeft = [1,4];
fixedNods = [ndim*indLeft-1, ndim*indLeft];
freeNods = setdiff(1:ndim*numNodes,fixedNods);

% Natural B.C.
% Not traction or compression forces are applied, so all the Q's on free 
% nodes are zero. These conditions are already set, since Q was initialised 
% to zero

% Essential B.C.: 
% set displacements along the hole to zero
u=zeros(ndim*numNodes,1); %initialize the solution to u=0
u(fixedNods)=0.0;

%Reduced system
%Remark: the linear system is not modified. This is only valid if the BC=0
Km=K(freeNods,freeNods);
Qm=Q(freeNods) + F(freeNods);

%solve the system
um=Km\Qm;
u(freeNods)=um;

fprintf('max UY = %.4e\n',max(abs(u(2:2:end))))
fprintf('Hint. UX = %.4e\n',max(abs(u(1:2:end))))

%Graphical output (displacements and stress)
esc=10; %scale for the displacements
plotPlaneNodElemDespl(nodes, elem, u, esc);