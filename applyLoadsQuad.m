function [Q]=applyLoadsQuad(nodes,elem,nodLoads,Q, forceLoad)
%
% (c) Numerical Factory 2019
%
%------------------------------------------------------------------------
% Apply a load BC on a boundary of a quadrilateral meshed domain. 
% The behabiour is very similar to the applyConvecQuad funtion used to
% apply convection BC.
%-------------------------------------------------------------------------
numElem=size(elem,1);
[numNod,ndim]=size(nodes);
numCo=size(nodLoads,2);
if (numCo==1) 
    error('applyLoadQuad: It can''t manage a single node'); 
end
for k=1:numElem
  aux=[0,0,0,0];
  for j=1:4
    r=find(nodLoads==elem(k,j)); %find if j node of element k is in the boundary nodes list 
      if(~isempty(r)), aux(j)=1; end
  end
  %See the Triang version for the explanation of the binary number    
  number=aux*[1;2;4;8];
  % Here we codify both edges and corners on the boundaries
  switch (number) 
    case  3, ij=[1,2,0,0];  % edge 1:   aux=[1,1,0,0];
    case  6, ij=[2,3,0,0];  % edge 2:   aux=[0,1,1,0];
    case  7, ij=[1,2,2,3];  % corner 2: aux=[1,1,1,0];
    case  9, ij=[4,1,0,0];  % edge 4
    case 11, ij=[4,1,1,2];  % corner 1 
    case 12, ij=[3,4,0,0];  % edge 3
    case 13, ij=[3,4,4,1];  % corner 4
    case 14, ij=[2,3,3,4];  % corner 3
    case 15, error('applyLoadQuad: It can''t manage two corners'); % two corners 
    otherwise, ij=[0,0,0,0];
  end
  for j=1:2:3
   if ~ij(j)==0
    n1=elem(k,ij(j)); 
    n2=elem(k,ij(j+1));
    h=norm(nodes(n1,:)-nodes(n2,:)); 
    posit=[n1*ndim-1, n1*ndim];
    Q(posit)=Q(posit)+0.5*h*forceLoad; %same (x,y) force on both nodes
    posit=[n2*ndim-1, n2*ndim];   
    Q(posit)=Q(posit)+0.5*h*forceLoad; %same (x,y) force on both nodes
   end
  end
end
end
