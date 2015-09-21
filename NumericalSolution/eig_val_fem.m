function lambda = eig_val_fem(P)
L = P.L;
S = P.S;
E = P.E;
ro = P.ro;
% definition initial values
% numbEl: number of elements
% xle: length of each element
numbEl=400;
nodeCoordinates=linspace(0,L,numbEl+1);
xle=L/numbEl;

% numberNodes: number of nodes
numberNodes=size(nodeCoordinates,2);

% stiffnessmatrix: stiffness matrix
% massmatrix: mass matrix
stiffnessmatrix=zeros(numberNodes,numberNodes);
massmatrix=zeros(numberNodes,numberNodes);
                 
% Ke:  stiffness matrix 
% Me:  mass matrix 
% A1d: product E1*S1/xle  
% A2d: product ro1*S1*xle 
A1d=E*S/xle;
A2d=ro*S*xle;
Ke=A1d*[1 -1; -1 1];
Me=A2d*[1/3 1/6; 1/6 1/3];
 
for in=1:numbEl
    stiffnessmatrix(in:in+1,in:in+1)=stiffnessmatrix(in:in+1,in:in+1)+Ke;
    massmatrix(in:in+1,in:in+1)=massmatrix(in:in+1,in:in+1)+Me;
end
% calculation of matrices D and V 
% D:    diagonal matrix containing the eigenvalues
% V:    matrix whose columns are the corresponding right eigenvectors
D=eig(stiffnessmatrix(2:numberNodes,2:numberNodes),massmatrix(2:numberNodes,2:numberNodes));
lambda = zeros(length(D),1);
for in=1:length(D)
 lambda(in) = D(in);
end

lambda = sort(lambda);
lambda = lambda(1:P.ModeCnt);
end