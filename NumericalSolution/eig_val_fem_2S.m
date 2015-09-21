% function lambda = eig_val_fem_2S(P)

function lambda = eig_val_fem_2S()
% definition initial values
% L = P.L;
% S1 = P.S1;
% S2 = P.S2;
% E = P.E;
% ro = P.ro;
% xc = P.xc;
L = 2;
R1 = 0.02; R2 = 0.01; S1=pi*R1^2; S2=pi*R2^2;
E = 2.06e11; ro = 7850;
xc = 0.9854; P.ModeCnt = 3;

% numbl: number of elements of the first section
% numbm: number of elements of the second section
% xle: length of each element
numbEl = 5;
% numbEl = 200;
nodeCoordinates = linspace(0,L,numbEl+1);
xle = L/numbEl;

% numberNodes: number of nodes
numberNodes = size(nodeCoordinates,2);

p = 1;
% minimal value and number of node enclosed to discontinuity
[val, id] = min(abs(nodeCoordinates - xc));
% correction of discontinuity position
if (val/xle)<(p/100)
xc = nodeCoordinates(id);
end;

% stiffnessmatrix: stiffness matrix
% massmatrix: mass matrix
stiffnessmatrix = zeros(numberNodes,numberNodes);
massmatrix = zeros(numberNodes,numberNodes);

% A11, A21: product E1*S1/xle and ro1*S1*xle 
% A12, A22: product E1*S2/xle and ro1*S2*xle 
% Ke:  stiffness matrix 
% Me:  mass matrix 
A11 = E*S1/xle;  A12 = E*S2/xle; A21 = ro*S1*xle; A22 = ro*S2*xle;
Ke1 = A11*[1 -1; -1 1]; Ke2 = A12*[1 -1; -1 1];
Me1 = A21*[1/3 1/6; 1/6 1/3]; Me2 = A22*[1/3 1/6; 1/6 1/3];
 
for in=1:id-1
    stiffnessmatrix(in:in+1,in:in+1)=stiffnessmatrix(in:in+1,in:in+1)+Ke1;
    massmatrix(in:in+1,in:in+1)=massmatrix(in:in+1,in:in+1)+Me1;
end
for in=id:numbEl
    stiffnessmatrix(in:in+1,in:in+1)=stiffnessmatrix(in:in+1,in:in+1)+Ke2;
    massmatrix(in:in+1,in:in+1)=massmatrix(in:in+1,in:in+1)+Me2;
end

% calculation of matrices D and V 
% D:    diagonal matrix containing the eigenvalues
% V:    matrix whose columns are the corresponding right eigenvectors
D=eig(stiffnessmatrix(2:numberNodes,2:numberNodes),massmatrix(2:numberNodes,2:numberNodes));
lambda = zeros(length(D),1);

for in=1:length(D)
%  lambda(in) = D(in);
 lambda(in) = sqrt(D(in))/(2*pi);
end

lambda = sort(lambda);
lambda = lambda(1:P.ModeCnt);
end