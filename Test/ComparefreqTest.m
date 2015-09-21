function ComparefreqTest
% test for comparing numerical and ahalytical calculation of natural frequencies
rootpath = startupTest;

P = Parameters;

% cycle goes through all tests
% xc:    coordinate of the discontinuity
% L:     length of bar
% R1,R2: radius of each section
% E:     modulus of elasticity: Young's modulus for each section
% ro:    density for each section
X = P.X;
L = P.L;
S1 = P.S1;
S2 = P.S2;
xc = P.xc;
E = P.E;
ro = P.ro;
ModeCnt = P.ModeCnt;
PointCnt = P.PointCnt;

alp = AlpShapes

% N_tau: the time's points
N_tau = 120;

% dTaumax: the maximal time
% dTaumax: the array of deltatau
% xinf: the position of interface
dTaumax = funcdtm(X);
dtau = linspace(0, dTaumax, N_tau)';
xinf = zeros(N_tau, 1);
xinf = xc * ones(N_tau, 1) + dtau;
% the graph of the position of interface
figure(1)
plot(dtau, xinf, 'LineWidth', 2)
grid on;
xlabel('\Delta\tau [s]')
figurename = ['position of interface: xc = ', num2str(xc), ' m.'];
title(figurename);
ylabel('interface [m]')

% freq_analyt: natural frequencies from analytical solution  
lambda_analyt = zeros(N_tau, ModeCnt);
% freq_analyt = zeros(N_tau, ModeCnt);
% freq_XFEM: natural frequencies from numerical solution
lambda_xfem = zeros(N_tau, ModeCnt);
for n = 1:N_tau
  xct = xinf(n);  
  lambda_analyt(n,:) = AnalytMethod2(xct,L,S1,S2,ModeCnt);
%   [lambda_xfem, eigvec] = eig_val_xfem(xct,L,S1,S2,E,ro,ModeCnt,PointCnt);
  lambda_xfem(n,:) = eig_val_xfem(xct,L,S1,S2,E,ro,ModeCnt,PointCnt);
end
 freq_analyt = lambda_analyt*sqrt(E/ro)/(2*pi);
 freq_xfem = sqrt(lambda_xfem)/(2*pi);
% relative deviation between numerical and analytical results
df = freq_xfem - freq_analyt;
df_e = df./ freq_analyt * 100.0;
% for n = 1:ModeCnt
%   lambda_num(:,n) =  DerivLambda(1,n) * dtau + lambda_xfem(n);
% end

for n = 1:ModeCnt
   figure(n);
   plot(xinf, freq_xfem(:,n), xinf, freq_analyt(:,n));
   grid on;
%  axis([min(xinf) max(xinf) min([freq_xfem(:,n) freq_analyt(:,n)]) max([freq_xfem(:,n) freq_analyt(:,n)])]);
   title(strcat('Mode #', num2str(n)));
   xlabel('discontinuity position, xc [m]');
   ylabel('Natural frequency, f [Hz]');
end

for n = 1:ModeCnt
   figure(n+3);
   plot(xinf, df(:,n));
   grid on;
   axis([min(xinf) max(xinf) -inf inf]);
   title(strcat('Absolute error #', num2str(n)));
   xlabel('discontinuity position, xc [m]');
   ylabel('Natural frequency error, df [Hz]');
end

for n = 1:ModeCnt
   figure(n+6);
   plot(xinf, df_e(:,n));
   grid on;
   axis([min(xinf) max(xinf) -inf inf]);
   title(strcat('Relative error #', num2str(n)));
   xlabel('discontinuity position, xc [m]');
   ylabel('Natural frequency error, df [%]');
end
stopdown(rootpath);
end

function freq = AnalytMethod2(xc,L,S1,S2,ModeCnt)

global dl;
global l;
global k;

l = L;
Pn = ModeCnt; % number of eigenvalues

% secondary parameters
l1 = xc;
l2 = L-xc;

dl = abs(l1 - l2);

k = (S1 - S2)/(S1 + S2);

% error definition
Err = 10^6;
e = 2*pi/(L*Err);

freq = zeros(Pn,1);

for p=0:Pn-1
    x1 = pi*p/L;
    x2 = (p + 1)*pi/L;
    
    y1 = func(x1);
    y2 = func(x2);
    
    while true
       % different signum verification
       if (y1*y2)>0 % function doesn't change the sign
          disp('function does not change the sign'); 
       end
       
       % error calculation
       if abs(x2 - x1)<e
          freq(p+1) = x;
          break;
       end
          
       x = (x1 + x2)/2.0;
       y = func(x);
       
       if y == 0
          freq(p+1) = x;
          break
       end
       
       if y1*y<0
          x2 = x;
          y2 = y;
          continue;
       end
       
       if y2*y<0
          x1 = x;
          y1 = y;
          continue;
       end
       
    end
    
end
 %freq = freq*sqrt(E/ro)/(2*pi);
end

function y = func(lam)

global dl;
global l;
global k;
 y = cos(lam * l) + k*cos(lam*dl);
end


function lambda = eig_val_xfem(xc,L,S1,S2,E,ro,ModeCnt,PointCnt)

E1 = E; E2 = E;
ro1 = ro; ro2 = ro;
 
% numbEl: number of elements
% xle: length of each element
% addterm: two additional terms
addterm = 2;
numbEl = PointCnt;
nodeCoordinates = linspace(0,L,numbEl + 1);
xle = L/numbEl;

p = 0;
% minimal value and number of node enclosed to discontinuity
[val, id] = min(abs(nodeCoordinates - xc));
% correction of discontinuity position
if (val/xle)<(p/100)
xc = nodeCoordinates(id);
end;

% N: number of nodes
N = length(nodeCoordinates);

% stiffnessmatrix: stiffness matrix
% massmatrix: mass matrix
% totalnodes: the sum of number of nodes and additinal nodes
totalnodes = N+addterm;

stiff=zeros(totalnodes,totalnodes);
mass=zeros(totalnodes,totalnodes);

% fct_ls: calculation the value of level set function for each node
fct_ls = nodeCoordinates - xc;

% Ke1:  stiffness matrix for normal element for 1-st domain
% Me1:  mass matrix for normal element for 1-st domain
% A11d: product E1*S1/xle  for 1-st domain
% A21d: product ro1*S1*xle for 1-st domain
A11d = E1*S1;
A21d = ro1*S1;
Ke1 = (A11d/xle)*[1 -1; -1 1];
Me1 = (A21d*xle)*[1/3 1/6; 1/6 1/3];

% Ke2:  stiffness matrix for normal element for 2-d domain
% Me2:  mass matrix for normal element for 2-d domain
% A12d: product E2*S2/xle  for 2-d domain
% A22d: product ro2*S2*xle for 2-d domain
A12d = E2*S2;
A22d = ro2*S2;
Ke2 = (A12d/xle)*[1 -1; -1 1];
Me2 = (A22d*xle)*[1/3 1/6; 1/6 1/3];

% calculation of level set's product of two nodes
for in = 1:numbEl
    if fct_ls(in)*fct_ls(in+1)>0 
        if fct_ls(in)<0
            stiff(in:in+1,in:in+1) = stiff(in:in+1,in:in+1) + Ke1;
            mass(in:in+1,in:in+1) = mass(in:in+1,in:in+1) + Me1;
        end
        if fct_ls(in)>0
            stiff(in:in+1,in:in+1) = stiff(in:in+1,in:in+1) + Ke2;
            mass(in:in+1,in:in+1) = mass(in:in+1,in:in+1) + Me2;
        end;
    end;
    
    if fct_ls(in)*fct_ls(in+1)==0
        if fct_ls(in+1)==0
            stiff(in:in+1,in:in+1) = stiff(in:in+1,in:in+1) + Ke1;
            mass(in:in+1,in:in+1) = mass(in:in+1,in:in+1) + Me1;    
        end;
        if fct_ls(in)==0
            stiff(in:in+1,in:in+1) = stiff(in:in+1,in:in+1) + Ke2;
            mass(in:in+1,in:in+1) = mass(in:in+1,in:in+1) + Me2;
        end;
     continue;
    end;
    
    if  fct_ls(in)*fct_ls(in+1)<0
        phi1 = fct_ls(in);
        phi2 = fct_ls(in+1);
        [xmass, xk] = stiffness_and_mass_matrix(xle,phi1,phi2,A11d,A21d,A12d,A22d);
        
        % Assembling global stiffness matrix for bar
        stiff(in:in+1,in:in+1) = stiff(in:in+1,in:in+1) + xk(1:2,1:2);
        stiff(N+1:N+2,N+1:N+2) = stiff(N+1:N+2,N+1:N+2) + xk(3:4,3:4);
        stiff(in:in+1,N+1:N+2) = stiff(in:in+1,N+1:N+2) + xk(1:2,3:4);
        stiff(N+1:N+2,in:in+1) = stiff(N+1:N+2,in:in+1) + xk(3:4,1:2);
        
        % Assembling global mass matrix for each element of bar
        mass(in:in+1,in:in+1) = mass(in:in+1,in:in+1) + xmass(1:2,1:2);
        mass(N+1:N+2,N+1:N+2) = mass(N+1:N+2,N+1:N+2) + xmass(3:4,3:4);
        mass(in:in+1,N+1:N+2) = mass(in:in+1,N+1:N+2) + xmass(1:2,3:4);
        mass(N+1:N+2,in:in+1) = mass(N+1:N+2,in:in+1) + xmass(3:4,1:2);
    end
end

% calculation of matrices D and V 
% D:    diagonal matrix containing the eigenvalues
% V:    matrix whose columns are the corresponding right eigenvectors
% freq: natural frequences
D = eig(stiff(2:totalnodes, 2:totalnodes), mass(2:totalnodes, 2:totalnodes));
% [V, D] = eig(stiff(2:totalnodes, 2:totalnodes), mass(2:totalnodes, 2:totalnodes));
lambda = zeros(length(D),1);
% eigvec = zeros(length(V),P.ModeCnt);
% condnumber: condition number of V (matrix of eigenvectors) 
% condnumber=cond(V);
for in = 1:length(D)
% lambda(in) = sqrt(D(in,in))/(2*pi);
% lambda(in) = D(in,in);  
lambda(in) = D(in);  
end
lambda = sort(lambda);
lambda = lambda(1:ModeCnt);

% for in=1:P.ModeCnt
%    eigvec(:,in) = V(:,in);
% end
% eigvec = eigvec';
% eigvec = eigvec(:,1:length(V) - 2);
end

function [P, K] = stiffness_and_mass_matrix(xle,phi1,phi2,A11d,A21d,A12d,A22d)

% definition size of mass matrix and stiffness matrix
xmatrixmass = zeros(4,4);
xmatrixk = zeros(4,4);
npg=5;
[x,w] = cal_point_gauss(npg);
d = abs(phi1);
% numerical integration (0 - d)
for ig = 1:npg
    s = (d/2)*(x(ig)+1);
    N2 = s/xle;
    N1 = 1-N2;
    
    F2 = calculation_f2(xle,phi1,phi2,s);
    N = [N1, N2, N1*F2, N2*F2];
    xmatrixmass = xmatrixmass+A21d*(N'*N)*w(ig);
    
    dF2 = calculation_df2(xle,phi1,phi2,s);
    
    dN(1) = -1/xle;
    dN(2) = 1/xle;
    dN(3) = dN(1)*F2+N(1)*dF2;
    dN(4) = dN(2)*F2+N(2)*dF2;
    dN = [dN(1), dN(2), dN(3), dN(4)];
    xmatrixk = xmatrixk+A11d*(dN'*dN)*w(ig);
    
end
    xmatrixmass = xmatrixmass*(d/2);
    xmatrixk = xmatrixk*(d/2);

    K_minus = xmatrixk;
    P_minus = xmatrixmass;
    
    xmatrixmass = zeros(4,4);
    xmatrixk = zeros(4,4);
    
 % numerical integration (d - xle)
for ig = 1:npg
    s =((xle-d)/2)*(x(ig)+1)+d;
    N2 = s/xle;
    N1 = 1-N2;
    
    F2 = calculation_f2(xle,phi1,phi2,s);
    N = [N1, N2, N1*F2, N2*F2];
    xmatrixmass = xmatrixmass+A22d*(N'*N)*w(ig);
    
    dF2 = calculation_df2(xle,phi1,phi2,s);
    
    dN(1) = -1/xle;
    dN(2) = 1/xle;
    dN(3) = dN(1)*F2+N(1)*dF2;
    dN(4) = dN(2)*F2+N(2)*dF2;
    dN = [dN(1), dN(2), dN(3), dN(4)];
    xmatrixk = xmatrixk+A12d*(dN'*dN)*w(ig);
    
end

    xmatrixmass = xmatrixmass*((xle-d)/2);
    xmatrixk = xmatrixk*((xle-d)/2);   
    
    K_plus = xmatrixk;
    P_plus = xmatrixmass;
    
    K = K_minus + K_plus;
    P = P_minus + P_plus;
    
end

function [F2] = calculation_f2(xle,phi1,phi2,variable)
% definition the shape function
    x2 = variable/xle;
    x1 = 1-x2;
    
% calculation function F2
    F2 = abs(phi1)*x1+abs(phi2)*x2-abs(phi1*x1+phi2*x2); 
end

function [dF2] = calculation_df2(xle,phi1,phi2,variable)
% definition of shape function
    x2 = variable/xle;
    x1 = 1-x2;
    aux = phi1*x1+phi2*x2;
    dF2 = abs(phi1)*(-1/xle)+abs(phi2)*(1/xle);   
    
% calculation derivative of the function dF2   
 if aux<=0 
    dF2 = dF2+phi1*(-1/xle)+phi2*(1/xle); 
 else dF2 = dF2-(phi1*(-1/xle)+phi2*(1/xle));
 end     
end

function [x,w] = cal_point_gauss(n)    

 x = zeros(n,1)                        ;
 w = zeros(n,1)                        ;
        
 if n == 1                             ;
    x(1)=0.0                           ;
    w(1)=2.0                           ;
 end;

 if n == 2                             ;
    x(1) = -0.5773502691896257         ;
    x(2) =   0.5773502691896257        ;
    w(1) =   1.0                       ;
    w(2) =   1.0                       ;
 end;

 if n == 3    
    x(1) = -0.7745966692414834         ;
    x(2) =   0.0                       ;
    x(3) =   0.7745966692414834        ;
    w(1) =   0.5555555555555552        ;
    w(2) =   0.8888888888888888        ;
    w(3) =   0.5555555555555552        ;   
 end;

 if n == 4                             ;
    x(1) = -0.8611363115940526         ;
    x(2) = -0.3399810435848563         ;
    x(3) =   0.3399810435848563        ;
    x(4) =   0.8611363115940526        ;
    w(1) =   0.3478548451374537        ;
    w(2) =   0.6521451548625464        ;
    w(3) =   0.6521451548625464        ;
    w(4) =   0.3478548451374537        ;    
 end;
       
 if n == 5                             ;    
    x(1) = -0.9061798459386640         ;
    x(2) = -0.5384693101056831         ;
    x(3) =   0.0                       ;
    x(4) =   0.5384693101056831        ;
    x(5) =   0.9061798459386640        ;
    w(1) =   0.2369268850561890        ;
    w(2) =   0.4786286704993665        ;
    w(3) =   0.5688888888888889        ;
    w(4) =   0.4786286704993665        ;
    w(5) =   0.2369268850561890        ;    
 end;

 if n == 6                             ;    
    x(1) = -0.9324695142031521         ;
    x(2) = -0.6612093864662646         ; 
    x(3) = -0.2386191860831969         ;
    x(4) =   0.2386191860831969        ;
    x(5) =   0.6612093864662646        ;        
    x(6) =   0.9324695142031521        ;
    w(1) =   0.1713244923791705        ;
    w(2) =   0.3607615730481386        ;
    w(3) =   0.4679139345726913        ; 
    w(4) =   0.4679139345726913        ; 
    w(5) =   0.3607615730481386        ;
    w(6) =   0.1713244923791705        ;   
 end;
   
 if n == 7                             ;    
    x(1) = -0.9491079123427585         ;
    x(2) = -0.7415311855993945         ;
    x(3) = -0.4058451513773972         ;
    x(4) =   0.0                       ;
    x(5) =   0.4058451513773972        ;
    x(6) =   0.7415311855993945        ;
    x(7) =   0.9491079123427585        ;
    w(1) =   0.1294849661688697        ;
    w(2) =   0.2797053914892767        ;
    w(3) =   0.3818300505051190        ;
    w(4) =   0.4179591836734694        ;
    w(5) =   0.3818300505051190        ;
    w(6) =   0.2797053914892767        ;
    w(7) =   0.1294849661688697        ;    
 end;

 if n == 8;    
    x(1) = -0.9602898564975363         ;
    x(2) = -0.7966664774136268         ;
    x(3) = -0.5255324099163290         ;
    x(4) = -0.1834346424956498         ;
    x(5) =   0.1834346424956498        ;
    x(6) =   0.5255324099163290        ;
    x(7) =   0.7966664774136268        ;
    x(8) =   0.9602898564975363        ;
    w(1) =   0.1012285362903768        ;
    w(2) =   0.2223810344533745        ;
    w(3) =   0.3137066458778874        ;
    w(4) =   0.3626837833783620        ;
    w(5) =   0.3626837833783620        ;
    w(6) =   0.3137066458778874        ;
    w(7) =   0.2223810344533745        ;
    w(8) =   0.1012285362903768        ;    
 end;

 if n == 9                             ;    
    x(1) = -0.9681602395076261         ;
    x(2) = -0.8360311073266359         ;
    x(3) = -0.6133714327005905         ;
    x(4) = -0.3242534234038089         ;
    x(5) = 0.0                         ;
    x(6) =   0.3242534234038089        ;
    x(7) =   0.6133714327005905        ;
    x(8) =   0.8360311073266359        ;
    x(9) =   0.9681602395076261        ;
    w(1) =   0.08127438836157463       ;
    w(2) =   0.1806481606948576        ;
    w(3) =   0.2606106964029355        ;
    w(4) =   0.3123470770400029        ;
    w(5) =   0.3302393550012598        ;
    w(6) =   0.3123470770400029        ;
    w(7) =   0.2606106964029355        ;
    w(8) =   0.1806481606948576        ;
    w(9) =   0.08127438836157463       ;
 end;

 if n == 10                            ;    
    x(1)   = -0.9739065285171716       ;
    x(2)   = -0.8650633666889845       ;
    x(3)   = -0.6794095682990244       ;
    x(4)   = -0.4333953941292472       ;
    x(5)   = -0.1488743389816312       ;
    x(6)   =   0.1488743389816312      ;
    x(7)   =   0.4333953941292472      ;
    x(8)   =   0.6794095682990244      ;
    x(9)   =   0.8650633666889845      ;
    x(10) =   0.9739065285171716       ;
    w(1)   =   0.06667134430868774     ;
    w(2)   =   0.1494513491505805      ;
    w(3)   =   0.2190863625159822      ;
    w(4)   =   0.2692667193099962      ;
    w(5)   =   0.2955242247147529      ;
    w(6)   =   0.2955242247147529      ;
    w(7)   =   0.2692667193099962      ;
    w(8)   =   0.2190863625159822      ;
    w(9)   =   0.1494513491505805      ;
    w(10) =   0.06667134430868774      ;    
 end;

 if n == 11                            ;
    x(1)   = -0.9782286581460570       ;
    x(2)   = -0.8870625997680953       ;
    x(3)   = -0.7301520055740494       ;
    x(4)   = -0.5190961292068118       ;
    x(5)   = -0.2695431559523450       ;
    x(6)   = 0.000000000000000         ;
    x(7)   =   0.2695431559523450      ;
    x(8)   =   0.5190961292068118      ;
    x(9)   =   0.7301520055740494      ;
    x(10) =   0.8870625997680953       ;
    x(11) =   0.9782286581460570       ;
    w(1)   =   5.566856711617354e-02   ; 
    w(2)   =   1.255803694649047e-01   ; 
    w(3)   =   1.862902109277343e-01   ; 
    w(4)   =   2.331937645919903e-01   ; 
    w(5)   =   2.628045445102466e-01   ;
    w(6)   =   2.729250867779006e-01   ;
    w(7)   =   2.628045445102466e-01   ;
    w(8)   =   2.331937645919903e-01   ;
    w(9)   =   1.862902109277343e-01   ;
    w(10) =   1.255803694649047e-01    ;
    w(11) =   5.566856711617354e-02    ;
 end;

 if n == 12                            ;   
    x(1)   = -9.815606342467192e-01    ;
    x(2)   = -9.041172563704748e-01    ;
    x(3)   = -7.699026741943047e-01    ;
    x(4)   = -5.873179542866175e-01    ;
    x(5)   = -3.678314989981802e-01    ;
    x(6)   = -1.252334085114689e-01    ;
    x(7)   =   1.252334085114689e-01   ;
    x(8)   =   3.678314989981802e-01   ;
    x(9)   =   5.873179542866175e-01   ;
    x(10) =   7.699026741943047e-01    ;
    x(11) =   9.041172563704748e-01    ;
    x(12) =   9.815606342467192e-01    ;
    w(1)   =   4.717533638651183e-02   ;
    w(2)   =   1.069393259953182e-01   ;
    w(3)   =   1.600783285433463e-01   ;
    w(4)   =   2.031674267230658e-01   ;
    w(5)   =   2.334925365383548e-01   ;
    w(6)   =   2.491470458134029e-01   ;
    w(7)   =   2.491470458134029e-01   ;
    w(8)   =   2.334925365383548e-01   ;
    w(9)   =   2.031674267230658e-01   ;
    w(10) =   1.600783285433463e-01    ;
    w(11) =   1.069393259953182e-01    ;
    w(12) =   4.717533638651183e-02    ;
 end;

 if n == 13                            ;
    x(1)   = -9.841830547185881e-01    ;
    x(2)   = -9.175983992229780e-01    ;
    x(3)   = -8.015780907333099e-01    ;
    x(4)   = -6.423493394403402e-01    ;
    x(5)   = -4.484927510364468e-01    ;
    x(6)   = -2.304583159551348e-01    ;
    x(7)   =   0.000000000000000e-00   ;
    x(8)   =   2.304583159551348e-01   ;
    x(9)   =   4.484927510364468e-01   ;
    x(10) =   6.423493394403402e-01    ;
    x(11) =   8.015780907333099e-01    ;
    x(12) =   9.175983992229780e-01    ;
    x(13) =   9.841830547185881e-01    ;
    w(1)   =   4.048400476531581e-02   ;
    w(2)   =   9.212149983772838e-02   ;
    w(3)   =   1.388735102197872e-01   ;
    w(4)   =   1.781459807619457e-01   ;
    w(5)   =   2.078160475368884e-01   ;
    w(6)   =   2.262831802628971e-01   ;
    w(7)   =   2.325515532308739e-01   ;
    w(8)   =   2.262831802628971e-01   ;
    w(9)   =   2.078160475368884e-01   ;
    w(10) =   1.781459807619457e-01    ;
    w(11) =   1.388735102197872e-01    ;
    w(12) =   9.212149983772838e-02    ;
    w(13) =   4.048400476531581e-02    ;
 end;
end