function u = EigenFunctionSectionBarExact(X,xc)
% Eigenfunctions for 1D bar with two Sections 
% Eigenvalues for 1D bar with two Sections 
% input:
% X - array of discretization points
% xc - point-coordinate of the discontinuity
% output:
% u - array of eigenfunction values (several modes)

global P;

P = Parameters;
S1 = P.S1;
S2 = P.S2;
N_mode = P.ModeEst;
N_tau = P.Pointdisc;
N_x = length(X);
xinf = P.xinf;

for n = 1:N_tau
 p = 0;
% the definition of index of discontinuity position
 if abs(xinf(n) - xc) <= (p/1000)
  index = n;
 end;
end

u = zeros(N_mode, N_x);
W = AnalytMethod(P);
% fi = W*sqrt(P.E/P.ro)/(2*pi);

B1 = 1;
A2 = zeros(N_mode, 1);
B2 = zeros(N_mode, 1);

 for n_mode = 1:N_mode
     A2(n_mode) = B1*sin(W(index, n_mode)*xc);
     B2(n_mode) = (S1/S2)*B1*cos(W(index, n_mode)*xc);
    for n_x = 1:N_x
     u(n_mode, n_x) = u_calc(X(n_x), xc, W(index, n_mode), A2(n_mode), B2(n_mode));
    end
 end

end

function u = u_calc(x, xc, W, A2, B2)
%u_exact - exact solution of u in scalar form
% if x <= xc u = 1*sin(W*x); if x > xc u = cos(W*x) + sin(W*x);
B1 = 1;

 if x <= xc
   u  = B1*sin(W*x);    
   else u  = A2*cos(W*(x-xc)) + B2*sin(W*(x-xc));
 end
end


