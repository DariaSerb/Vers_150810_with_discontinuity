function CompFreqTest
% test for comparing numerical and ahalytical calculation of natural frequencies
rootpath = startupTest;

P = Parameters;

% cycle goes through all tests
% xc:    coordinate of the discontinuity
% L:     length of bar
% R1,R2: radius of each section
% E:     modulus of elasticity: Young's modulus for each section
% ro:    density for each section

E = P.E;
% X = P.X;
ro = P.ro;
N_tau = P.Pointdisc;
ModeEst = P.ModeEst;
dTau = P.dTau;
dTaumax = P.dTaumax;
% xc = P.xc;
% xinf = P.xinf;

% dTaumax: the maximal time
% xinf: the position of interface
% dTaumax = funcdtm(X);

% the graph of the position of interface

% figure(1)
% plot(dTau, xinf, 'LineWidth', 2)
% grid on;
% xlabel('\Delta\tau [s]')
% figurename = ['position of interface: xc = ', num2str(xc), ' m.'];
% title(figurename);
% ylabel('interface [m]')

% xlswrite('testlvalueout.xlsx',dTau,'Feuil2','A2:A130');
% xlswrite('testlvalueout.xlsx',xinf,'Feuil2','B2:B130');

% freq_analyt: natural frequencies from analytical solution  
% lambda_analyt = zeros(N_tau, ModeEst);
  lambda_analyt = AnalytMethod(P);
  freq_analyt = lambda_analyt*sqrt(E/ro)/(2*pi);

% the definition eigenvalues using XFEM
[lambda, uf, V, stiffdq, massdq, mui] = eig_val(P);
  lambda = lambda(:,1:ModeEst);
  freq_XFEM = sqrt(lambda)/(2*pi);

% xlswrite('testlvalueout.xlsx',freq_analyt,'Feuil2','D2:F130');
% xlswrite('testlvalueout.xlsx',freq_XFEM,'Feuil2','H2:J130');

% relative deviation between numerical and analytical results
  df = abs(freq_XFEM - freq_analyt);
  dfe = df./ freq_analyt * 100.0;

% xlswrite('testlvalueout.xlsx',dfe,'Feuil2','L2:N130');  
  
% the definition Directional Derivatives of eigenvalues
[alp_num,DerivLambda,Du] = AlpShapes_num;

% DerivLambda(1,:) = [-3.185634e006 4.999783e007 -1.699867e008];

% xlswrite('testlvalueout.xlsx',DerivLambda,'Feuil2','A133:C133');

% freqEst: natural frequencies using Directional Derivatives
% lambda(1,n) is initiql approximation for eigenvalues 
LambdaEst = zeros(N_tau,ModeEst);
  for n = 1:ModeEst
   LambdaEst(:,n) = dTau*DerivLambda(1,n) + lambda(1,n);    
  end
freqEst = sqrt(LambdaEst)/(2*pi);

% xlswrite('testlvalueout.xlsx',freqEst,'Feuil2','P2:R130');

% relative deviation between numerical and analytical results
dfEst = abs(freqEst - freq_analyt);
dfeEst = dfEst./ freq_analyt * 100.0;

% xlswrite('testlvalueout.xlsx',dfeEst,'Feuil2','T2:V130');  

% -------------------------------------------------------------------------

for n = 1:ModeEst
   figure(n+1);
   plot(dTau, freq_analyt(:,n),':', dTau, freq_XFEM(:,n), dTau, freqEst(:,n), 'LineWidth', 2);
   legend('freqExact','freqXFEM','freqEst',3)
   grid on;
%  axis([min(xinf) max(xinf) min([freq_XFEM(:,n) freq_analyt(:,n)]) max([freq_XFEM(:,n) freq_analyt(:,n)])]);
   title(strcat('Mode # ', num2str(n), '  \Delta\tau_m_a_x = ', num2str(dTaumax), ' [s]'));
   xlabel('\Delta\tau   [s]');
   ylabel('Natural frequency, f [Hz]');
end

% for n = 1:ModeModeEst
%    figure(n+3);
%    plot(xinf, df(:,n));
%    grid on;
%    axis([min(xinf) max(xinf) -inf inf]);
%    title(strcat('Absolute error #', num2str(n)));
%    xlabel('discontinuity position, xc [m]');
%    ylabel('Natural frequency error, df [Hz]');
% end

for n = 1:ModeEst
   figure(n+4);
   plot(dTau, dfeEst(:,n), 'LineWidth', 3);
%    legend('freqExact','freqXFEM','freqEst',3)
   grid on;
%    axis([min(dTau) max(dTau) -inf inf]);
   title(strcat('Relative error #', num2str(n), '  \Delta\tau_m_a_x = ', num2str(dTaumax), ' [s]'));
   xlabel('\Delta\tau   [s]');
   ylabel('deviation frequency, df [%]');
end

stopdown(rootpath);
end

