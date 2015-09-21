function ExactEigenFunctionTest
% test for EigenFunctionSingleBarExact
rootpath = startupTest;

P = Parameters;
X = P.X;
u = EigenFunctionSectionBarExact(X);

figure(1);
plot(X,u);
legend('Uexact1','Uexact2','Uexact3',3)
grid on;
title(['plot of analytical solution u']);
xlabel('x');
ylabel('u');

stopdown(rootpath);
end