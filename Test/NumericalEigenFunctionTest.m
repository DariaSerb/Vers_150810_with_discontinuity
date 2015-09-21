function NumericalEigenFunctionTest
% test for numerical calculation of eigen function
rootpath = startupTest;

dTau = 0.05;
[X,u] = EigenFunction(dTau);

figure(1);
plot(X, u);
legend('Uapprox1','Uapprox2','Uapprox3',3)
grid on;
title(['plot of numerical solution \Delta\tau =',num2str(dTau),'(s)']);
xlabel('x');
ylabel('u');

stopdown(rootpath);
end