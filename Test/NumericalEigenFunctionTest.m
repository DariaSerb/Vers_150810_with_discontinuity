function NumericalEigenFunctionTest
% test for numerical calculation of eigenfunction
rootpath = startupTest;

dTau = 0.01;
% [X,u] = EigenFunction(dTau);
[X,u] = EigenFunction_num(dTau);

figure(1);
plot(X, u);
legend('Uapprox1','Uapprox2','Uapprox3',3)
grid on;
title(['plot of numerical solution \Delta\tau =',num2str(dTau),'(s)']);
xlabel('x');
ylabel('u');

stopdown(rootpath);
end
