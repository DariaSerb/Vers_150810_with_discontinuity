function plotTest()
% test for plot of numerical and exact solutions
rootpath = startupTest;

dTau = 0.05;
[x,u_num] = EigenFunction(dTau);

P = Parameters;
X = P.X;
u_exact = EigenFunctionSectionBarExact(X);

figure(1);
plot(x, u_num, X, u_exact);
legend('Uapprox1','Uapprox2','Uapprox3','Uexact1','Uexact2','Uexact3',3)
grid on;
title(['plot of numerical and exact solutions \Delta\tau =',num2str(dTau),'(s)']);
xlabel('x');
ylabel('u');
% ifig = 0    ;
% for im = 1 : P.ModeCnt
%     ifig = ifig + 1  ;
%     figure(ifig)
%     plot(u_exact(im,:),u_num(im,:),'LineWidth',2)
%     grid
% end
stopdown(rootpath);
end

