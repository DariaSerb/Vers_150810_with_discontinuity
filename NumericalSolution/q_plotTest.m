function q_plotTest()
%plot graph q from X
L = 2.0;
N=1000;
X=linspace(0,L,N);

[q, qs] = q_calc(X);

figure(1)
plot(X,q);
grid on;
title('the dependence q(X)-function');
xlabel('X [m]');
ylabel('q');

figure(2)
plot (X,qs);
grid on;
title('differentiation of the q-function');
xlabel('X [m]');
ylabel('qs');
end