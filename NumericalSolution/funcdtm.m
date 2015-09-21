function deltatmax = funcdtm(X)
%deltatmax the definition maximal deltatime

[q, qs] = q_calc(X);
deltatmax = 1./max(abs(qs));

end

