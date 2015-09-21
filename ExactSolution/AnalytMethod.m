function freq = AnalytMethod(P)

xc = P.xc;
L = P.L;
S1 = P.S1;
S2 = P.S2;
% E = P.E;
% ro = P.ro;

global dl;
global l;
global k;

l = L;
Pn = P.ModeCnt; % number of eigenvalues

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