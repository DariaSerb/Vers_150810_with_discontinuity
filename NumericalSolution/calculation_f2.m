function [F2] = calculation_f2(xle,phi1,phi2,variable)
% definition the shape function
    x2 = variable/xle;
    x1 = 1-x2;
    
% calculation function F2
    F2 = abs(phi1)*x1+abs(phi2)*x2-abs(phi1*x1+phi2*x2); 
end