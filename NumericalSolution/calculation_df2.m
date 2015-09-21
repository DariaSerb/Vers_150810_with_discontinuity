function [dF2] = calculation_df2(xle,phi1,phi2,variable)
% definition of shape function
    x2 = variable/xle;
    x1 = 1-x2;
    aux = phi1*x1+phi2*x2;
    dF2 = abs(phi1)*(-1/xle)+abs(phi2)*(1/xle);   
    
% calculation derivative of the function dF2   
 if aux <= 0 
    dF2 = dF2+phi1*(-1/xle)+phi2*(1/xle); 
 else dF2 = dF2-(phi1*(-1/xle)+phi2*(1/xle));
 end     
end