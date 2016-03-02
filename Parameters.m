classdef Parameters
% system parameteters value
    properties
        L         = 2;
        R1        = 0.02;
        R2        = 0.016;
        E         = 2.06e11;
        ro        = 7860;
        ModeCnt   = 20;
        ModeEst   = 3;
        PointCnt  = 200;
        Pointdisc = 129;
        xc        = 1.2;
        dTaumax   = 0.32;
    end
    
    methods
        function out = S1(obj)
            out = pi*obj.R1^2;
        end
        function out = S2(obj)
            out = pi*obj.R2^2;
        end
        function out = X(obj)
            out = linspace(0, obj.L, obj.PointCnt);
        end
        function out = dTau(obj)
            out = linspace(0, obj.dTaumax*0.9, obj.Pointdisc)';
        end
        function out = xinf(obj)
            out = obj.xc*ones(obj.Pointdisc, 1) + obj.dTau;
        end
    end
end
