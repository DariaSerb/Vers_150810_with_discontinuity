classdef Parameters
% system parameteters value
    properties
        L = 2;
        R1 = 0.02;
        R2 = 0.01;
        E = 2.06e11;
        ro = 7850;
        ModeCnt = 16;
        ModeEst = 3;
        PointCnt = 200;
%       Pointdisc = 401;
        Pointdisc = 1;
                xc = 0.9854;
        T = 10;
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
%         function out = xc(obj)
%             out = linspace(0, obj.L, obj.Pointdisc);
%         end
    end
end