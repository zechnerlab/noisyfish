function [m, s] = SolveConditionalMoments(c1, c2, m0, s0, t)

% analytical expression for the conditional first and second order RNA moments
expc2 = exp(-c2*t);
expc2sq = expc2.^2;

m = c1/c2-(c1*expc2)/c2+expc2*m0;
s = c1^2/c2^2 + c1/c2 + (c1^2*expc2sq)/c2^2 - (2*c1^2*expc2)/c2^2 ...
    - (c1*expc2)/c2 - expc2sq*m0 - (2*c1*expc2sq*m0)/c2 + expc2*m0 + (2*c1*expc2*m0)/c2 + expc2sq*s0;
end