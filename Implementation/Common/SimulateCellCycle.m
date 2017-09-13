function [Z, t] = SimulateCellCycle(T, tm, randFrac, cycleLength)

maxN = length(cycleLength);

tIdx = 2;

K = 15;
t = zeros(1, K);
Z = zeros(1, K);

tauFun = @(n) cycleLength(min(n, maxN));

i = 1;

while t(tIdx-1)<T
    cellCycleLength = tauFun(i)-tm;
    
    tau0 = randFrac*cellCycleLength;
    ta = (1-randFrac)*cellCycleLength;
    
    Z(tIdx) = 0;
    t(tIdx) = t(tIdx-1)+ta;
    
    tIdx = tIdx + 1;
    
    
    tau = exprnd(tau0);
    Z(tIdx) = 1;
    t(tIdx) = t(tIdx-1) + tau;

    
    tIdx = tIdx + 1;
    
    
    Z(tIdx) = 0;
    t(tIdx) = t(tIdx-1) + tm;
    

    tIdx = tIdx + 1;
    
    i = i + 1;
    
end

t = t(1:tIdx-1);
Z = Z(1:tIdx-1);
