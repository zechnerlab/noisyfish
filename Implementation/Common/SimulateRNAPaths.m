function [XS, VS, ZS] = SimulateRNAPaths_OneState(M, grid, tm, randFrac, cycleLength, a, b, v0T, c1, c2, m0s, s0s, splitFactor, sameAct)

if (nargin<14)
   sameAct = 0; 
end

T = max(grid);
numBins = 100;

if (sameAct)
    [tAct0, Act0] = SimulateActivationTime(T, a, b);
end

for k=1:M
    
    [Z, t] = SimulateCellCycle2(T, tm, randFrac, cycleLength);
    if (sameAct==0)
       [tAct, Act] = SimulateActivationTime(T, a, b); 
       %tAct = tAct2;
    else
       tAct = tAct0;
       Act = Act0;
    end

    gIdx = find(tAct(2)<t);
    if (isempty(gIdx))
       ActProfile = zeros(size(t)); 
    else
    gIdx = gIdx(1);
        t = [t(1:gIdx-1) tAct(2) t(gIdx:end)];
        Z = [Z(1:gIdx-1) Z(gIdx-1), Z(gIdx:end)];
        ActProfile = zeros(size(t));
        ActProfile(gIdx:end) = 1;
    end
    
    v0 = normrnd(v0T, 12000);
    
    
    m0 = m0s*v0;
    s0 = s0s*v0^2;

    
    
    X = m0;
    V = v0;
    tM = 0;
    
    [Pre, Post, c, X0] = CreateBirthDeathSystem();
    c(2) = c2;
    
    variance0 = s0 - m0^2;
    beta = m0 / variance0;
    alpha = m0*beta;
    X0(1) = floor(gamrnd(alpha, 1/beta));
    
    divideVolume = 1;
    
    for l=1:length(Z)-1
       
        activeRate = (1-Z(l))*ActProfile(l);
        

        currT = t(l);
        nextT = t(l+1);
        
        tGrid = linspace(0, nextT-currT, numBins);
        
        %activeRate = 1;
        
        c(1) = c1*activeRate*v0;
        
        
        [x, tx] = SimulateSSA(X0, c, Pre, Post, 10000, 0, max(tGrid));
        
        [xS] = SampleCTMPPathGrid_mex(x, tx, tGrid);
        
        %[m, s] = SolveConditionalMoments(c1*activeRate*v0, c2, m0, s0, tGrid);
        
        
        X = [X, xS(2:end)];
        V = [V, repmat(v0, 1, numBins-1)];
        tM = [tM, tGrid(2:end)+currT];
        
        %splitFactor = 1/2;
        if (Z(l)==1 && divideVolume<5)%if mitosis, split cell in two halfs until the 5th cell cycle
            factor = splitFactor;
            divideVolume = divideVolume + 1;
        else
            factor = 1;
        end
        X0 = x(end)*factor;
        
        v0 = v0*factor;
    end
    
    %ComputeConditionalRNASolution(t, Z);
    
    XTot = [X;V];
    XTotS = SampleCTMPPathGrid_mex(XTot, tM, grid);
    
    XS(k, :) = XTotS(1, :);
    VS(k, :) = XTotS(2, :);
    
    ZS(k, :) = SampleCTMPPathGrid_mex(Z, t, grid);
    
    %stairs(t/60, Z); hold on;
    %plot(tM/60, M); hold off; pause;
end