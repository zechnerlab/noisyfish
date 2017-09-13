function [MS, SS, VS, ZS] = SimulateRNAMoments(M, grid, tm, randFrac, cycleLength, a, b, c1, c2, m0s, s0s, sameAct)

if (nargin<12)
   sameAct = 0; 
end

T = max(grid);
numBins = 40;

% Divide cell volume by a factor of #splitFactor upon cell division
splitFactor = 1/2; 

if (sameAct)
    [tAct0, ~] = SimulateActivationTime(T, a, b);
end

for k=1:M
    
    [Z, t] = SimulateCellCycle(T, tm, randFrac, cycleLength);
    if (sameAct==0)
       [tAct, ~] = SimulateActivationTime(T, a, b); 
    else
       tAct = tAct0;
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
    
    % draw random initial cell volume
    % values taken from experimental data (Stapel et al., Genes Dev (2017))
    v0 = normrnd(75000, 12000);
    
    m0 = m0s*v0;
    s0 = s0s*v0^2;

    S = s0;
    M = m0;
    V = v0;
    tM = 0;
    
    divideVolume = 1;
    
    for l=1:length(Z)-1
       
        activeRate = (1-Z(l))*ActProfile(l);
        

        currT = t(l);
        nextT = t(l+1);
        
        tGrid = linspace(0, nextT-currT, numBins);
        
        % Solve conditional Moments
        [m, s] = SolveConditionalMoments(c1*activeRate*v0, c2, m0, s0, tGrid);
        
        
        M = [M, m(2:end)];
        S = [S, s(2:end)];
        V = [V, repmat(v0, 1, numBins-1)];
        tM = [tM, tGrid(2:end)+currT];
        
        %splitFactor = 1/2;
        if (Z(l)==1 && divideVolume<5)%if mitosis, split cell in two halfs until the 5th cell cycle
            factor = splitFactor;
            divideVolume = divideVolume + 1;
        else
            factor = 1;
        end
        m0 = m(end)*factor;
        s0 = s(end)*factor^2;
        
        v0 = v0*factor;
    end
    
    %ComputeConditionalRNASolution(t, Z);
    
    X = [M;S;V];
    XS = SampleCTMPPathGrid_mex(X, tM, grid);
    
    MS(k, :) = XS(1, :);
    SS(k, :) = XS(2, :);
    VS(k, :) = XS(3, :);
    
    ZS(k, :) = SampleCTMPPathGrid_mex(Z, t, grid);
    
    %stairs(t/60, Z); hold on;
    %plot(tM/60, M); hold off; pause;
end