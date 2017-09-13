clear;
close all;


addpath('Common/');
addpath('../Data/');



ending = 'DEL'; 


ccResults = load(['CellCycleResults' ending '.mat']);
txData = load(['DataProcessed' ending '.mat']);


load(['ResultsTranscriptionModel' ending '.mat']);

counter = 1;

% Evaluate / plot genes in the same order as in Stapel et al., Genes Dev
% (2017)
permMat{1} = [1, 1];
permMat{2} = [3, 1];
permMat{3} = [2, 1];
permMat{4} = [4, 1];
permMat{5} = [1, 2];
permMat{6} = [3, 2];
permMat{7} = [2, 2];
permMat{8} = [4, 2];

for u=1:length(txData.Genes)
    
    
    for n=1:2
        
        permV = permMat{counter};
        geneIdx = permV(1);
        chIdx = permV(2);
        
        breakNames = split(txData.Genes{geneIdx}.Name, '_');
        GeneNames{counter} = breakNames{chIdx}; 
        [~, data, stat] = GetSamples(txData.Genes{geneIdx}.Stages, {'TranscriptDens'}, chIdx);
        [~, volData] = GetCellProperties(txData.Genes{geneIdx}.Stages, {'AreaCell'});
        exclIdx = [1,2,3];
        allIdx = 1:length(data);
        remIdx = setdiff(allIdx, exclIdx);
        data = data(remIdx);
        volData = volData(remIdx);
        
        t0 = ccResults.bestParams(end);
        measurementTimes = ccResults.Time(1:end) - min(ccResults.Time) + t0;
        grid = linspace(0, measurementTimes(end), 200);
        longGrid = measurementTimes;
        
        
        %% cell-cycle parameters
        randFrac = ccResults.bestParams(1);
        tm = ccResults.bestParams(2);
        cycleLength = ccResults.bestParams(3:end-1);

        % Extract mean and variance from FISH data
        [means, vars, cvs] = BootstrapUncertainties(data, 1000);

        targetMoments.Mean = mean(means, 2);
        targetMoments.Variance = mean(cvs, 2);
        targetUncertainties.Mean = var(means, [], 2);
        targetUncertainties.Variance = var(cvs, [], 2);
        
        
        %% Select best inference run
        Runs = Results{geneIdx, chIdx}.Runs;
        
        for j=1:length(Runs)
            LOptVec(j) = Runs{j}.LOpt;
        end
        
        [maxL, maxIdx] = max(LOptVec);
        bestRun = Runs{maxIdx};
        bestParams = bestRun.bestParams;
        
        a = bestParams(1);
        b = bestParams(2);
        c1 = bestParams(3);
        c2 = bestParams(4);
        m0 = bestParams(5);
        var0 = bestParams(6);
        s0 = m0^2 + var0;
       
        %Set number of samples used to evaluate inner and outer conditional
        %expectations (see Supplement of Stapel et al, Genes Dev (2017))
        M = 500;
        L = 200;
        
        
        %Simulate M conditional moments for the same gene activation time for L
        %different activation times
        for k=1:L
            [MS, ~, VS, ~] = SimulateRNAMoments(M, longGrid, tm, randFrac, cycleLength, a, b, c1, c2, m0, s0, 1);
            condVars(k, :) = var(MS./VS);
            condMeans(k, :) = mean(MS./VS);
            fprintf('Iteration %d\n', k);
        end
        
        [MS, SS, VS, ZS] = SimulateRNAMoments(2000, longGrid, tm, randFrac, cycleLength, a, b, c1, c2, m0, s0, 0);
        
        
        % compute cell-cycle, intrinsic and activation noise
        cellcycleVar = mean(condVars);
        varsCondBoth = SS./VS.^2 - (MS./VS).^2;
        RNAVar = mean(varsCondBoth);
        activationVar = var(condMeans);
        meanTot = mean(condMeans);
        
        ccSCV = cellcycleVar ./ meanTot.^2;
        RNASCV = RNAVar ./ meanTot.^2;
        actSCV = activationVar ./ meanTot.^2;
        
        figure(1);
        subplot(2,4,counter);
        bar(longGrid/60, [RNASCV; actSCV; ccSCV]', 'Stacked'); hold on;
        xlabel('Time');
        ylabel('SCV of transcript density');
        
        SCVVec(counter) = totSCV(end);
        
        drawnow;

        counter = counter + 1;
    end
end

save NoiseDecomposition.mat;
    
    

