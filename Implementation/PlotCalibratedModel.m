clear;
close all;


addpath('Common/');
addpath('../Data/');



ending = 'DEL'; 


ccResults = load(['CellCycleResults' ending '.mat']);
txData = load(['DataProcessed' ending '.mat']);


load(['ResultsTranscriptionModel' ending '.mat']);

counter = 1;
M = 500;


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
        

        randFrac = ccResults.bestParams(1);
        tm = ccResults.bestParams(2);
        cycleLength = ccResults.bestParams(3:end-1);
        
        % Extract mean and variance from FISH data using bootstrapping
        [means, vars, scvs] = BootstrapUncertainties(data, 3000);
        
        targetMoments.Mean = mean(means, 2);
        targetMoments.SCV = mean(scvs, 2);
        targetMoments.Variance = mean(vars, 2);
        targetUncertainties.Mean = var(means, [], 2);
        targetUncertainties.SCV = var(scvs, [], 2);
        targetUncertainties.Variance = var(vars, [], 2);
        
        %% Select best inference run
        Runs = Results{geneIdx, chIdx}.Runs;
        
        for j=1:numRuns

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
       
        
        [MS, SS, VS, ZS] = SimulateRNAMoments(M, grid, tm, randFrac, cycleLength, a, b, c1, c2, m0, s0);
        
        
        
        predictedMoments.Mean = mean(MS./VS);
        predictedMoments.Variance = mean(SS./VS.^2) - predictedMoments.Mean.^2;
        
        
        figure(1);
        subplot(2,4,counter);
        title(GeneNames{n});
        plot(grid/60, predictedMoments.Mean, '-b', grid/60, predictedMoments.Mean - sqrt(predictedMoments.Variance), '--b', grid/60, predictedMoments.Mean + sqrt(predictedMoments.Variance), '--b'); hold on;
        plot(measurementTimes/60, targetMoments.Mean, '-r.', measurementTimes/60, targetMoments.Mean - sqrt(targetMoments.Variance), '--r.', measurementTimes/60, targetMoments.Mean + sqrt(targetMoments.Variance), '--r.'); hold off;
        xlabel('Time');
        ylabel('Transcript density mean and std');
        
        figure(2);
        subplot(2,4,counter);
        title(GeneNames{n});
        plot(grid/60, predictedMoments.Mean, '-b'); hold on;
        plot(measurementTimes/60, targetMoments.Mean, 'ro'); hold off;
        xlabel('Time');
        ylabel('Mean transcript density');
        
        figure(3);
        subplot(2,4,counter);
        title(GeneNames{n});
        plot(grid/60, predictedMoments.Variance./predictedMoments.Mean.^2, '-b'); hold on;
        plot(measurementTimes/60, targetMoments.Variance ./ targetMoments.Mean.^2, 'ro'); hold off;
        xlabel('Time');
        ylabel('SCV of transcript density');
      
        drawnow;

  
        fprintf('%s &%f &%f &%f &%f &%f &%f\n', GeneNames{counter}, a/b/(60*60), sqrt(a)/b/(60*60), c1*60*60, c2*60*60, bestParams(5), bestParams(6)/bestParams(5)^2);
        
        counter = counter + 1;
    end
end
    
    

