%Script for calibrating the cell cycle model
%Author: Christoph Zechner

clear;
close all;



addpath('Common/');
addpath('../Data/');

% Parameter inference only performed for DEL cells because data in EVL
% cells seems to be very noisy
ending = 'DEL'; 

ccResults = load(['CellCycleResults' ending '.mat']);
txData = load(['DataProcessed' ending '.mat']);

%half lives from PCR measurements in Stapel et al., Genes Dev (2017)
HLs = [1065 516
        185 56
        74 109
        1016 204]*60;

M = 300;
L = 600;
numRuns = 5;

% Run %numRuns independent MCMC inference for each of the genes
for u=1:length(txData.Genes)
    
    for n=1:2
        
        geneIdx = u;
        chIdx = n;
        [~, data, stat] = GetSamples(txData.Genes{geneIdx}.Stages, {'TranscriptDens'}, chIdx);
        [~, volData] = GetCellProperties(txData.Genes{geneIdx}.Stages, {'AreaCell'});
        
        % remove first three time points before transcription activation
        exclIdx = [1,2,3];
        allIdx = 1:length(data);
        remIdx = setdiff(allIdx, exclIdx);
        data = data(remIdx);
        volData = volData(remIdx);
        
        v0 = 75000;
        
        t0 = ccResults.bestParams(end);
        measurementTimes = ccResults.Time(1:end) - min(ccResults.Time) + t0;
        grid = measurementTimes;
        
        
        % set cell cycle parameters from calibrated cell cycle model
        randFrac = ccResults.bestParams(1);
        %randFrac = 0.2;
        tm = ccResults.bestParams(2);
        cycleLength = ccResults.bestParams(3:end-1);
        
        % Extract mean and variance from FISH data using bootstrapping
        [means, vars, scvs] = BootstrapUncertainties(data, 3000);
        
        targetMoments.Mean = mean(means, 2);
        targetMoments.SCV = mean(scvs, 2);
        targetUncertainties.Mean = var(means, [], 2);
        targetUncertainties.SCV = var(scvs, [], 2);
        
        
        c2Vec = logspace(-5, -2, numRuns);
        for j=1:numRuns
            
            % Set starting values
            meanC2 = log(2)/HLs(u, n);
            aGam = 2;
            bGam = aGam / meanC2;

            c2 = aGam / bGam;
            meanRNA = 3;
            c1 = meanRNA*c2;
            
            meanTau = 35*60;
            a = 4;
            b = a*1/meanTau;
            m0 = targetMoments.Mean(1);
            s0 = targetMoments.SCV(1)*m0^2 + m0^2;

            numParams = 6;
            chain = zeros(numParams, L);
            chain(:, 1) = [a, b, c1, c2, m0, s0-m0^2];
            
            % randomize initial condition
            chain(:, 1) = lognrnd(log(chain(:, 1)), 0.3);
            
            sigmaProp = 0.05;
            LOpt = -inf;
            
            for i=1:L
                currParams = chain(:, i);
                
                % computer forward and backward probabilities of the
                % proposal
                forwardProb = 1;
                backwardProb = 1;
                
                for l=1:length(currParams)
                    propParams = lognrnd(log(currParams), sigmaProp);
                    forwardProb = forwardProb * lognpdf(propParams(l), log(currParams(l)), sigmaProp);
                    backwardProb = backwardProb * lognpdf(currParams(l), log(propParams(l)), sigmaProp);
                end
                
                a = propParams(1);
                b = propParams(2);
                c1 = propParams(3);
                c2 = propParams(4);
                
                m0 = propParams(5);
                var0 = propParams(6);
                s0 = m0^2 + var0;
                
                % Forward-simulate model
                [MS, SS, VS, ZS] = SimulateRNAMoments(M, grid, tm, randFrac, cycleLength, a, b, c1, c2, m0, s0);
                
                
                predictedMoments.Mean = mean(MS./VS);
                predictedMoments.SCV = (mean(SS./VS.^2) - predictedMoments.Mean.^2) ./ predictedMoments.Mean.^2;
                
                % Evaluate parameter log-posterior
                LNew = EvaluateLogLikelihood(targetMoments, targetUncertainties, predictedMoments) + log(gampdf(c2, aGam, 1/bGam));
                
                
                if (i==1)
                    accept = 1;
                else
                    accept = min(1, exp(LNew - LOld) * backwardProb / forwardProb);
                end
                
                if (rand<accept)
                    chain(:, i+1) = propParams;
                    LOld = LNew;
                    
                    if (LNew>LOpt)
                        LOpt = LNew;
                        bestParams = propParams;
                        bestMoments = predictedMoments;
                    end
                else
                    chain(:, i+1) = chain(:, i);
                end
                
                LVec(i) = LOld;
                
                if (mod(i, 10)==0)
                    subplot(2,2,1);
                    m0 = bestParams(5);
                    var0 = bestParams(6);
                    cv0 = var0 / m0^2;
                    
                    plot([0, measurementTimes]/60, [m0, bestMoments.Mean], '-b'); hold on;
                    plot(measurementTimes/60, targetMoments.Mean, 'or'); hold off;
                    xlabel('Time');
                    ylabel('Mean transcript density');
                    title('Best fit');
                    
                    subplot(2,2,2);
                    plot([0, measurementTimes]/60, [cv0, bestMoments.SCV], '-b'); hold on;
                    plot(measurementTimes/60, targetMoments.SCV, 'or'); hold off;
                    xlabel('Time');
                    ylabel('SCV of transcript density');
                    title('Best fit');
                    
                    subplot(2,2,3);
                    plot(LVec(1:i));
                    xlabel('Iterations');
                    ylabel('Log-Posterior');
                    
                    subplot(2,2,4);
                    plot(log10(chain(:, 1:i)'));
                    xlabel('Iterations');
                    ylabel('Log-Parameters');
                    drawnow;
                end
  
            end
            
            Runs{j}.bestParams = bestParams;
            Runs{j}.LOpt = LOpt;
            Runs{j}.bestMoments = bestMoments;
            
            fprintf('Finished Run %d (%d %d)\n', j, u, i);
        end
      
        Results{u, n}.Runs = Runs;
        
    end
    
end


save(['ResultsTranscriptionModel' ending '.mat']);

