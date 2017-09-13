%Script for calibrating the cell cycle model
%Author: Christoph Zechner

clear;
close all;


addpath('Common/');
addpath('Data');


% Specify if analysis should be performed on DEL or EVL cells
maskIdx = 1;

if maskIdx==3
    ending = 'DEL'; 
elseif maskIdx==1
    ending = 'EVL';
end

load(['MitosisData' ending '.mat']);


T = max(Time);
grid = linspace(0, T, 500);

%some starting values
ta = 6*60;
tm = 8*60;
t0 = 8*60;
randFrac = 0.4;
cycleLength = [17, 30, 40, 120, 140]*60; 

TimeTmp = Time-min(Time)+t0;

% number MC samples used for computing the cell-cycle model
M = 1000;

% number of MCMC iterations used for estimating the parameters
L = 3000;

chain = zeros(3+length(cycleLength), L);
chain(:, 1) = [randFrac, tm, cycleLength, t0];

%set sigma of the lognormal proposal distribution
sigmaProp = 0.08;

PM_data = NumMitosis ./ NumCells;


LOpt = -inf;


%% MCMC scheme to fit the cell-cycle data
for i=1:L

    currParams = chain(:, i);
    
    forwardProb = 1;
    backwardProb = 1;
    
    for l=1:length(currParams)
       propParams = lognrnd(log(currParams), sigmaProp);
       forwardProb = forwardProb * lognpdf(propParams(l), log(currParams(l)), sigmaProp);
       backwardProb = backwardProb * lognpdf(currParams(l), log(propParams(l)), sigmaProp);
    end
    
    % evalulate likelihood
    for k=1:M
        

        randFrac = propParams(1);
        tm = propParams(2);
        cycleLength = propParams(3:end-1);
        t0 = propParams(end);
        
        TimeTmp = Time-min(Time)+t0;
        fullGrid = linspace(0, max(TimeTmp), 100);
        
        [Z, t] = SimulateCellCycle2(T, tm, randFrac, cycleLength);
        ZS(k, :) = SampleCTMPPathGrid_mex(Z, t, TimeTmp);
        ZSFull(k, :) = SampleCTMPPathGrid_mex(Z, t, fullGrid);
    end
    
    PM = mean(ZS);
    PMFull = mean(ZSFull);
    
    for k=1:length(PM)
        
        kFact = sum(log(1:NumMitosis(k)));
        nFact = sum(log(1:NumCells(k)));
        nmkFact = sum(log(1:(NumCells(k)-NumMitosis(k))));
        
        L(k) = nFact - kFact - nmkFact + NumMitosis(k).*log(PM(k))+(NumCells(k) - NumMitosis(k)).*log(1-PM(k)); %binomial likelihood
    end
    
    LNew = sum(L);
    if (i==1)
        accept = 1;
    else
        accept = min(1, exp(LNew - LOld)*backwardProb/forwardProb);
    end
    
    if (rand<accept)
       fprintf('Accepted sample\n');
       
       chain(:, i+1) = propParams;
       LOld = LNew;
       
       if (LNew>LOpt)
           bestPM = PMFull;
           bestTimeFull = fullGrid;
           bestTime = TimeTmp;
           bestParams = propParams;
           LOpt = LNew;
       end
    else
       fprintf('Rejected sample\n');
       
       chain(:, i+1) = chain(:, i);
    end
    

    LVec(i) = LOld;
    
    if (mod(i, 10)==0)
        subplot(1,3,1);
        plot(bestTimeFull/60, bestPM, 'o-', bestTime/60, PM_data, 'xr'); 
        title('Best Fit');
        xlabel('Time');
        ylabel('P(mitosis)');
        
        subplot(1,3,2);
        plot(LVec(1:i));
        xlabel('MCMC iterations');
        ylabel('Log-Likelihood');
        
        subplot(1,3,3);
        plot(chain(2:end, 1:i)');
        xlabel('MCMC iterations');
        ylabel('Parameters');
        
        drawnow;
    end

end

save(['CellCycleResults' ending '.mat']);




