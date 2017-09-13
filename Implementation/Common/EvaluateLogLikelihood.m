function [L] = EvaluateLogLikelihood(TargetMoments, TargetUncertainties, PredictedMoments)

    MomentsT(1, :) = TargetMoments.Mean;
    MomentsT(2, :) = TargetMoments.SCV;

    MomentsP(1, :) = PredictedMoments.Mean;
    MomentsP(2, :) = PredictedMoments.SCV;
    
    Uncertainties(1, :) = TargetUncertainties.Mean;
    Uncertainties(2, :) = TargetUncertainties.SCV;

    exclIdx = isnan(Uncertainties(1, :));
    
    L = -1/2 * nansum(nansum((MomentsT - MomentsP).^2 ./ Uncertainties));

end