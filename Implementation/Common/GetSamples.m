function [Time, Samples, Stats] = GetSamples(Stages, fields, channelIdx, embryoIdx)

if (nargin<4)
    embryoIdx = [];
end

NumStages = length(Stages);

Stats = struct;

for k=1:NumStages
    
    Time(k) = Stages{k}.Time;
    
    if isempty(embryoIdx)
        validIdx = 1:length(Stages{k}.Embryo);
    else
        validIdx = embryoIdx;
    end
    
    
    Samples{k} = [];
    
    for u=1:length(fields)
        samples = [];
        
        
        for l=validIdx
            emb = Stages{k}.Embryo{l};
            val = getfield(emb.Channel{channelIdx}, fields{u});
            samples = [samples; val];
        end
        
        Samples{k} = [Samples{k}, samples];
        Stats = setfield(Stats, ['Mean_' fields{u}], {k}, mean(samples));
        Stats = setfield(Stats, ['Var_' fields{u}], {k}, var(samples));
    end
    
    
end
