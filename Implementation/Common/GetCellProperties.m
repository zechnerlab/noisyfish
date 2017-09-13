function [Time, Samples] = GetCellProperties(Stages, fields, embryoIdx)

if (nargin<3)
    embryoIdx = [];
end

NumStages = length(Stages);


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
            val = getfield(emb, fields{u});
            samples = [samples; val];
        end
        
        Samples{k} = [Samples{k}, samples];
    end
    
    
end
