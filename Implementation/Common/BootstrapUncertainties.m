function [means, vars, scvs] = BootstrapUncertainties(data, M)

L = length(data);

for l=1:L
    K = length(data{l});
    N = K;
    rndIdx = randi(K, M, N);
    
    dat = data{l};
    
    for i=1:M
        rndIdx = randi(K, 1, K);
        means(l, i) = mean(dat(rndIdx));
        vars(l, i) = var(dat(rndIdx));
        scvs(l, i) = vars(l, i)./means(l, i).^2;
    end
end


end