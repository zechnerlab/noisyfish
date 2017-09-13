function [t, Act] = SimulateActivationTime(T, a, b)
    tau = gamrnd(a, 1/b);
    
    t = [0, tau, T];
    Act = [0, 1, 1];
    
end