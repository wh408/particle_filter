function [T steps F f sigma H h Yt] = GenerationOfObservation(T, steps, F, f, sigma, H, h)

% the signal: dXt = (F*Xt+f)*dt + sigma*dVt
% the observation: dYt = (H*Xt+h)*dt + dWt

dt = T / steps;

Xt = ones(1,steps+1);

Yt = ones(1,steps+1); 

for i=2:(steps+1)
    
    Xt(i) = Xt(i-1) + (F*Xt(i-1)+f)*dt + sigma*randn(1,1)*sqrt(dt); 
    
    Yt(i) = Yt(i-1) + (H*Xt(i-1)+h)*dt + randn(1,1)*sqrt(dt);
    
end