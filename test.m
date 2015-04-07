
% the signal: dXt = (F*Xt+f)*dt + sigma*dVt

T=1;
steps=10000;
F=2;
f=0.5;
sigma=0.25;
dt = T / steps;

H=0.5;
h=0.2;

Xt=ones(steps,1);
Yt=ones(steps,1);



for i=2:(steps+1)
    
    Xt(i) = Xt(i-1) + (F*Xt(i-1)+f)*dt + sigma*randn(1,1)*sqrt(dt); 
    
    Yt(i) = Yt(i-1) + (H*Xt(i-1)+h)*dt + randn(1,1)*sqrt(dt);
    
end

%%explicit solution, p152 in Dan's book  
[xthat] = ExplicitSolution(T, steps, F, f, sigma, H, h, Yt);

%%simple Monte-Carlo
%N=100;
%[pi_t] = ParticleFilterWithoutBranching(T, steps, N, F, f, sigma, H, h, xthat, Yt);

%%particle filter
N=500;
M=5;
[pi_t a] = ParticleFilter(T, steps, M, N, F, f, sigma, H, h, xthat, Yt);

figure

plot(Yt,'r')
hold on
plot(Xt);
hold on
plot(xthat,'c')
hold on
plot(pi_t,'g')

