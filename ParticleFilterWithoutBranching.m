function [pi_t] = ParticleFilterWithoutBranching(T, steps, N, F, f, sigma, H, h, xthat, Yt)

dt = T / steps;

v = ones(N,steps+1); % sample the initial position from pi_0

b = zeros(N,steps+1);
    
a = ones(N,steps+1); % initial weights 

a_bar = zeros(N,steps+1);
    

for j = 1:N
    
    for i = 2:(steps+1)
        
        v(j,i) = v(j,i-1) + (F*v(j,i-1)+f)*dt + sigma*randn(1,1)*sqrt(dt);
        
        b(j,i) = (H*v(j,i-1)+h)*(Yt(i) - Yt(i-1)) - (dt/2)*((H*v(j,i-1)+h)^2);
        
        a(j,i) = a(j,i-1)*exp(b(j,i));
        
    end
    
end


pi_t = zeros(1,steps+1);
 
for i = 1:(steps+1)

    for j=1:N
        
        a_bar(j,i) = a(j,i) / sum(a(:,i));
        
        pi_t(i) = pi_t(i) + a_bar(j,i)*v(j,i);
        
    end
    
end


figure(1)
plot(1:(steps+1),xthat,'r', 1:(steps+1),pi_t,'g', 1:(steps+1),Yt,'b');
xlabel('time step'); ylabel('value');
legend('E[Xt|Yt]', 'Parcitle Approximation', 'Observation');
title('Three processes');

figure(2)
plot(1:(steps+1),(xthat-pi_t)./xthat,'b-d');
title('Ratio distance between the signal and its approximation');
