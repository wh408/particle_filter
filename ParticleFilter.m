
function [pi_t a] = ParticleFilter(T, steps, M, N, F, f, sigma, H, h, xthat, Yt)

dt = T / steps;

delta = M * dt;

Nbranching = int32(fix(T/delta));

v = ones(N,steps+1); % sample the initial position from pi_0

b = zeros(N,steps+1);
    
a = ones(N,steps+1); % initial weights

pi_t = zeros(1,steps+1);

a_bar = ones(N,steps+1);

a_bar(:,1) = 1/N;

pi_t(1) = 1;


for k = 1:Nbranching 
    
% Following is the code for the evolution of the particle system after 
% (k-1)th branching and before kth branching.

    for j = 1:N 
    
        for i = ((k-1)*M+2):(k*M+1)
        
            v(j,i) = v(j,i-1) + (F*v(j,i-1)+f)*dt + sigma*randn(1,1)*sqrt(dt);
        
            b(j,i) = (H*v(j,i-1)+h)*(Yt(i) - Yt(i-1)) - (dt/2)*((H*v(j,i-1)+h)^2);
        
            a(j,i) = a(j,i-1)*exp(b(j,i));
        
        end
    
    end
    
    for i = ((k-1)*M+2):(k*M+1)
        
        for j = 1:N
           
            a_bar(j,i) = a(j,i) / sum(a(:,i));
            
            pi_t(i) = pi_t(i) + a_bar(j,i)*v(j,i);
            
        end
        
    end
    


% The branching algorighm described in section 9.2.1, in order to obtain o:
    
    u = rand(1,N-1);
    
    o = zeros(1,N);

    g = N; h1 = N;

    for j = 1:(N-1)
    
        if (fracp(N*a_bar(j,k*M+1))) + (fracp(g-N*a_bar(j,k*M+1))) < 1
        
            if u(j) < 1 - ((fracp(N*a_bar(j,k*M+1))) / (fracp(g)))
            
                o(j) = intp(N*a_bar(j,k*M+1));
            
            else
            
                o(j) = intp(N*a_bar(j,k*M+1)) + (h1 - intp(g));
            
            end 
        
        else
        
            if u(j) < 1 - (1-(fracp(N*a_bar(j,k*M+1)))) / (1-(fracp(g)))
            
                o(j) = intp(N*a_bar(j,k*M+1)) + 1;
            
            else
            
                o(j) = intp(N*a_bar(j,k*M+1)) + (h1 - intp(g));
            
            end
            
        end 
    
        g = g - N*a_bar(j,k*M+1);
    
        h1 = h1 - o(j);
    
    end

    o(N) = h1;
    



% Re-index procedure: to set values to v(j,k*M+1), 1<=j<=N, according to
% the vector o:

    d = find(o >0);

    v2 = zeros(N,1);

    for jj = 1:length(d)
        
        for jk = 1:o(d(jj))
            
            if jj == 1
                
                v2(jk) = v(d(jj),k*M+1);
                
            else
       
                v2(sum(o(d(1):d(jj-1)))+jk) = v(d(jj),k*M+1);
                
            end
        
        end
    
    end   
    
% re-assign values for the positions and un-normalised weights:

    for j = 1:N
        
        v(j,k*M+1) = v2(j);
        
    end
        
    a(:,k*M+1) = 1;
    
end


% the evolution after the final branching
%TBA

A = 1:(steps+1);

figure(1)
plot(A,xthat,'r', A,pi_t,'g', A,Yt,'b');
xlabel('time step'); ylabel('value');
legend('E[Xt|Yt]', 'Parcitle Approximation', 'Observation');
title('Three processes');

figure(2)
plot(A,(xthat-pi_t)./xthat,'b-d');
xlabel('time step'); ylabel('value');
title('Ratio error between the signal and its approximation');
