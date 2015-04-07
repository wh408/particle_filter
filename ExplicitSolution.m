function [xthat] = ExplicitSolution(T, steps, F, f, sigma, H, h, Yt)

Rt = zeros(1,steps+1);

xthat = ones(1,steps+1);

dt = T / steps;

for i = 2:(steps+1)
   
    Rt(i) = Rt(i-1) + (sigma^2+2*H*Rt(i-1)-(H*Rt(i-1))^2)*dt;
    
    xthat(i) = xthat(i-1) + (F*xthat(i-1)+f)*dt + H*Rt(i-1)*((Yt(i)-Yt(i-1)) - (H*xthat(i-1)+h)*dt);
    
end