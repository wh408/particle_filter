%a_bar = ones(1,N);

%syms y;
        
%pdf = exp(-0.5*y^2)/sqrt(6.28);

%Gaussian = int(pdf,-inf,inf);

%omega = zeros(N,steps+1);

%omega(j,i) = omega(j,i-1) + (1/sqrt(N))*sigma*sigma*dt;

%pi_t = pi_t + a(j,i)*(v(j,i)*Gaussian);