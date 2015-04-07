function [ml1,epsilon,dist1] =  pf(x,y,te,sig,nu,n,N,alpha,beta,x0,ep)

    u = zeros(N,1);
    
    wei = zeros(N,1);
    
    dist = zeros(N,1);
    
    vec = zeros(N,1);
    
    tee = 1.0/te;
    
    ml1 = 0.0;
    
    vec1 = zeros(100,1);
    
    for i=1:100
       
        pp = beta*x0 + tee*nu*randn(1,1);
        
        vec1(i) = abs(y(1)-alpha*pp + sig*randn(1,1));
    end
    
    ep1 = max(vec1);
   
    for i=1:n
        
       if i==1 % step 0
       
           epsilon(i) = ep1;
           
           %epsilon(i) = ep;
           
           l1 = 0;
           
           for j=1:N
              
               x(j) = beta*x0 + tee*nu*randn(1,1);
               
               u(j) = alpha*x(j) + sig*randn(1,1);
               
               dist(j) = abs(u(j)-y(i));
               
               dist1(i,j) = dist(j);
               
               if  dist(j) < epsilon(i)
                   
                    wei(j) = 0.0;
                    
                    l1 = l1 + 1;
               else
                    
                    wei(j) = -300.0;
               end
           end

           if l1>0 % what does this if-else do?
               
               ml1 = ml1 + log(l1/N);
           else
              
               ml1 = ml1 - 300.0;
           end
       else % step 1: i >= 2

           vec = sort(dist); 
           
           epsilon(i) = vec(floor(0.95*N));

           sa = max(wei);
            
            nc = 0.0;
            
            for j=1:N
               
                wei(j) = exp(wei(j)-sa);
                
                nc = nc + wei(j);
            end
            
            for j=1:N
                
                wei(j) = wei(j)/nc; % to normalise the weights
            end
            
            for j= 1:N
               
                l = inver(wei,N);
                
                vec(j) = x(l);
                
                if i > 2
                
                    for l1=1:(i-2)
                   
                        dist1(l1,j) = dist1(l1,l);
                    end
                end
                
                dist1(i-1,j) = dist(l);
            end
            
            l1 = 0;
            
            for j =1:N % step 3
                
                x(j) = beta*vec(j) + tee*nu*randn(1,1);
               
                u(j) = alpha*x(j) + sig*randn(1,1);
               
                dist(j) = abs(u(j)-y(i));
                
                dist1(i,j) = dist(j);
               
                if  dist(j) < epsilon(i)
                   
                     wei(j) = 0.0;
                     
                     l1 = l1 + 1;
                else
                    
                     wei(j) = -300.0;
                end
            end
            
            if l1 >0 
                
                ml1 = ml1 + log(l1/N);
            else
                
               ml1 = ml1 - 300.0; 
            end
       end
    
    end
%%end of function