
function [o] = BranchingAlgorithm (a_bar, N)

MC = 1000;

o = zeros(MC,N);

for k = 1:MC
    
    u = rand(1,N-1);

    g = N; h = N;

    for j = 1:(N-1)
    
        if (fracp(N*a_bar(j))) + (fracp(g-N*a_bar(j))) < 1
        
            if u(j) < 1 - ((fracp(N*a_bar(j))) / (fracp(g)))
            
                o(k,j) = intp(N*a_bar(j));
            
            else
            
                o(k,j) = intp(N*a_bar(j)) + (h - intp(g));
            
            end 
        
        else
        
            if u(j) < 1 - (1-(fracp(N*a_bar(j)))) / (1-(fracp(g)))
            
                o(k,j) = intp(N*a_bar(j)) + 1;
            
            else
            
                o(k,j) = intp(N*a_bar(j)) + (h - intp(g));
            
            end
            
        end 
    
        g = g - N*a_bar(j);
    
        h = h - o(k,j);
    
    end

    o(k,N) = h;
    
end

