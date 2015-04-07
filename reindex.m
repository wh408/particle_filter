function [v2] = reindex(v, o)

    d = find(o >0);

    v2 = zeros(1,10);

    for jj = 1:length(d)
        
        for jk = 1:o(d(jj))
            
            if jj == 1
                
                v2(jk) = v(d(jj));
                
            else
       
                v2(sum(o(d(1):d(jj-1)))+jk) = v(d(jj));
                
            end
        
        end
    
    end   