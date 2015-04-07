function [intpart] = intp(number)

if number >= 0 
    
    intpart = int32(fix(number));
    
else
        
    intpart = int32(fix(number)) - 1;
    
end 