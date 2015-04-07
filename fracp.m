function [fracpart] = fracp (number)

d = double(intp(number));

if number >= 0.0
    
    fracpart = number - d;
    
else
    
    fracpart = number - (d+1.0);
    
end