function [ f ] = logBarrier( x )
%LOGBARRIER Summary of this function goes here
%   Detailed explanation goes here
alpha = 0.011;
cutoff = 1 - alpha;
if size(x,1)>1

    for ix = 1:size(x,1)
        f(ix) = logBarrier( x(ix) );
    end
else
    

 
 
L = 2*(1 - cutoff);
k = 4/L;
dY = L/2-cutoff ;

 if x>cutoff
   f =  L/(1+exp(-k*(x-cutoff)))-dY;
 elseif x<alpha
   f =  L/(1+exp(-k*(x-alpha)));  
 else
     f = x;
 end
end
end

