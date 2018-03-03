 
tic 
theta = [1 2 0.1 2 3 0.3]';
parfor ix=1:100000
    a.evalIntegral(theta);
end
toc