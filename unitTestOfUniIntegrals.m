
beta1 = [1;2];
theta(1:2,1) =  beta1 ;
delta1 = 2;
theta(3,1) = delta1 ;

beta2 = [1;1];
theta(4:5,1) =  beta2 ;
delta2 = 1;
theta(6,1) = delta2 ;

x = [  1 0 ];
xB_1 = [1 x(1)]* beta1;
n = size(x,1);
sieves = DOP(x);
order = sieves.order;
analyticIntegrals = sieves.uniIntegrals(theta);


mu_v_1 = 1 - [ones(n,1) x(:,1)]*beta1 /delta1;
sigma_v_1 = 1/delta1^2;


nMC = 1000000;
v_0 = randn(nMC,1);
v = mu_v_1 + sigma_v_1 * v_0;

% test of integrals of powers of standard normal
for i = 0:order
    MCIntegrals_1(i+1,1) = mean(v_0.^i.*(v_0>0).*(v_0<1));
end
error_1 = MCIntegrals_1 - (sieves.stdIntegral( 1) - sieves.stdIntegral(0));


% test of integrals of powers of non-standard normal
for i = 0:order
    MCIntegrals_2(i+1,1) = mean(v.^i.*(v>0).*(v<1));
end
error_2 = MCIntegrals_2 - sieves.SPfromI_ix(delta1,xB_1)*(sieves.stdIntegral( (1-mu_v_1)/sigma_v_1) - sieves.stdIntegral(-mu_v_1/sigma_v_1));


% test of integrals of bernstein polynomials of non-standard normal

for i = 0:order
    MCIntegrals(i+1,1) = mean(sieves.bernstein_i(i,v).*(v>0).*(v<1));
end

error_3 = analyticIntegrals(:,1)-MCIntegrals;
