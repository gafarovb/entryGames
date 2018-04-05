nMC = 10000;


beta1 = [.4;.5];
theta(1:2,1) =  beta1 ;
delta1 = 1;
theta(3,1) = delta1 ;

beta2 = [0.4;.5];
theta(4:5,1) =  beta2 ;
delta2 = 1;
theta(6,1) = delta2 ;
seed = 222353;
[Y1,X] = DGP(nMC,seed,[beta1;beta2;delta1;delta2],01,0);
clear w;
for i = 1 : nMC
    x = X(i,:);
    switch find(Y1(i,:))
        case 1
            outcome = [0;0];
        case 2
            outcome = [1;1];
        case 3
            outcome = [1;0];
        otherwise
            outcome = [0;1];
    end
    w(i,1) = DOP(X(i,:),outcome) ;
end
order =1;
nOrdersq = (order+1)^2;

hMat =[ zeros((order+1)^2,2) ones(nOrdersq,1) ];

gamma_0 =[theta;reshape(hMat,3*nOrdersq,1)];





options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');
options.MaxFunctionEvaluations = 6000;

x0 = gamma_0;
fun = @(gamma) -loglikelihood(w ,gamma)/nMC;
A = [zeros(nOrdersq,6) -eye(nOrdersq)  -eye(nOrdersq)  -eye(nOrdersq);
    zeros(nOrdersq,6) eye(nOrdersq)  eye(nOrdersq)  eye(nOrdersq) ];

b = [zeros(nOrdersq,1);...
    ones(nOrdersq,1)];
Aeq = [];
beq = [];
eps = 0.2;
lb = [eps*ones(6,1);zeros(3*nOrdersq,1)];
ub = [3*ones(6,1);ones(3*nOrdersq,1)];




[gamma_MLE, fval,exitflat,~,lambda,~]  = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options);
%%
options = optimoptions('fmincon','Algorithm','active-set','Display','iter');
eps = 0.2;
lb = [-eps*ones(6,1);-eps*ones(3*nOrdersq,1)];
[gamma_MLE_1, fval,exitflat,~,lambda,~]  = fmincon(fun,gamma_MLE,[],[],Aeq,beq,lb,ub,[],options);

