nMC = 1000;


beta1 = [.4;.5];
theta(1:2,1) =  beta1 ;
delta1 = 1;
theta(3,1) = delta1 ;

beta2 = [0.4;.5];
theta(4:5,1) =  beta2 ;
delta2 = 1;
theta(6,1) = delta2 ;
seed = 222353;
[Y1,X] = DGP(nMC,seed,[beta1;beta2;delta1;delta2],0.333,0.333);
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
order =2;

 
hMat =[  ones((order+1)^2,3)/4];

 gamma_0 =[theta;reshape(hMat,3*(order+1)^2,1)];

 
 
 
 
options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');

 
 x0 = gamma_0;
fun = @(gamma) -loglikelihood(w ,gamma)/nMC;
 A = [zeros(9,6) -eye(9)  -eye(9)  -eye(9);
             zeros(9,6) eye(9)  eye(9)  eye(9) ];
 
 b = [zeros(9,1);...
             ones(9,1)];
 Aeq = [];
 beq = [];
 eps = 0.3;
 lb = [eps*ones(6,1);zeros(3*9,1)];
 ub = [3*ones(6,1);ones(3*9,1)];

 


[gamma_MLE, fval,exitflat,~,lambda,~]  = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options);
 
 
