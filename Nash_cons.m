function [ C_eq,C1,Cnen,Cnne,c7 ] = Nash_cons(d1)
points=halton(2,d1);
t1=points(:,1);
t2=points(:,2);

G=[ones(d1,1) t1 t2  t1.*t2 t1.^2 t2.^2];
C1=[G zeros(d1,6) zeros(d1,6); zeros(d1,6) G zeros(d1,6);
    zeros(d1,6) zeros(d1,6) G];
%C2=[G G G];
%C01=[C1;C2;-C1];

Sigma1=t1.*t2;
Sigma2=(1-t1).*(1-t2);
Sigma3=(1-t2).*t1;
%%%%%%%% Eq  %%%%%
C_eq=[repmat(Sigma2,1,6).*G -repmat(Sigma1,1,6).*G zeros(d1,6)];
%c_eq=zeros(d1,1);
%%%%%%%% Ineq not EN %%%%
Cnen=[zeros(d1,6) repmat(Sigma3,1,6).*G -repmat(Sigma2,1,6).*G];
%%%%%%%% Ineq not NE %%%%
Cnne=[repmat(1-Sigma2-Sigma3,1,6).*G repmat(Sigma1,1,6).*G repmat(Sigma1,1,6).*G];
c7=Sigma1;
%%%%%%%%%%%%%%%%%%%%%%%%%%
 
end

