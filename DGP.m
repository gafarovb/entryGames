%%%%%%%%%%%% DGP %%%%%%%%%%%
function [Y1,X]=DGP(n,seed,b,a,be)
rng(seed,'twister');
E=randn(n,4);
e1=E(:,1);
e2=E(:,2);
X=E(:,3:4);
x1=X(:,1);
x2=X(:,2);

u1u=b(1)+b(2)*x1;
u2u=b(3)+b(4)*x2;
u1l=u1u-b(5);
u2l=u2u-b(6);


%%%%%%%%%%%%%% 2LRational %%%%%%%%%%%%%%%%%%%%% 
Iee=and(e1<=u1l,e2<=u2l);
Inn=and(e1>=u1u,e2>=u2u);
Ine=or(and(e1>u1l,e2<=u2l),and(e1>=u1u,e2<u2u));
Ien=or(and(e1<=u1l,e2>u2l),and(e1<u1u,e2>=u2u));
Im=not(Iee|Inn|Ine|Ien);

%%%%%%%%%%%%%% NE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  (NN,EE,EN,NE)
a1=(e1-u1l)./(u1u-u1l);
a2=(e2-u2l)./(u2u-u2l);
sigma1=repmat([0,0,1,0],n,1);
sigma2=repmat([0,0,0,1],n,1);
sigma3=[a1.*a2 (1-a2).*(1-a1) a1.*(1-a2) (1-a1).*a2];
%%%%%%%% Fixed selection %%%%%%%%%%%%%%%
%a=.3; be=.3;
sigma=a*sigma1+be*sigma2+(1-a-be)*sigma3;
%%%%%%%% Continious selection %%%%%%%%%%%%%
%sigma=repmat((x1/2),1,4).*sigma1+repmat((x2/2),1,4).*sigma2+repmat((1-(x1+x2)/2),1,4).*sigma3;

%%%%%%%%%%%%%% Distribution %%%%%%%%%%%%%%
P=zeros(n,4);
P(Inn,:)=repmat([1,0,0,0],sum(Inn),1);
P(Ien,:)=repmat([0,0,1,0],sum(Ien),1);
P(Ine,:)=repmat([0,0,0,1],sum(Ine),1);
P(Iee,:)=repmat([0,1,0,0],sum(Iee),1);
P(Im,:)=sigma(Im,:);

%%%%%%%%%%%%%% Otcomes %%%%%%%%%%%%%%%%%
Y1=P;
rng(seed+2,'twister');
z=rand([sum(Im),1]);
Y1(Im,:)=[(z<=P(Im,1)) ...
    and(z>P(Im,1),z<=(P(Im,1)+P(Im,2))) ...
    and(z>(P(Im,1)+P(Im,2)),z<=(P(Im,1)+P(Im,2)+P(Im,3))) ...
    and(z>(P(Im,1)+P(Im,2)+P(Im,3)),z<=1)];
%Y=[Y1(:,1) Y1(:,4) Y1(:,2) Y1(:,3)];
%save('data.mat','YX','n', 'dx', 'beta_0');
%sum(Im)
%sum(YX,1)

end