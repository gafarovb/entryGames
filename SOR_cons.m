function [ C,c ] = SOR_cons(d1)
points=halton(2,d1);
t1=points(:,1);
t2=points(:,2);

f=[ones(d1,1) t1 t2  t1.*t2 t1.^2 t2.^2];
C1=[f zeros(d1,size(f,2)) zeros(d1,size(f,2)); zeros(d1,size(f,2)) f zeros(d1,size(f,2));
    zeros(d1,size(f,2)) zeros(d1,size(f,2)) f];
C2=[f f f];
C=[C2;-C1];
c=[ones(d1,1); zeros(3*d1,1)];
end

