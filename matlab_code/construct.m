function set = construct(A,B,v)
%% 吸引不变集构建
eps=0.00001;
[V,J]=jordan(A);
V=real(V);
J=real(J);
[m,n]=size(J);
b=inv(eye(m)-abs(J))*abs(inv(V)*B)*v+eps;
A_=inv(V);
A=[A_;-A_];
b=[b;b];
set=Polyhedron(A,b);
end

