clearvars;
close all;

beta=0.95;
alpha=0.33;
s=2;
delta=0.1;
rho=0.95;

kstar=((1-(1-delta)*beta)/(alpha*beta))^(1/(alpha-1));
cstar=kstar^alpha-delta*kstar;
ystar=kstar^alpha;

A(1,1)=1-(1-delta)*beta;
A(1,2)=(1-(1-delta)*beta)*(alpha-1);
A(1,3)=-s;
A(2,2)=kstar;
A(3,1)=1;
B(1,3)=-s;
B(2,1)=ystar;
B(2,2)=alpha*ystar+(1-delta)*kstar;
B(2,3)=-cstar;
B(3,1)=rho;

C=inv(B)*A;

[V,lam] = eig(C);
P=inv(V);
if (lam(1,1)<1 &lam(2,2)<1) | (lam(1,1)<1 &lam(3,3)<1) | (lam(2,2)<1 & lam(3,3)<1) | (lam(1,1)>1 & lam(2,2)>1 & lam(3,3)>1)
display('No saddle path stability.');
else
display('Blanchard-Khan conditions hold.')
dlam= diag(lam);
t1=flipud(sortrows([dlam V']));
lam=diag(t1(:,1));
V=t1(:,2:4);
P=inv(V');

lam1=lam(1:2,1:2);
lam2=lam(3,3);
P11=P(1:2,1:2);
P12=P(1:2,3);
P21=P(3,1:2);
P22=P(3,3);

display('Coefficients for chat.')
z1 = -inv(P22)*P21
display('Coefficients for [zhat; khat]')
z2 = inv(P11-P12*inv(P22)*P21)*inv(lam1)*(P11-P12*inv(P22)*P21)
end

