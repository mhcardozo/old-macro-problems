clearvars;
close all;
tic

global v0 beta delta alpha kmat k0 z0 Zprob zmat s kgrid j i L sigma

%% set parameters
alpha = 0.33; % capital's share
beta = 0.95;
delta= 0.1; % depreciation rate (annual)
s = 2;

gamma = (1/(1-beta*alpha)); %parameter for the analytical value function


tol = 0.01;
maxits = 300;
dif = tol+1000;
its = 0;
kgrid = 99; % grid points + 1

%% Solve for steady state
 kstar = (alpha/(1/beta - (1-delta)))^(1/(1-alpha)); % steady state k
 cstar = kstar^(alpha)-delta*kstar;
 istar = delta*kstar;
 ystar = kstar^(alpha);

 kappa = 0.6;
 kmin = (1-kappa)*kstar; % minimum capital
 kmax = (1+kappa)*kstar; % maximum capital
 grid = (kmax-kmin)/kgrid; % grid

 kmat = kmin:grid:kmax;

 kmat = kmat';

 [N,n] = size(kmat);

 %% Discretize stochastic process
mu = 0;
rho = 0.95;
L = 7;
m = 3;
sigma = 0.01;
 
 [Z, Zprob] = tauchen(L,mu,rho,sigma,m);
 zmat = exp(Z);
 
 
%% Iterate Value Function
 
 v0 = zeros(N,L);

 while dif>tol && its < maxits
  for j = 1:L
        for i = 1:N
        k0 = kmat(i,1);
        z0 = zmat(j,1);
        k1 = fminbnd(@valfun2s,kmin,kmax);
        v1(i,j) = -valfun2s(k1);
        k11(i,j) = k1;
        end
 end
 dif = norm(v1-v0)
 v0 = v1;
 its = its+1
 end

%% Compute policy functions from ValFunct It. and Log Linear

 cons = zeros(N,1); % policy function for consumption
 consl = zeros(N,1);
 vfun = zeros(N,1);
 for j =1:L
 for i = 1:N
 cons(i,1)= zmat(j,1)*kmat(i,1)^(alpha)-k11(i,j)+(1-delta)*kmat(i,1); % consumption
 consl(i,1) = cstar + cstar*0.6793*(zmat(j,1)-1)+ cstar*0.4270*(kmat(i,1)-kstar)/kstar; %log linear consumption
 kl(i,1) = kstar + kstar*0.2163*(zmat(j,1)-1)+ 0.8979*(kmat(i,1)-kstar); % log linear capital
 vfun(i,1)= gamma*log(zmat(j,1)*kmat(i,1)^(alpha))+ (log(1-alpha*beta)+beta*gamma*log(beta))/(1-beta); %analytical value function
 end
 end
  
 
 
 toc
 
%% Plot
 figure
 plot(kmat,k11)
 title('capital policy function')
 xlabel('k_t')
 ylabel('k_(t+1)')
 saveas(gcf,'cap','png')
 
 figure
 plot(kmat,v1)
 title('value function')
 xlabel('k_t')
 ylabel('V(k_t)')
 saveas(gcf,'valfun','png')
 
 

 figure
 plot(kmat,cons)
 title('consumption policy function')
 xlabel('k_t')
 ylabel('c_t')
 
 saveas(gcf,'cons','png')
 
figure
hold on
plot(kmat,cons)
plot(kmat,consl)
title('Comparison of consumption')
hold off
saveas(gcf,'consL','png')
 
figure
hold on
plot(kmat,k11)
plot(kmat,kl)
title('Comparison of capital')
hold off
saveas(gcf,'capL','png')