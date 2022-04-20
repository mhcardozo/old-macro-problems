
 function val=valfun2(k)

 % This program gets the value function for a neoclassical growth model with
 % no uncertainty and CRRA utility

 global v0 beta delta alpha kmat k0 sigma kgrid


 klo=max(sum(k>kmat),1); % identify the gridpoint that falls just below . .
 % . . the choice for k
 khi=klo+1;

 % do the interpolation
 gg = v0(klo) + (k-kmat(klo))*(v0(khi) -v0(klo))/(kmat(khi) - kmat(klo));

 c = k0^alpha - k + (1-delta)*k0; % consumption
 if c<0
 val = -9999999 - 999*abs(c);
 else
 val = (1/(1-sigma))*(c^(1-sigma)-1) + beta*gg;
 end
 val = -val; % make it negative since we're maximizing and code is to minimize.
