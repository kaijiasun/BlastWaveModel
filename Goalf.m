function [Err,SPT]=Goalf(param)
global datax datay dataerr
global R0 s2 c1 c2 as faip  m0 tao0 gscal

global Iter

pt = datax;
best_T=param(1);best_rho=param(2);xi=param(3);
SPT=Cooper_noas(m0,pt,best_T,best_rho,tao0);
SPT=SPT.'*xi*gscal;

Err =sqrt(mean((SPT-datay).^2./dataerr.^2));
v=tanh(best_rho);
Teff = best_T*sqrt((1+v)/(1-v));

Iter = Iter+1;
[Iter Teff param Err]
