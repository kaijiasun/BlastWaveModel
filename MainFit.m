function [param,sp,Err]=MainFit(data,pt)
global datax datay dataerr
global R0 s2 c1 c2 as faip  m0 tao0 gscal
global Iter
R0=13.6;
tao0 =11;
s2=0.0;   % 0,no anisotropy in coordinate space
as = 0; % no dispersion
faip = 0; %
c1= 0 ;  % 0,no anisotropy in momentum space
c2 = 1.8;  %1, suppression
gscal = 3;
m0 = 0.5;%2.2865;%1.86484;

datax=data{1};
datay=data{2};
dataerr=data{3};
x0=[0.15,1,1.3]; %  T, rho, xi ;
ximin = 1E-2;ximax=1E3;
Tmin=0.06; Tmax=0.25;
rhomin=0.2;rhomax=1.1;
%rhomin=1.03;rhomax=1.03;
Iter = 0;
options = psoptimset('TolFun',1E-6,'TolCon',1E-6,'MaxIter',1000,'MaxFunEvals',1e4);
%x0 = [2 0.175 0 0 1 1];
[x,fval,exitflag,output]=patternsearch(@Goalf,x0,[],[],[],[],[ Tmin rhomin ximin],[Tmax rhomax ximax],[],options);
x
fval
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(@Goalf,x,[],[],[],[],[ Tmin rhomin ximin],[Tmax rhomax ximax],[],options);
param=x;
param

[Err,y]=Goalf(param);
figure;semilogy(datax,datay,'-r.')
hold on;semilogy(datax,y,'-b.')

best_T=param(1);best_rho=param(2);xi=param(3);
SPT=Cooper_noas(m0,pt,best_T,best_rho,tao0);
sp=SPT*xi*gscal;
