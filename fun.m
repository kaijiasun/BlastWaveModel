function y = fun(xin,param)

xi=param(1);
T=param(2);
m0=param(3);
q=param(4);
x= sqrt(m0^2+xin.^2)/T;
y=xi*(1+(q-1)*x).^(-1/(q-1));