function y = fermenter_f2(t,C_g,D,Sf)

X = C_g(1); %initial guess for biomass concentration
S = C_g(2); %initial guess for substrate concentration
P = C_g(3); %initial guess for product concentration

mu_m = 0.48/60; %min^-1 
P_m = 50; %g/l
K_m = 1.2; %g/l
K_i = 22; %g/l


mu = mu_m*(1-P/P_m)*S/(K_m+S+S^2/K_i);

alpha = 2.2; %g/g
beta = 0.2/60; %min^-1
Y_xs = 0.4; %g/g

y(1) = -D*X+ mu*X;
y(2) = D*(Sf-S)-(mu*X)/(Y_xs);
y(3) = -D*P+ (alpha*mu+beta)*X;

y = y';

end