function[A_x] = fermenter_jacob(X_k)

global u1 u2 

A_x =[];

X = X_k(1);
S = X_k(2);
P = X_k(3);

Sf = u2;
D = u1;

%Initialisation of model parameters
mu_m = 0.48/60; %min^-1 
P_m = 50; %g/l
K_m = 1.2; %g/l
K_i = 22; %g/l

%Other parameters
alpha = 2.2; %g/g
beta = 0.2/60; %min^-1
Y_xs = 0.4; %g/g

mu = mu_m*(1-P/P_m)*S/(K_m+S+S^2/K_i);

dmu_ds = mu_m*(1-P/P_m)*(K_m-S^2/K_i)/(K_m + S + S^2/K_i)^2;

dmu_dp = -mu_m*S/(P_m*(K_m + S + S^2/K_i));

A_x(1,1) = -D + mu;
A_x(1,2) = X*(dmu_ds);
A_x(1,3) = X*(dmu_dp);

A_x(2,1) = -mu/Y_xs;
A_x(2,2) = -D - X*(dmu_ds)/Y_xs;
A_x(2,3) = -X*(dmu_dp)/Y_xs;

A_x(3,1) = alpha*mu + beta;
A_x(3,2) = alpha*X*(dmu_ds);
A_x(3,3) = -D + alpha*X*(dmu_dp);



