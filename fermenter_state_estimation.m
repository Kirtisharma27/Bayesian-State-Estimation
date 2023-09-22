clear all
clc

global u1 u2 Pk Yk Xk_pred_hat

sample_T = 1; %min

n_op = 3; %no. of outputs(Biomass, Substrate and Product concentration)
n_ip = 2; %no. of inputs(Dilution Rate, Substrate feed concentration)

n_st = 3; %no. of states

N_samples = 1500;

%Initialisation of model parameters
mu_m = 0.48/60; %min^-1 
P_m = 50; %g/l
K_m = 1.2; %g/l
K_i = 22; %g/l

%Other parameters
alpha = 2.2; %g/g
beta = 0.2/60; %min^-1
Y_xs = 0.4; %g/g

%Operating conditions
state_sigma = [0.002 0.002 0.002]; %standard deviation of noise in X,S,P

meas_sigma = [0.02 0.02]; %standard deviation of noise in measurements(X,S)

%State and Measurement Noise Covariance Matrices

R = diag([meas_sigma(1)^2;meas_sigma(2)^2]); %measurement noise matrix

Q = diag([state_sigma(1)^2,state_sigma(2)^2,state_sigma(3)^2]); %state noise covariance matrix

%steady state values of state and manipulated inputs

x_ss = [7.038, 2.404, 24.87]'; %initial values of X,S,P

u_ss = [0.15/60, 20]'; %Manipulated input values(D,Sf)

Xk = x_ss; %absolute states
Uk(:,1) = u_ss; %absolute inputs


%Noise generation in states

for i=1:1:1500
    
    wk(1,i) = state_sigma(1)*randn; %noise generation in X
    wk(2,i) = state_sigma(2)*randn; %noise generation in S
    wk(3,i) = state_sigma(3)*randn; %noise generation in P
    
end


%Noise generation in measurements

for i=1:1:1500
    
    vk(1,i) = meas_sigma(1)*randn; %noise generation in X
    vk(2,i) = meas_sigma(2)*randn; %noise generation in S
    
end

%PBRS Input Signal (Pseudo Random Binary Input Signal)

%u1_D = 0.15/60 + idinput(1500,'PRBS',[0 0.166],[-0.0005 0.0005]); %Dilution rate
%u2_Sf = 20 + idinput(1500,'PRBS',[0 0.166],[-0.75 0.75]); %Substrate feed concentration


%Input signal generation

D_ss = 0.15/60; %min^-1
Sf_ss = 20; %g/l


for i=1:1:1500
    
    if i<50
        
        u1_D(i) = D_ss;
        u2_Sf(i) = Sf_ss;
        
    elseif i>=50 && i<200
        
        u1_D(i) = D_ss+0.1*D_ss;
        u2_Sf(i) = Sf_ss;
        
    elseif i>=200 && i<=1500
        
        u1_D(i) = D_ss+0.1*D_ss;
        u2_Sf(i) = Sf_ss+0.1*Sf_ss;
        
    end
    
end







%Measurement matrix

C = [1 0 0; 0 1 0]; %only X and S values are measured

%Initial inputs in absolute form

Uk(1:2,1) = u_ss;
Yk(:,1) = C*Xk(1:3,1) + vk(:,1);


%Initial state of estimator

Xkhat = x_ss + 0.01*x_ss; % 1% deviation from steady state is assumed


%Initial covariance matrix
%Since initial state is unknown, initial covariance is assumed to be large

Pk = 5*Q;

%Initial inputs

u1 = 0.15/60; %Dilution rate(min^-1)
u2 = 20; %Substrate feed concentration(g/l)


%Plant and Estimator simulation

for i = 1:1:N_samples
    
    %Plant simulation
    
    %Input profiles
    
    u1 = u1_D(1,i);
    u2 = u2_Sf(1,i);
    
    %u1 = u1_D(i,1);
    %u2 = u2_Sf(i,1);
    
    
    [t,Xt] = ode45(@fermenter_f2,[0, sample_T], Xk',[],u1,u2);
    
    Xk = Xt(length(t),:)'+ wk(:,i); %generating X(k+1) 
    
    Yk(:,i) = C*Xk + vk(:,i);
    
    
    %Extended kalman filter(state estimator)
    
    %Prediction step
    
    [T, Xj] = ode45(@fermenter_f2,[0, sample_T], Xkhat',[],u1,u2);
    
    Xk_pred_hat = Xj(end,:)' ;
    
    %generate linear pertubation model
    
    z_vec = (Xkhat)';
    
    [A_mat] = fermenter_jacob(z_vec);
    
    phy = expm(A_mat*sample_T);
    
    Pk = phy*Pk*phy' + Q;
    
    var = R + C*Pk*C';
    
    Lk = Pk*C'*inv(var); %kalman gain computation
    
    %updated step
    Xkhat = Xk_pred_hat + Lk*(Yk(:,i)-C*Xk_pred_hat);
    
    Pk = (eye(3,3) - Lk*C)*Pk;
    
    
    %simulation results
    
    res_matrix(i,:) = [i Xk(1) Xk(2) Xk(3) Xkhat(1) Xkhat(2) Xkhat(3) Yk(1,i) Yk(2,i) u1 u2];
    
    
end


figure(1)
plot(res_matrix(:,1),res_matrix(:,4),'b', res_matrix(:,1),res_matrix(:,7),'r')
xlabel('Time(min)')
ylabel('Product Concentration(g/l)')
legend('True state', 'Estimated state')

figure(2)
plot(res_matrix(:,1),res_matrix(:,8),'k',res_matrix(:,1),res_matrix(:,5),'r',res_matrix(:,1),res_matrix(:,2),'b' )
xlabel('Time(min)')
ylabel('Biomass Concentration(g/l)')
legend('Measurement','Estimated state','True state')

figure(3)
plot(res_matrix(:,1),res_matrix(:,9),'k',res_matrix(:,1),res_matrix(:,6),'r',res_matrix(:,1),res_matrix(:,3),'b')
xlabel('Time(min)')
ylabel('Substrate Concentration(g/l)')
legend('Measurement','Estimated state','True state')

    
    
%Metrics to compute Sum of Squared Errors (SSE)

esterr1 = res_matrix(50:end,2) - res_matrix(50:end,5);
esterr2 = res_matrix(50:end,3) - res_matrix(50:end,6);
esterr3 = res_matrix(50:end,4) - res_matrix(50:end,7);

SSE1 = esterr1'*esterr1;
SSE2 = esterr2'*esterr2;
SSE3 = esterr3'*esterr3;

figure(4)
plot(res_matrix(50:end,1),esterr1)
xlabel('Time(min)')
ylabel('Biomass Concentration(g/l)')
legend('Error Plot')

figure(5)
plot(res_matrix(50:end,1),esterr2)
xlabel('Time(min)')
ylabel('Substrate Concentration(g/l)')
legend('Error Plot')


figure(6)
plot(res_matrix(50:end,1),esterr3)
xlabel('Time(min)')
ylabel('Product Conncentration(g/l)')
legend('Error Plot')


SSEa = [SSE1 SSE2 SSE3];

    
    
    




















