function [Y_A,Y_H,f_p,i_XB,i_XE,i_XP,mu_H,K_s,K_OH,K_NO,b_A,b_H,eta_g,eta_h,k_h,K_x,mu_A,K_NH,K_OA,k_a] = ASM_param(T)

% parameters for AS system at 15 degC, based on BSM1
mu_H = 4.0;  %6.0;
K_s = 10.0;  %20;
K_OH = 0.2;
K_NO = 0.5;
b_H = 0.3;  %0.62;
mu_A = 0.5;  %0.8;
K_NH = 1.0;
K_OA = 0.4;
b_A = 0.05;  %0.2;
eta_g = 0.8;
k_a = 0.05;  %0.08;
k_h = 3.0;
K_x = 0.1;  %0.03;
eta_h = 0.8;  %0.4;
Y_H = 0.67;
Y_A = 0.24;
f_p = 0.08;
i_XB = 0.08;  %0.086;
i_XP = 0.06;
i_XE = i_XP;

mu_H = mu_H*exp((log(mu_H/3.0)/5.0)*(T-15.0)); %/* Compensation from the current temperature in the reactor */
b_H = b_H*exp((log(b_H/0.2)/5.0)*(T-15.0));
mu_A = mu_A*exp((log(mu_A/0.3)/5.0)*(T-15.0));
b_A = b_A*exp((log(b_A/0.03)/5.0)*(T-15.0));
k_h = k_h*exp((log(k_h/2.5)/5.0)*(T-15.0));
k_a = k_a*exp((log(k_a/0.04)/5.0)*(T-15.0));
end