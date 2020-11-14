function [A1,B1,C1,D1,rho] = mean_squares_aer(x0,Q,V,T)
% Generates the state-space linearized model for given state x0

[Y_A,Y_H,f_p,i_XB,i_XE,i_XP,mu_H,K_s,K_OH,K_NO,b_A,b_H,eta_g,eta_h,k_h,K_x,mu_A,K_NH,K_OA,k_a] = ASM_param(T);

% define switch functions
Sub_switch1 = x0(2)/(K_s+x0(2));
NO_switch1 = x0(9)/(K_NO+x0(9));
O2_switch1 = x0(8)/(K_OH+x0(8));
O2_switch2 = K_OH/(K_OH+x0(8));
NH4_switch1 = x0(10)/(K_NH+x0(10));

% Process rates
rho1 = mu_H*Sub_switch1*O2_switch1*x0(5); % aerobic growth of heterotrophs rate 
rho2 = mu_H*Sub_switch1*O2_switch2*NO_switch1*eta_g*x0(5); % anoxic growth of heterotrophs rate
rho3 = mu_A*NH4_switch1*O2_switch1*x0(6); % aerobic growth of autotrophs rate
rho4 = b_H*x0(5); % decay of heterotrophs rate
rho5 = b_A*x0(6); % decay of autotrophs rate
rho6 = k_a*x0(11)*x0(5); % ammonification of soluble organic nitrogen rate
rho7 = k_h*((x0(4)/x0(5))/(K_x+(x0(4)/x0(5))))*(O2_switch1+ eta_h*O2_switch2*NO_switch1)*x0(5); % hydrolysis of entrapped organics rate
rho8 = rho7*x0(12)/x0(4); % hydrolosis of entrapped organic nitrogen rate

rho = [rho1 rho2 rho3 rho4 rho5 rho6 rho7 rho8]';

M1 = [x0(2), x0(5), x0(8)];
M2 = [x0(2), x0(5), x0(8), x0(9)];
M3 = [x0(6), x0(8), x0(10)];
M6 = [x0(5), x0(11)];
M7 = [x0(4), x0(5), x0(8), x0(9)];
M8 = [x0(5), x0(8), x0(9), x0(12)];

theta1 = M1.\rho1; % finds coefficients [a12 a15 a18]
theta2 = M2.\rho2; % coefficients [a22 a25 a28 a29]
theta3 = M3.\rho3; % [a36 a38 a310]
theta4 = b_H; % a45
theta5 = b_A; % a56
theta6 = M6.\rho6; % [a65 a611]
theta7 = M7.\rho7; % [a74 a75 a78 a79]
theta8 = M8.\rho8; % [a85 a88 a89 a812]

phi = [0 theta1(1) 0 0 theta1(2) 0 0 theta1(3) 0 0 0 0 0; ...
    0 theta2(1) 0 0 theta2(2) 0 0 theta2(3) theta2(4) 0 0 0 0; ...
    0 0 0 0 0 theta3(1) 0 theta3(2) 0 theta3(3) 0 0 0; ...
    0 0 0 0 theta4 0 0 0 0 0 0 0 0; ...
    0 0 0 0 0 theta5 0 0 0 0 0 0 0; ...
    0 0 0 0 theta6(1) 0 0 0 0 0 theta6(2) 0 0; ...
    0 0 0 theta7(1) theta7(2) 0 0 theta7(3) theta7(4) 0 0 0 0; ...
    0 0 0 0 theta8(1) 0 0 theta8(2) theta8(3) 0 0 theta8(4) 0];

W = [zeros(1,8); ...
    -1/Y_H -1/Y_H zeros(1,4) 1 0; ...
    zeros(1,8); ...
    zeros(1,3) 1-f_p 1-f_p 0 -1 0; ...
    1 1 0 -1 zeros(1,4); ...
    0 0 1 0 -1 zeros(1,3); ...
    zeros(1,3) f_p f_p zeros(1,3); ...
    -(1-Y_H)/Y_H 0 (1- (4.57/Y_A)) zeros(1,5); ...
    0 -(1-Y_H)/2.86*Y_H 1/Y_A zeros(1,5); ...
    -i_XB -i_XB -i_XB-(1/Y_A) 0 0 1 0 0; ...
    zeros(1,5) -1 0 1; ...
    zeros(1,3) -i_XB-f_p*i_XB -i_XB-f_p*i_XB 0 0 -1; ...
    -i_XB/14 ((1-Y_H)/142.86*Y_H)-i_XB/14 -(i_XB/14)-1/7*Y_A 0 0 1/14 0 0];
W(8,:) = [];
dim = length(W(:,1));

[ A1,B1,C1,D1 ] = state_space_gen_aer(Q,V,W,phi,dim );

end