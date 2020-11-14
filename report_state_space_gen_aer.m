function [ A1,B1,C1,D1 ] = state_space_gen_aer(Q,V,W,phi,dim )
%State space generation
% kL_a = 10.26; % [1/hr] 240 [1/d] Mass transfer coefficient of the aeration equipment under the operating conditions imposed. 10 degrees values
% J = [zeros(7,1)' kL_a zeros(5,1)']';
J = [W*phi(:,8)]; % oxygen is taken to be the final input
phi(:,8) = [];
A = W*phi - ((Q/V)*eye(dim)); B = [(Q/V)*eye(dim) J];
C = eye(dim); D = zeros(dim,dim+size(J,2));

tank_ss = ss(A,B,C,D);
A1 = tank_ss.A;
B1 = tank_ss.B;
C1 = tank_ss.C;
D1 = tank_ss.D;

end