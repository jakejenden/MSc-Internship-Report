%% Prep
close all; clc; clear all;

%% Data processing 

% Fields of dataset: C1:t [6 hrs] C2:Si [g COD m^-3] C3:Ss [g COD m^-3] C4:Xi [g COD m^-3] 
% C5:Xs [g COD m^-3] C6:Xbh [g COD m^-3] C7:Xba [g COD m^-3] C8:Xp [g COD
% m^-3] C9:So [g -COD m^-3]
% C10:Sno [g N m^-3] C11:Snh [g N m^-3] C12:Snd [g N m^-3] C13:Xnd [g N
% m^-3] C14:Salk [mole m^-3] C15:Q [m^3 d^-1]
% Notes: S_O units -COD is equivalent to DO consumption
% Benchmark WWTP data extraction- Inf_rain_2006.txt
content = fileread('Inf_rain_2006.txt');
data_inf = textscan(content,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');

data_inf = [data_inf{2} data_inf{3} data_inf{4} data_inf{5} data_inf{6} data_inf{7}, ...
    data_inf{8} data_inf{9} data_inf{10} data_inf{11} data_inf{12} data_inf{13} data_inf{14} data_inf{15}];

% extend data for more datailed steps
data_inf_ext = zeros(6*length(data_inf)-5,14);
data_inf_ext(1,:) = data_inf(1,:);
up =1;
data_pts = 6;
for j = 1:(length(data_inf)-1)
    for i = 1:data_pts
        data_inf_ext(i+up,:) = (i/data_pts).*(data_inf(j+1,:)-data_inf(j,:)) + data_inf(j,:);
    end
    up = up + 6;
end

%% Using BSM1 simulation set up
% Initial state space for aerobic tank
V = 1*1333; % [m^3]
Q = 53377.6074;% Influent flow to AS %open-loop = 92230, closed-loop = 53377.6074
% Q/V = D (dilution rate)
T =15; % [degrees] temperature of system
cut = 4000; % number of points to be erased

% Create time vector
content = fileread('Inf_rain_2006.txt');
data_inf1 = textscan(content,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
data_inf1 = data_inf1{1};
t0 = zeros(6*length(data_inf1)-5,1);
t0(1) = data_inf1(1);
up =1;
data_pts = 6;
for j = 1:(length(data_inf1)-1)
    for i = 1:data_pts
        t0(i+up,:) = (i/data_pts).*(data_inf1(j+1,1)-data_inf1(j,1)) + data_inf1(j,1);
    end
    up = up + data_pts;
end
t0 = t0(1:(end-cut));
label = ["S_I","S_s","X_I","X_s","X_{BH}","X_{BA}","X_p","S_{NO}","S_{NH}","S_{ND}","X_{ND}","S_{ALK}"];
state_dim = length(label);

% Aerobic tank 5

% values inside reactor 5 from the BSM1 benchmark from WWTP simulation
% (steady state values)
SI_prox = 30; SS_prox = 0.80801; XI_prox = 1149.1683; % [mg COD/l]
XS_prox = 44.4828; XBH_prox = 2562.8514; XBA_prox = 154.163; % [mg COD/l]
XP_prox = 452.7367; % [mg COD/l] 
SO_prox = 2; % [mg -COD/l]
SNO_prox = 13.5243; SNH_prox = 0.67193; SND_prox = 0.6645; XND_prox = 3.2605; % [mg N/l]
SALK_prox = 3.8277; % [mol HCO3/m3]

x0_aer5 = [SI_prox SS_prox XI_prox XS_prox XBH_prox XBA_prox XP_prox SO_prox SNO_prox SNH_prox SND_prox XND_prox SALK_prox]'; % tank initial condition (ensure no zero elements)
x0 = [30; 6.9; 46.95; 108.3; 133.2; 4.9; 0.036; 19.7; 11.5; 4; 4.8; 6.83];

[A1,B1,C1,D1,rho_aer5] = mean_squares_aer(x0_aer5,Q,V,T); % mean squares linearisation method
aer_tank5 = ss(A1,B1,C1,D1);

isstable(aer_tank5) % quick stability check [1=stable,0=not stable]

figure(180); clf;
pzmap(aer_tank5) % pole-zero map

set_pt = 2; % [mg/l] O2 set-point
data_len = length(data_inf_ext); % data length

inp = [data_inf_ext(:,1:7) data_inf_ext(:,9:13) set_pt*ones(data_len,1)];

adj = length(inp) - cut+1; % cut off excess data
inp(adj:end,:) = [];
y5 = lsim(aer_tank5,inp,t0',x0);

figure(85); %clf;
for i = 1:state_dim 
    subplot(4,3,i)
    plot(t0,y5(:,i))
    title(label(i))
    hold on
end
hold on