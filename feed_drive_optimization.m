clear; clc; close all;

%% Import UBC
Tr_set = load('UBC Data/Data_Trial_4.mat');
Te_set = load('UBC Data/Data_Trial_13.mat');
%% Set a range around the wanted regions

% we set the default values for start and end indices. If you want to select a specific region of 
% data for identification, uncomment the code below, and comment the
% default values
Ranges = load('UBC Data/Ranges.mat');
Ns4 = Ranges.Ns(4);
Ne4 = Ranges.Ne(4);
Ns13 = Ranges.Ns(13);
Ne13 = Ranges.Ne(13);

%Set time for each data set
omega = 60/3000;
t = Tr_set.t(Ns4:Ne4);
dt = t(end)/(length(t)-1);  % sample rate of data

%% Convert data
%Downsample filter
dS = 30; %Downsampling Rate
N = 10000; %Sampling window
fc = 1/(2*dt*dS); %Cutoff Freq
fs = 1/dt; %Sampling Freq
[b,a] = butter(6,fc/(fs/2));

t_tr = downsample(Tr_set.t(Ns4:Ne4),dS);
t_te = downsample(Te_set.t(Ns13:Ne13),dS);

v_r0 = Tr_set.VelocityXNom(Ns4:Ne4);
v_m0 = Tr_set.VelocityXAct(Ns4:Ne4);
vr_tr = downsample(filtfilt(b,a,v_r0*0.001/60),dS); %% mm/min --> m/sec
vm_tr = downsample(filtfilt(b,a,v_m0*0.001/60),dS); %% mm/min --> m/sec

v_r0 = Te_set.VelocityXNom(Ns13:Ne13);
v_m0 = Te_set.VelocityXAct(Ns13:Ne13);
vr_te = downsample(filtfilt(b,a,v_r0*0.001/60),dS); %% mm/min --> m/sec
vm_te = downsample(filtfilt(b,a,v_m0*0.001/60),dS); %% mm/min --> m/sec

Fx_tr = downsample(Tr_set.Fx(Ns4:Ne4),dS);
Fx_te = downsample(Te_set.Fx(Ns13:Ne13),dS);

y_tr = vm_tr(2:end);  % output
phi_tr = [vm_tr(1:end-1), vr_tr(1:end-1)]; % regressor array
t_tr = t_tr(1:end-1); % time

y_te = vm_te(2:end);  % output
phi_te = [vm_te(1:end-1), vr_te(1:end-1)]; % regressor array
t_te = t_te(1:end-1); %time

%% Create data vecotors for PL-LSSVM method
% %As per Espinoza,
% omega = 60/3000;
% N = 1000;
% t = 0:0.003:(N-1)*0.003;
% noise = 0.001*rand(N,1);
% u = rand(N,1);
% y(1) = 0.3*u(1) + noise(1);
% y(2) = 0.3*u(2) - 0.2*u(1) + noise(2);
% for i = 3:N
%     y(i) = -0.5*y(i-2) + 0.3*u(i) - 0.2*u(i-1) + noise(i);
% end
% y2 = (y + (0.3*sin((2*pi()/omega)*t)))';
% 
% %As per Espinoza, discard first 350 data
% phi_tr = [y2(350:749),u(352:751),u(351:750)];
% y_tr = y2(352:751);
% t_tr = t(352:751)';
% 
% %As per Espinoza,
% N = 1000;
% t = 0:0.003:(N-1)*0.003;
% noise = 0.001*rand(N,1);
% u = rand(N,1);
% y(1) = 0.3*u(1) + noise(1);
% y(2) = 0.3*u(2) - 0.2*u(1) + noise(2);
% for i = 3:N
%     y(i) = -0.5*y(i-2) + 0.3*u(i) - 0.2*u(i-1) + noise(i);
% end
% y2 = (y + 0.2*sin((2*pi/2.4)*t).*sin((2*pi()/omega)*t))';
% 
% %As per Espinoza, discard first 350 data
% phi_te = [y2(350:749),u(352:751),u(351:750)];
% y_te = y2(352:751);
% t_te = t(352:751)';

%% Import Munich Data
% %Load Data From Files
% T_0 = 'Munich Data\Feedstep150_N0_H55.sco';
% T_con = 'Munich Data\Feedstep150_N300_H55.sco';
% T_sweep = 'Munich Data\Feedstep150_Nsweep300_H55.sco';
% 
% tmp1 = load([T_con(1:end-3) 'mat']);
% t_con = tmp1.t_in;
% tmp1 = tmp1.AI;
% con_shaker = tmp1(:,1)*1000;
% [Fq_con_shak, Gy_con_shak] = fftspec(con_shaker', 500, 4e-4);
% 
% tmp1 = load([T_sweep(1:end-3) 'mat']);
% t_sweep = tmp1.t_in;
% tmp1 = tmp1.AI;
% sweep_shaker = tmp1(:,1)*1000;
% [Fq_sweep_shak, Gy_sweep_shak] = fftspec(sweep_shaker', 500, 4e-4);
% 
% [T_0_data] = read_tncscope_file(T_0);
% [T_con_data] = read_tncscope_file(T_con);
% [T_sweep_data] = read_tncscope_file(T_sweep);
% 
% %Set Inital Range Around The Trajectory
% Ns = 1880;
% Ne = 4600;
% 
% %Collect And Store Data From Inital Trajectory
% v_r = T_0_data.channel(2).values;
% v_m = T_0_data.channel(1).values;
% v_r_0 = v_r(Ns:Ne)*0.001/60; %% mm/min --> m/sec
% v_m_0 = v_m(Ns:Ne)*0.001/60; %% mm/min --> m/sec
% 
% %Set Time Across The Range
% dt = 0.003;
% t = ([0:(Ne-Ns)]*dt)';
% omega = 1/55;
% 
% %Select The Same Range In The Other Datasets
% v_r1 = T_con_data.channel(2).values;
% [~,offset_con] = max(xcorr(v_r,v_r1));
% offset_con = length(v_r1) - offset_con;
% v_m1 = T_con_data.channel(1).values;
% v_r_con = v_r1(Ns+offset_con:Ne+offset_con)*0.001/60; %% mm/min --> m/sec
% v_m_con = v_m1(Ns+offset_con:Ne+offset_con)*0.001/60; %% mm/min --> m/sec
% 
% v_r2 = T_sweep_data.channel(2).values;
% [~,offset_sweep] = max(xcorr(v_r,v_r2));
% offset_sweep = length(v_r2) - offset_sweep;
% v_m2 = T_sweep_data.channel(1).values;
% v_r_sweep = v_r2(Ns+offset_sweep:Ne+offset_sweep)*0.001/60; %% mm/min --> m/sec
% v_m_sweep = v_m2(Ns+offset_sweep:Ne+offset_sweep)*0.001/60; %% mm/min --> m/sec
% 
% %Create data vecotors as per the method setup
% y_tr = v_m_con(2:end);
% phi_tr = [v_m_con(1:end-1), v_r_con(1:end-1), -1*pv(v_m_con(1:end-1), 0.0001), -1*nv(v_m_con(1:end-1), 0.0001)];
% t_tr = t(1:end-1);
% 
% y_te = v_m_sweep(2:end);
% phi_te = [v_m_sweep(1:end-1), v_r_sweep(1:end-1), -1*pv(v_m_sweep(1:end-1), 0.0001), -1*nv(v_m_sweep(1:end-1), 0.0001)];
% t_te = t(1:end-1);

%% Cross Validation
sig1 = linspace(-1,3,4); 
sig2 = linspace(-1.5,3,4);
gam = 150;%linspace(1,50,10); gam(3) = 9.5;
RMS_b = zeros(length(sig1), length(sig2), length(gam));

c = 1;

for a = 1:length(sig1)
    for b = 1:length(sig2)
        for d = 1:length(gam)
            c
            alpha = 10^-sig1(a);
            bravo = 10^-sig2(b);
            echo = gam(d);
            ohm = omega;
            [RMS_b(a,b,d),~] = CV_Cost([y_te,phi_te,t_te],[alpha, bravo, echo, ohm]);
            c = c+1;
        end
    end
end

[v,loc] = min(RMS_b(:));
[ii,jj,k] = ind2sub(size(RMS_b),loc);

X = 10.^-sig1;
Y = 10.^-sig2;
for d = 1:length(gam)
    figure(d);
    surf(Y,X,RMS_b(:,:,d))
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca,'zscale','log')
    set(gca, 'Colorscale','log')
    ylabel("\sigma_1")
    xlabel("\sigma_2")
    zlabel("Fitness of Parameters")
    title("Cross Validation of \sigma_1 and \sigma_2 for PL-LSSVM")
end