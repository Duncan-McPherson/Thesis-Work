clear; clc; close all;

rng('default')
rng(1);

%Rec
sig1 = 0.717;
sig2 = 1;
gam = 10;
omega = 60/3000;
hyp = [omega, sig1, sig2, gam];
%batch
sig1 = 0.077;
sig2 = 1.6;
gam = 10;

kaf_param.eta = 0.01; % learn rate
kaf_param.M = 20;
kaf_param.kerneltype = 'poly';
kaf_param.kernelpar = hyp;
kaf_class = @sknlms;

%% UBC
Range = load("UBC Data\Ranges.mat");
i = 1;
%% Load MAL Data
name = "UBC Data\Data_Trial_"+i+".mat";
load(name);

%Set Range Around Time When Cutting Occured
Ns = Range.Ns(i);%+15000;
Ne = Range.Ne(i);%-5000;

%Collect And Store Data From Measurements
t1 = t(Ns:Ne)-t(Ns); % [sec]
Vm = VelocityXAct(Ns:Ne)*(1e-3/60); %Convert to [m/s]
Vr = VelocityXNom(Ns:Ne)*(1e-3/60);
Fnl = Fx(Ns:Ne); % [Newtons]

%Downsampling
dt = t1(end)/(length(t1)-1); %Sampling time
dS = 30; %Downsampling Rate
N = 10000; %Sampling window
fc = 1/(2*dt*dS); %Cutoff Freq
fs = 1/dt; %Sampling Freq
[b,a] = butter(6,fc/(fs/2));

%Apply Downsample to Data
Vr_ds = downsample(filtfilt(b,a,Vr),dS);
Vm_ds = downsample(filtfilt(b,a,Vm),dS);
Fx_ds = downsample(filtfilt(b,a,Fnl),dS);
t2 = downsample(t1,dS); t2 = t2(1:end-1);

%Make Regressor and Output Vectors
phi_ds = [Vm_ds(1:end-1),Vr_ds(1:end-1)];
y_ds = Vm_ds(2:end);

%% Espinoza
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
% i = 13;
% switch i
%     case 1
%         y_di = 0;
%     case 4
%         y_di = (0.3*sin((2*pi()/omega)*t));
%     case 13
%         y_di = (0.3*sin(2*pi/2.5*t).*sin((2*pi()/omega)*t));
%     otherwise
%         y_di = 0;
% end
% y2 = (y + y_di)';
% 
% %As per Espinoza, discard first 350 data
% phi_ds = [y2(350:749),u(352:751),u(351:750)];
% y_ds = y2(352:751);
% t2 = t(352:751)';
    

tic
[y_rec, y_nl_rec, rbeta] = SKRLST(phi_ds, t2, y_ds, 0, 80, kaf_param, 0.9);
                            %regressor, time, outptu, lin buff, nonlin
                            %buff, kernel hyp, rls learning/forgetting 0.99995
toc

rbeta

% tic
% [y_batch, y_nl_batch, kbeta, c, alpha] = KLS(y_ds, phi_ds, t2, omega, sig1, sig2, gam);
% toc
% kbeta
beta = (phi_ds'*phi_ds)\phi_ds'*y_ds

rms(y_ds - y_rec)
%rms(y_ds - y_batch);

figure(i);
subplot(2,1,1);
plot(t2, y_ds); hold on;
plot(t2, phi_ds*beta); 
plot(t2, y_rec); 
ylabel("Velocity [m/s]"); xlim([t2(1) t2(end)]);
xlabel("Time [s]"); legend("Measured Velocity","Predicted Kernel Output");
subplot(2,1,2);
plot(t2, y_ds - phi_ds*[0.3473;0.6527]); hold on;
plot(t2, y_nl_rec); 
ylabel("Simulated Output []"); xlim([t2(1) t2(end)]);
xlabel("Time [s]");