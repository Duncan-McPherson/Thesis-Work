clear; clc; close all;
%% Hyperparameters
%Good Result Parameters
gam = 150; %Regularization Factor 
sig1 = 0.075; % %Kernel Parameter: Sigma 1 for amp
sig2 = 1; % %Kernel Parameter: Sigma 2 for freq 
omega = 60/3000; %UBC/SIM period time [1/RevPerSec]
%omega = 1/55; %Munich 

ls_betas = [0;0];
kernel_betas = [0;0];
%% Load Data and Confirm Ranges for Each Data Selection
Range = load("UBC Data\Ranges.mat");
for i = [1,4,13]
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

%% Generate Espinoza
%     %As per Espinoza,
%     N = 1000;
%     t = 0:0.003:(N-1)*0.003;
%     noise = 0.001*rand(N,1);
%     u = rand(N,1);
%     y(1) = 0.3*u(1) + noise(1);
%     y(2) = 0.3*u(2) - 0.2*u(1) + noise(2);
%     for k = 3:N
%         y(k) = -0.5*y(k-2) + 0.3*u(k) - 0.2*u(k-1) + noise(k);
%     end
%     y1 = (y)';
%     switch i 
%         case 1
%             y_sin = zeros(1,N);
%         case 4
%             y_sin = (0.2*sin((2*pi()/omega)*t));
%         case 13
%             y_sin = (0.2*sin((2*pi/2.4)*t).*sin((2*pi()/omega)*t));
%         otherwise
%             y_sin = zeros(1,N);
%     end    
%     y2 = (y + y_sin)';
%     
%     %As per Espinoza, train on 400 of the data points
%     phi_ds = [y2(350:749),u(352:751),u(351:750)];
%     y_ds = y2(352:751);
%     t2 = t(352:751);
% 
%     phi_te = [y2(750:end-2),u(752:end),u(751:end-1)];
%     y_te = y2(752:end);
%     t_te = t(752:end);

    %% Load TUM Data
    %[phi_ds, y_ds, t2, yd] = loadMunich(i);
    
    %% Conduct LS Anaylsis on Downsampled Data
    beta = (phi_ds'*phi_ds)\phi_ds'*y_ds;
    ls_betas = [ls_betas,beta];

    ls_rms = rms(y_ds - phi_ds*beta)
    
    %% Conduct LS-SVM Anaylsis on Downsampled Data
    [y_kls, y_nl, kbeta, c, alpha] = KLS(y_ds, phi_ds, t2, omega, sig1, sig2, gam);
    kernel_betas = [kernel_betas,kbeta];
    
    Krn = Kernel(t2, t2, omega, sig1, sig2);
    kls_rms = rms(y_ds - (phi_ds*kbeta + Krn'*alpha + c))
    %f_rms = rms(y_sin(352:751)' - (Krn'*alpha + c))

    figure(i);
    %subplot(2,1,1);
    plot(t2, y_ds); hold on;
    plot(t2, y_kls); 
    ylabel("Velocity [m/s]"); xlim([t2(1) t2(end)]);
    xlabel("Time [s]"); legend("Measured Velocity","Predicted Kernel Output");
%     subplot(2,1,2);
%     plot(t2, y_sin(352:751)'); hold on;
%     plot(t2, y_nl); 
%     ylabel("Simulated Output []"); xlim([t2(1) t2(end)]);
%     xlabel("Time [s]"); legend("Simulated Disturbance","Predicted Kernel Disturbance");

%% Kernel Recursive Least Squares on Downsampled Data
%     hp = [omega, 1, 1, gam, 0.999, 1e-6];
%     [Fc, Y] = RKLS(phi_ds, y_ds, t2', 75, [-0.5;0.3;-0.2], hp);
    
%% Cutting Force Reconstruction
    %LS
    pd_ls = beta(1)+beta(2);
    p_ls = log(pd_ls)/(dt*dS);
    Ferror = (p_ls/(1-pd_ls))*(y_ds - (phi_ds*beta));
    
    %KLS
    pd_kls = kbeta(1)+kbeta(2);
    p_kls = log(pd_kls)/(dt*dS);
    Fkls = (p_kls/(1-pd_kls))*y_nl;
    
%     %KRLS
%     Fkrls = (p_kls/(1-pd_kls))*Fc;
    
    %Detrend Force Data from Dyano
    Fx_nl = detrend(Fx_ds(1:end-1));

%% Plot Results
%     figure(i+1);
%     yyaxis left; 
%     plot(t2, Ferror); ylim([-2e-1 2e-1]); ylabel("Mass Normalized Force (LS) [m/s^2]");
%     %xlim([10 10.5]);
%     yyaxis right; 
%     plot(t2, Fx_nl); ylim([-200 200]); ylabel("Force [N]");
%     %xlim([10 10.5]);
%     xlabel("Time [s]"); title("Least Squares Error (y-\phi\beta) vs. Cutting Force ({F_c})");
%     
    figure(i+2);
    yyaxis left;
    plot(t2, Fkls);  ylabel("Mass Normalized Force (Kernel) [m/s^2]");
    ylim([-0.2 0.2]); xlim([0 t2(end)]);
    yyaxis right;
    plot(t2, Fx_nl);  ylabel("Force [N]");
    ylim([-200 200]); xlim([0 t2(end)]);
    xlabel("Time [s]"); title("Kernel Predicted Force (\alphaK+c) vs. Cutting Force ({F_c})");
%     
%     figure(i+3);
%     yyaxis left; 
%     plot(t2(1:end-1), Fkrls);  ylabel("Mass Normalized Force (Recursive Kernel) [m/s^2]");
%     ylim([-0.2 0.2]); xlim([0 t2(end)]);
%     yyaxis right; 
%     plot(t2(1:end-1), Fx_nl);  ylabel("Force [N]");
%     ylim([-200 200]); xlim([0 t2(end)]);
%     xlabel("Time [s]"); title("Recursive Kernel Predicted Force vs. Cutting Force ({F_c})");
%     
%     figure(i+4);
%     yyaxis left; 
%     plot(t2(1:end-1), Fkrls);  ylim([-2e-1 2e-1]); ylabel("Mass Normalized Force (Kernel) [m/s^2]");
%     %xlim([10 10.5]);
%     yyaxis right; 
%     plot(t2(1:end-1), Ferror);  ylim([-2e-1 2e-1]); ylabel("Mass Normalized Force (LS) [m/s^2]");
%     %xlim([10 10.5]);
%     xlabel("Time [s]"); title("Recursive Kernel Predicted Force vs. Least Squares Error (y-\phi\beta)");
end