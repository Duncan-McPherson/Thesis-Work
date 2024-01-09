clear; clc; close all;
Zero = load("Analysis\T_1_Dec_1_DS30.mat");

load("UBC Data\Ranges.mat")
for j = [4,13,14,1]
    load("Analysis\T_"+j+"_Dec_1_DS30.mat");
    load("UBC Data\Data_Trial_"+j+".mat");
    j
    diff_LS = 100*abs((Zero.beta-beta)./Zero.beta)
    diff_KLS = 100*abs((Zero.kbeta-kbeta)./Zero.kbeta)

    t1 = t(Ns(j):Ne(j))-t(Ns(j));
    Vm = VelocityXAct(Ns(j):Ne(j))*(1e-3/60); %Convert mm/min to m/s

    dt = t1(2)-t1(1); %Sampling time
    dS = 30; %Downsampling Rate
    N = 10000; %Sampling window
    fc = 1/(2*dt*dS); %Cutoff Freq
    fs = 1/dt; %Sampling Freq
    [b,a] = butter(6,fc/(fs/2));

    Vm_ds = downsample(filtfilt(b,a,Vm),dS);
    t2 = downsample(t1,dS);

    figure(j)
    plot(t2(2:end), Vm_ds(2:end)); hold on
    %plot(t2(2:end), -y_ls);
%     plot(t2(2:end), y_nl);
    plot(t2(2:end), y_kls);
    xlim([17.9 18.35])
%     ylim([-7e-4 7e-4])
%     legend("LS Error to Measured", "Nonlinear Components of KLS")
    legend("Measured","Kernel Least Squares (KLS)")
    title("Test #"+j+":Velocity Output"); xlabel("Time [s]");ylabel("Velocity [m/s]")
end