Tr_set = load('UBC Data/Data_Trial_4.mat');
Te_set = load('UBC Data/Data_Trial_13.mat');

%Set Range Around Time When Cutting Occured
Ranges = load('UBC Data/Ranges.mat');
Ns4 = Ranges.Ns(3);
Ne4 = Ranges.Ne(3);
Ns13 = Ranges.Ns(13);
Ne13 = Ranges.Ne(13);

%Set time for each data set
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

%% Create data vecotors for PL-LSSVM method
beta = [0.3473;0.6527];

y_observed1 = vm_tr(2:end);  % output
x_observed1 = [vm_tr(1:end-1), vr_tr(1:end-1)]; % regressor array
X1 = t_tr(1:end-1);
Y1 = y_observed1-x_observed1*beta;

y_observed2 = vm_te(2:end);  % output
x_observed2 = [vm_te(1:end-1), vr_te(1:end-1)]; % regressor array
X2 = t_te(1:end-1);
Y2 = y_observed2-x_observed2*beta;

%% GPR

%Custom fitrgp
kfcn = @(XN,XM,theta) exp(-(pdist2(XN,XM).^2)/(2*(sig1^2)))*exp(-(sin(pi.*pdist2(XN,XM)./omega).^2)/(2*(sig2^2)));
hp = [1e-6,0.46,60/3000];

gprMdl1 = fitrgp(X1,Y1,'KernelFunction',kfcn,'KernelParameters',hp,'SigmaLowerBound',1e-10);
%gprMdl2 = fitrgp(x_observed2,y_observed2);

[ypred1,~,yint1] = predict(gprMdl1,X2);
%[ypred2,~,yint2] = predict(gprMdl2,x_observed3);

fig = figure;
fig.Position(3) = fig.Position(3)*2;

hold on
scatter(t_te(2:end),Y2)  % Observed data points
plot(t_te(2:end),ypred1,'g')                  % GPR predictions
patch([t_te(2:end);flipud(t_te(2:end))],[yint1(:,1);flipud(yint1(:,2))],'k','FaceAlpha',0.1); % Prediction intervals
hold off
title('GPR Fit of Noise-Free Observations')
%legend({'Noise-free observations','g(x) = x*sin(x)','GPR predictions','95% prediction intervals'},'Location','best')
