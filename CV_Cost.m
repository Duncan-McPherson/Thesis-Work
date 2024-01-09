function [RMSE_batch, RMSE_rec, RMSE_Cross] = CV_Cost(Test_Data1,theta)
%% Hyperparameters:
sig1 = theta(1);
sig2 = theta(2);
gam = theta(3);
omega = theta(4);

%% Seperate Data into training and test (LS-SVM)
%Parameter Identification Using Least Squares on Test 4 (Cutting)
y = Test_Data1(:,1);
Z = Test_Data1(:,2:end-1);
t1 = Test_Data1(:,end);

L = length(t1); %Dimension of t
L_tr = round(0.8*L);
R = randperm(L,L_tr);

y_tr = y(R); 
Z_tr = Z(R,:);
t_tr = t1(R);

y(R) = [];
Z(R,:) = [];
t1(R) = [];

y_te = y; 
Z_te = Z;
t_te = t1;

%% Partially Linear LS-SVM
%Regression on training data
N = size(Z,2); %dimension of phi
Krn_tr = Kernel(t_tr, t_tr, omega, sig1, sig2);
Ohm = (Krn_tr + gam^-1*eye(L_tr));
Phi = [Ohm,ones(L_tr,1),Z_tr;ones(1,L_tr),0,zeros(1,N);Z_tr',zeros(N,1),zeros(N)];
Y = [y_tr;zeros(N+1,1)];
KLS_parm = (Phi'*Phi)\Phi'*Y;

%These are the parameters of the model
alpha_kls = KLS_parm(1:end-N-1);
c_kls = KLS_parm(end-N);
beta_kls = KLS_parm(end-N+1:end);

%Make predictions with the test data
Krn_te = Kernel(t_tr, t_te, omega, sig1, sig2);

y_hat = Z_te*beta_kls + Krn_te'*alpha_kls + c_kls;

error_PL = rms(y_te - y_hat);

%% Fitness Test
%Prediction Error of LS-SVM on its training data
RMSE_batch = error_PL;

%Empty Variable
RMSE_rec =  0;

%Empty Variable
RMSE_Cross = 0; 
end