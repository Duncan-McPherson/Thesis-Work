function [y, y_nl, beta,c] = Duncan_KLS(y_ds, phi_ds, t, omega, sig1, sig2, gam)
    %Setup up regressor matrix
    jit = 1e-6;

    %Least Squares Anaylsis
    N = size(phi_ds,2); %dimension of phi
    L = length(t(1:end-1)); %dimension of t

    %Calculate
    Krn = Kernel(t(1:end-1), t(1:end-1), omega, sig1, sig2);
    Ohm = (Krn + gam*eye(L));
    Phi = [zeros(N), zeros(N,1), phi_ds'; zeros(1,N), 0, ones(1,L); phi_ds, ones(L,1), Ohm];
    Y = [zeros(N+1,1); y_ds];
    KLS_parm = (Phi'*Phi)\Phi'*Y;

    %Set Parameters
    alpha = KLS_parm(N+2:end);
    c = KLS_parm(N+1);
    beta = KLS_parm(1:N);

    %Determine Prediction and Normalized Error
    y = phi_ds*beta + (Krn+jit*eye(L))*alpha + c;
    y_nl = (Krn+jit*eye(L))*alpha + c;
end
    
    
    
    
  