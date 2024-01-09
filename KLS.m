function [y, y_nl, beta, c, alpha] = KLS(y, Z, t, omega, sig1, sig2, gam)
    %Least Squares Anaylsis
    N = size(Z,2); %dimension of phi
    L = length(t); %dimension of t

    %Calculate
    Krn = Kernel(t, t, omega, sig1, sig2);
    Ohm = (Krn + gam^-1*eye(L));
    Phi = [Ohm,ones(L,1),Z;ones(1,L),0,zeros(1,N);Z',zeros(N,1),zeros(N)];
    Y = [y;zeros(N+1,1)];
    KLS_parm = (Phi'*Phi)\Phi'*Y;

    %Set Parameters
    alpha = KLS_parm(1:end-N-1);
    c = KLS_parm(end-N);
    beta = KLS_parm(end-N+1:end);

    %Determine Prediction and Normalized Error
    y = Z*beta + Krn*alpha + c;
    y_nl = Krn*alpha + c;
end
    
    
    
    
  