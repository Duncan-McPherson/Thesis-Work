%  The RSPKLS - Recursive Semi-Parametric Kernel Least Squares
%
%  Alternate goal of predicting force given all hyperparameters including
%  pyhsics parameters
%
%  [Fc, Y] = RKLS(Vm, Vr, M, beta, hp, lambda, gam, del)
%    IN
%    Vm, Vr:      Output and Input identification data respectivly
%    M:           Max size of the kernel model
%    beta:        Pyhsical Parameters of the System
%    lambda:      Forgetting factor                                             
%    gam:         Regularization factor
%    del:         Error factor/jitter
%
%    OUT
%    Fc:          Force Predictions
%    Y:           Velocity Outputs
%
%   For further details, send email to jdmcpher@uvic.ca
%   2023 - Duncan McPherson

function [Ynl, Yl, beta_log] = Split_RKLS(phi, y, t, M, hp)
    %% Setup before Run Time
    %Hyperparameters
    omega = hp(1); %Kernel HP, freqnuecy of disturbance
    sig1 = hp(2);  %Kernel HP #2
    sig2 = hp(3);  %Kernel HP #3
    gam = hp(4);   %Regularization Parameter
    lambda = hp(5);%Forgetting Parameter
    eta = hp(6);   %Error Parameter
    
    %Data and Final Output Variables
    Xi = [phi,t(1:end-1)];  %Regressor with Kernel Input
    n = size(Xi,1); c = size(Xi,2) - 1;         %Rank and Column of Regressor
    Yi = y;                             %Output Vector
    
    Yl = zeros(n,1);     %Output Prediction
    Ynl = zeros(n,1);    %Force Prediction
    beta = [0.3684;0.6316];%zeros(c,1);   %Parameter Prediction
    beta_log = zeros(c,n+1);
    
    %Variables for Linear Portion (The egg starts here)
    P = 6*[10, 0; 0, 18.8];                           %Covariance Matrix for beta
    L = P*Xi(1,1:c)'/(lambda + Xi(1,1:c)*P*Xi(1,1:c)'); %Calculate Linear Gain Matrix
    beta = beta + L*(Yi(1) - Xi(1,1:c)*beta);           %Update Beta Prediction
    P = (eye(2) - L*Xi(1,1:c))*(1/lambda)*P; 
    beta_log(:,1) = beta;
    
    %Variables for Nonlinear Portion (The chicken is then comes next)
    kss = 1 + eta;                              %Covariance of first data
    mu = (Yi(1)-Xi(1,1:c)*beta)*kss/(gam+kss);  %Mean alpha values
    Sigma = kss - kss^2/(gam+kss);              %Alpha covaraince
    Q = 1/kss;                                  %Inverted Kernel Matrix
    D = [1]; %Index of basis vectors
    m = 1;   %Size of model
    
    %% Loop through data and determine output
    for j = 1:n
        %% New Data
        Xnl = Xi(j,c+1);%Input for Nonlinear from j time
        Xl = Xi(j,1:c); %Input for linear from j time
        Xb = Xi(D,c+1); %Input for Nonlinear from Dictionary
        y = Yi(j,:); %Output from the j+1 time step
    
        % Forget a bit of the Nonlinear Prediction
        Sigma = lambda*Sigma + (1-lambda)*(Kernel(Xb,Xb,omega,sig1,sig2) + eta*eye(m));
        mu = sqrt(lambda)*mu;
       
        %% Prediction (This happens before y is discovered)
        % Predict NonLinear element using KRLS-T
        kbs = Kernel(Xb,Xnl,omega,sig1,sig2); kss = 1 + eta;
        q = Q*kbs; 
        Ynl(j) = q'*mu;
        % Update variances based on prediction
        h = Sigma*q;
        gammat = kss - q'*kbs; gammat(gammat<0)=0;
        sf2 = gammat + q'*h; sf2(sf2<0)=0;
        sy2 = gam + sf2;
    
        %Predict Linear element using RLS
        Yl(j) = Xl*beta;
    
        %% Update (This happens ater y is discovered)
        % Update Linear Portion using (Y - Ynl)
        error_lin = y - Ynl(j);                 
        L = P*Xl'/(lambda + Xl*P*Xl');          %Calculate Linear Gain Matrix
        beta = beta + L*(error_lin - Xl*beta); %Update Beta Prediction
        P = (eye(2) - L*Xl)*(1/lambda)*P;       %Covariance Matrix for beta
        %(1/lamd)*(P - (P*phi'*phi*P)/(lamd+phi*P*phi')); EQ from Indian Paper
        beta_log(:,j+1) = beta;
    
        % Update NonLinear Portion using (Y - Yl)
        %Add the new data to the model
        Qold = Q;
        p = [q;-1];
        Q = [Q zeros(m,1);zeros(1,m) 0] + 1/gammat*(p*p');
      
        p = [h;sf2];
        error_nlin = y - Yl(j);
        mu = [mu;Ynl(j)] + ((error_nlin - Ynl(j))/sy2)*p;
        Sigma = [Sigma h; h' sf2] - 1/sy2*(p*p'); 
        D = [D;j]; m = m + 1; 
        
        %% Dictionary Management
        if m > M || gammat < eta
            if gammat < eta % To avoid roundoff error
                if gammat < eta/10
                    %warning('Numerical roundoff error too high, you should increase jitter noise')
                end
                criterium = [ones(1,m-1) 0];
            else
                errors = (Q*mu)./diag(Q);
                criterium = abs(errors);
            end
            % Remove element r, which incurs in the minimum error 
            [~, r] = min(criterium);            
            small = 1:m; small(r) = [];
    
            if  r == m   % If we must remove the element we just added, perform reduced update instead
                Q = Qold;
            else
                Qs = Q(small, r); qs = Q(r,r); Q = Q(small, small);
                Q = Q - (Qs*Qs')/qs; 
            end
            %Reduce other components
            mu = mu(small);
            Sigma = Sigma(small, small);
            D = D(small); m = m - 1;
        end
    end
end