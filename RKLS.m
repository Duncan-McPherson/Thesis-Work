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

function [Fc, Y] = RKLS(phi, y, t, M, beta, hp)
    %% Constants
    omega = hp(1); %Kernel HP, freqnuecy of disturbance
    sig1 = hp(2);  %Kernel HP #2
    sig2 = hp(3);  %Kernel HP #3
    gam = hp(4);   %Regularization Parameter
    lambda = hp(5);%Forgetting Parameter
    eta = hp(6);   %Error Parameter
    
    Xi = [phi,t];
    Yi = y;
    
    n = size(Xi,1);     %Rank of data
    c = size(Xi,2) - 1; %Col of data
    Y = zeros(n,1);     %Output prediction
    Fc = zeros(n,1);    %Force prediction
    
    kss = 1 + eta;                              %Covariance of first data
    mu = (Yi(1)-Xi(1,1:c)*beta)*kss/(gam+kss);  %mean alpha values
    Sigma = kss - kss^2/(gam+kss);              %alpha varaince
    Q = 1/kss;                                  %Inverted Kernel Matrix
    D = [1];                                    %Index of basis vectors
    m = 1;                                      %Size of model
    %% Loop through data and determine output
    for j = 1:n 
        X = Xi(j,:);
        Xl = X(1:c);
        Xb = Xi(D,c+1);
        y = Yi(j,:);
        
        % Forget a bit
        Sigma = lambda*Sigma + (1-lambda)*(Kernel(Xb,Xb,omega,sig1,sig2) + eta*eye(m));
        mu = sqrt(lambda)*mu;
        
        % Predict distrubance at sample t
        kbs = Kernel(Xb,X(c+1),omega,sig1,sig2); kss = 1 + eta;
        q = Q*kbs; 
        Fc(j) = q'*mu;
        Y(j) = Fc(j) + Xl*beta;
        h = Sigma*q;
        gammat = kss - q'*kbs; gammat(gammat<0)=0;
        sf2 = gammat + q'*h; sf2(sf2<0)=0;
        sy2 = gam + sf2;
        
        % Include sample t and add a new basis
        Qold = Q;
        p = [q;-1];
        Q = [Q zeros(m,1);zeros(1,m) 0] + 1/gammat*(p*p'); 
        
        p = [h;sf2];
        mu = [mu;Fc(j)] + ((y - Xl*beta - Fc(j))/sy2)*p;
        Sigma = [Sigma h; h' sf2] - 1/sy2*(p*p'); 
        D = [D;j]; m = m + 1; 
        
        %----- Delete a basis if necessary
        if m > M || gammat < eta
            if gammat < eta % To avoid roundoff error
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