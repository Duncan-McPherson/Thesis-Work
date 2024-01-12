function [y_hat, y_nl, beta] = SKRLST(phi, x, y, M_lin, M_nl, kpar, mu_L)
% Split Kernel Recursive Least Squares Tracker.
%
% The purpose of this function is to identify the parameters of a linear
% system under the effect of nonlinear disturbance using recursive kernel based methods
% 
% Inputs:   - phi: input, regressor matrix
%           - x: input, kernel input
%           - y: output of system
%           - M_lin: length of linear filter
%           - M_nl: memory of nonlinear (kernel) filter
%           - kaf_param: hyperparameters for kernel
%           - mu_L: NLMS learning rate
% Outputs:  - y_nl:estimated disturbance
%           - y: estimated output signal
%           - beta: linear parameters
%
% Duncan McPherson, 2024.
    
    L = size(phi,1); %Number of samples
    C = size(phi,2);

    y_hat = zeros(L,1); %Output signal initialization
    y_nl = zeros(L,1); %Disturbance signal initialization
    beta = zeros(C,1); %Parameters initialization

    P = (10^6)*eye(C); %Covariance Matrix for beta
    
    %Initialize the Recursive Kernel
    eta = 1e-2; % learn rate
    M = M_nl; % max dictionary size
    lambda = 1; % forgetting (0 all, 1 none)
    dict = []; % dictionary
    m = 0;

    %Kernel Parm
    sig1 = kpar.kernelpar(2);
    sig2 = kpar.kernelpar(3);
    gam = (1/kpar.kernelpar(4));
    omega = kpar.kernelpar(1);

    for k = 1:L
        if ~mod(k,floor(L/30)), fprintf('.'); end %Progress bar of 30 dots
               
        % compute linear output
        y_lin = phi(k,:)*beta; 
        
        % compute nonlinear output
        if m > 0
            %Forget part of the kernel
            Sigma = lambda*Sigma + (1-lambda)*(Kernel(dict,dict,omega,sig1,sig2) + eta*eye(m));
            mu = sqrt(lambda)*mu;

            kbs = Kernel(dict,x(k),omega,sig1,sig2); kss = 1 + eta;
            q = Q*kbs; 
            y_nl(k) = q'*mu;
            h = Sigma*q;
            gammat = kss - q'*kbs; gammat(gammat<0)=0;
            sf2 = gammat + q'*h; sf2(sf2<0)=0;
            sy2 = gam + sf2;
        else
            y_nl(k) = 0;
        end
        
        % compute overall output
        y_hat(k) = y_lin + y_nl(k); % linear combination of individual outputs
        d = y(k) - y_lin; %disturbance error
        e = y(k) - y_hat(k); %total error
        
        %Update Linear Parameters
        Lg = P*phi(k,:)'/(mu_L + phi(k,:)*P*phi(k,:)'); %Calculate Linear Gain Matrix
        P = (eye(C) - Lg*phi(k,:))*(1/mu_L)*P; 
        beta = beta + Lg*(e); %Update Beta Prediction
        
        %Update Kernel Dictionary (only after linear parameters are found) 
        if k > M_lin
            if m == 0 % initialize
                dict = x(k); m = m + 1;
                kss = 1 + eta;                             %Covariance of first data
                mu = (d)*kss/(gam + kss);            %Mean alpha values
                Sigma = kss - kss^2/(gam + kss);     %Alpha varaince
                Q = 1/kss;                                 %Inverted Krn
            else
                Qold = Q;
                p = [q;-1];
                Q = [Q zeros(m,1);zeros(1,m) 0] + 1/gammat*(p*p'); % O(n^2)
                
                p = [h;sf2];
                mu = [mu;y_nl(k)] + ((d - y_nl(k))/sy2)*p;
                Sigma = [Sigma h; h' sf2] - 1/sy2*(p*p'); 
                dict = [dict;x(k)]; m = m + 1;

                if m > M || gammat < eta
                    if gammat < eta % To avoid roundoff error
                        if gammat < eta/10 end
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
                    dict = dict(small); m = m - 1;
                end
            end
        end
    end
    fprintf('\n')
end