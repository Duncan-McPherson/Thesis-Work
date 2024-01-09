% [f, Gy] = fftspec(y, N, dt)
%Estimates the spectral density matrix as a function of frequency using the
%Welch averaging technique with a Hanning window and 50 % overlap.
%y:     Signal to be processed, measurement channels arranged in the rows    
%N:     Number of data points in the data segments (window size)
%dt:    Sampling time step
%f:     Vector containing the positive frequency axis
%Gy:   Spectral density matrices for the positive frequency axis
%The function is using the data segment routine fftseg.m. If no output
%arguments are used in the call, then the singular values of the spectral
%density matrices are plotted - or the auto spectral density in case of
%uni-variate signal.
%Example plotting the singular values of the Heritage court data set 1:
%   load dherita1.asc
%   y = dherita1';
%   fftspec(y, 1024, 0.025);

%All rights reserved, Rune Brincker, May 2012.

function [f, Gyy] = fftspec(y, N, dt)
[Y, df] = fftseg(y, N, dt);
nf = N/2+1;
f = (0:nf-1)'*df;
[nr,nc,ns] = size(Y);
Gyy = zeros(nc,nc,nf);
for k=1:nf;
    for s = 1:ns,
        Y1 = reshape(Y(k,:,s), nc, 1);
        Gyy(:,:,k) = Gyy(:,:,k) + conj(Y1)*Y1.';
    end
end
W = [0;hanning(N)];
Gyy = Gyy/(ns*W'*W);
if (nc == 1),
    Gyy = reshape(Gyy, 1, nf);
end
Gs = zeros(nf, nc);
if (nargout == 0)
    if (nc == 1),
        plot(f, 10*log10(Gyy'), 'k', 'LineWidth',2)
        title('Auto spectral density')
    else,
        for k=1:nf,
            [U,S,V] = svd(Gyy(:,:,k));
            Gs(k, :) = diag(S)';
        end
        plot(f, 10*log10(Gs), 'k', 'LineWidth',2)
        title('Singular values of spectral matrix')
    end
    ylabel('dB rel. to unit')
    xlabel('Frequency [Hz]')
    figure(gcf)
    grid
end


