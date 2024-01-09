function KMN = Kernel(XN, XM, omega, sig1, sig2)
    %Squared Exponential Kernel (A.K.A Gaussian)
    Num = pdist2(XN(:),XM(:)).^2;
    Dem = 2*(sig1^2);
    KMN1 = exp(-Num/Dem); 

    %Periodic Kernel
    Num = (2*sin(pi.*pdist2(XN(:),XM(:))./omega).^2);
    Dem = (sig2^2);
    KMN2 = exp(-Num/Dem); 

    %Combine to make the Locally Periodic Kernel, which has a frequency but
    %can change amplitudes over time.
    KMN = KMN1.*KMN2;
end