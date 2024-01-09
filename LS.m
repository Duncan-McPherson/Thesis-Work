function [y, J, B, MAE, beta] = LS(Vm, Vr, dt, dB)
    %Least Squares Anaylsis
    phi = [Vm(1:end-1),Vr(1:end-1),-pv(Vm(1:end-1),dB),-nv(Vm(1:end-1),dB)];
    beta = (phi'*phi)\phi'*Vm(2:end);

    %Determine Prediction and Normalized Error
    y = phi*beta;
    MAE = sum(abs(y-Vm(2:end)))/length(y);

    %Find Interia and Vicouse Dampening
    pd = beta(1) + beta(2);
    J = ((pd-1)*dt)/(beta(2)*log(pd));
    B = (1-pd)/(beta(2));
end