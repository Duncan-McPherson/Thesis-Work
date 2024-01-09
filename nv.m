function y = nv(w, Ohm)
    w(abs(w) <= Ohm) = 0;
    y = sign(w);
    y = -1/2*y.*(ones(length(w),1)-y);
end