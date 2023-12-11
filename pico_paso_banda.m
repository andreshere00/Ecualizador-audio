function [Bp, Ap] = pico_paso_banda(dBgain, f0, BW)
    Fs = 48000;
    A = sqrt(10^(dBgain/20));
    w0 = 2*pi*f0/Fs;
    alpha = sin(w0)*sinh(log(2)/2 * BW * w0/sin(w0));
    
    b0 =   1 + alpha*A;
    b1 =  -2*cos(w0);
    b2 =   1 - alpha*A;
    a0 =   1 + alpha/A;
    a1 =  -2*cos(w0);
    a2 =   1 - alpha/A;
    
    Bp = [b0 b1 b2];
    Ap = [a0 a1 a2];
end
