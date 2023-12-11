function [Bl, Al] = repisa_paso_bajo(dBgain, f0, BW)
    Fs = 48000;
    A  = sqrt(10^(dBgain/20));
    w0 = 2*pi*f0/Fs;
    alpha = sin(w0)*sinh(log(2)/2 * BW * w0/sin(w0));
    
    b0 =    A*( (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha );
    b1 =  2*A*( (A-1) - (A+1)*cos(w0)                   );
    b2 =    A*( (A+1) - (A-1)*cos(w0) - 2*sqrt(A)*alpha );
    a0 =        (A+1) + (A-1)*cos(w0) + 2*sqrt(A)*alpha;
    a1 =   -2*( (A-1) + (A+1)*cos(w0)                   );
    a2 =        (A+1) + (A-1)*cos(w0) - 2*sqrt(A)*alpha;
    
    Bl = [b0 b1 b2];
    Al = [a0 a1 a2];
end