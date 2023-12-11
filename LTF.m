%
% TDSÑ-G35
% Laboratorio: Trabajo final
%
% alumno1: a.herencia@alumnos.upm.es
% alumno2: jose.calderon.villar@alumnos.upm.es
%
%% Un poco de limpieza
clear; clc; close all;

%% Sistema a ecualizar Sala-sistema distorsionador

[Z, P, K] = room('a.herencia@alumnos.upm.es', 'jose.calderon.villar@alumnos.upm.es');
[B, A] = zp2tf(Z, P, K);
[H, w] = freqz(B, A, 4096);
mH = 20*log10(abs(H));
subplot(1,2,2);
plot(w, mH);
title('Modulo de la respuesta en frecuencia de la sala');
xlabel('\omega'); ylabel('dB'); grid on;
axis([0, pi, min(mH)-1, max(mH)+1]);
subplot(1,2,1);
zplane(Z,P);
title('Diagrama de polos y ceros de la distorsión de la sala')

%% Compensador perfecto
% Hmin(z)
    Zmin = Z;
    Pmin = P(abs(P)<=1); % Aquellos polos que estén dentro de la circunferencia unidad 
    % los conservamos
    z0 = Z(abs(Z)>1); % Cogemos los ceros que estén fuera de la circunferencia unidad
    zrc = 1/conj(z0); % Calculamos su recíproco conjugado
    Zmin(Zmin==z0)=zrc; % Y sustituimos los nuevos ceros dentro de la c.u. por los que había fuera.
    figure(1);
    subplot(2,2,1);
    zplane(Zmin, Pmin);
    title('Diagrama de ceros y polos H_{min}(\omega)');

% Hap(z)
    Zap = Z(abs(Z)>1);
    Pap = 1/conj(Zap);
    Kap = (1+z0)/(1+zrc); % Si Hap es real, como es el caso, 
    % podemos calcular su constante de esta manera.
    subplot(2,2,2);
    zplane(Zap, Pap);
    title('Diagrama de ceros y polos H_{ap}(\omega)');

% Hc(z)
    Zc = Pmin;
    Pc = Zmin;
    Kmin = 1 / Kap;
    subplot(2,2,3);
    zplane(Zc, Pc);
    title('Diagrama de ceros y polos H_{c}(\omega)');
    [Bc, Ac] = zp2tf(Zc, Pc,Kmin); 
    [Hc,wc] = freqz(Bc,Ac,4096);
    mHc = 20*log10(abs(Hc));
    subplot(2,2,4);
    plot(w, mHc);
    title('Modulo de Hc(w)');
    xlabel('\omega'); ylabel('dB'); grid on;
    axis([0, pi, min(mHc)-1, max(mHc)+1]);

%% Ecualizador parametrico discreto

%Filtro repisa paso bajo
    [Bl, Al] = repisa_paso_bajo(9.2,1600,2);
    [Hl, wl] = freqz(Bl, Al, 4096);
    mHl = 20*log10(abs(Hl));
    figure(2);
    subplot(3,2,1);
    plot(wl, mHl);
    title('Módulo de la respuesta en frecuencia de H_{l}(z)');
    xlabel('\omega'); ylabel('dB'); grid on;
    axis([0, pi, min(mHl)-1, max(mHl)+1]);

%Filtro repisa paso alto
    [Bh, Ah] = repisa_paso_alto(-1.1,20000,0.3);
    [Hh, wh] = freqz(Bh, Ah, 4096);
    mHh = 20*log10(abs(Hh));
    subplot(3,2,2);
    plot(wh, mHh);
    title('Módulo de la respuesta en frecuencia de H_{h}(z)');
    xlabel('\omega'); ylabel('dB'); grid on;
    axis([0, pi, min(mHh)-1, max(mHh)+1]);

%Filtro pico paso banda 1
    [Bp1, Ap1]=pico_paso_banda(5.6,3150, 2.2);
    [Hp1, wp1]= freqz(Bp1, Ap1, 4096);
    mHp1= 20*log10(abs(Hp1));
    subplot(3,2,3);
    plot(wp1, mHp1);
    title('Módulo de la respuesta en frecuencia de H_{p1}(\omega)');
    xlabel('\omega'); ylabel('dB'); grid on;
    axis([0, pi, min(mHp1)-1, max(mHp1)+1]);

%Filtro pico paso banda 2
    [Bp2, Ap2]=pico_paso_banda(-1.6,10000,0.8);
    [Hp2, wp2]= freqz(Bp2, Ap2, 4096);
    mHp2= 20*log10(abs(Hp2));
    subplot(3,2,4);
    plot(wp2, mHp2);
    title('Módulo de la respuesta en frecuencia de H_{p2}(\omega)');
    xlabel('\omega'); ylabel('dB'); grid on;
    axis([0, pi, min(mHp2)-1, max(mHp2)+1]);

%Filtro pico paso banda 3
    [Bp3, Ap3]=pico_paso_banda(-8.6,16000,1.1);
    [Hp3, wp3]= freqz(Bp3, Ap3, 4096);
    mHp3= 20*log10(abs(Hp3));
    subplot(3,2,5);
    plot(wp3, mHp3);
    title('Módulo de la respuesta en frecuencia de H_{p3}(\omega)');
    xlabel('\omega'); ylabel('dB'); grid on;
    axis([0, pi, min(mHp3)-1, max(mHp3)+1]);

%Filtro ecualizador parametrico discreto
    He = Hl  .* Hp1 .* Hp2 .*Hp3 .* Hh;
    mHe = 20*log10(abs(He));
    subplot(3,2,6);
    plot(wp3,mHe);
    title('Módulo de la respuesta en frecuencia de H_{e}(\omega)');
    xlabel('\omega'); ylabel('dB'); grid on;
    axis([0, pi, min(mHe)-1, max(mHe)+1]);

%% Filtro final

% Obtención del sistema ecualizador final, multiplicando el sistema
% distorsionador por el ecualizador paramétrico. Representación y
% caracterización.

Heq = H .* He;
mHeq = 20 * log10(abs(Heq));
figure(4);
plot(w, mHeq+1.5);
title('Módulo de la respuesta en frecuencia de H_{eq}(\omega)');
xlabel('\omega'); ylabel('dB'); grid on;
axis([0, pi, min(mHeq)-3, max(mHeq)+5]);

figure(3)
plot(wp3, mHe);
hold on
grid on;
plot(w, mHc+0.5);
plot(w, mHc-0.5);
legend('H_{e}(\omega)', 'Límite superior', 'Límite inferior')
hold off

%% Señal de prueba chirp y montaje del sistema final

fs = 48000;
f0 = 20;
f1 = 20000;
t1 = 10;
t = 0:1/fs:t1;
x = chirp(t, f0, t1, f1,'logarithmic');

% Pasamos la señal chirp por cada uno de los filtros caracterizados por sus
% parámetros A y B.

x1 = filter(Bl, Al, x); 
x2 = filter(Bp1, Ap1, x1);
x3 = filter(Bp2, Ap2, x2);
x4 = filter(Bp3, Ap3, x3);
% Señal procedente de la salida de todos los filtros del sistema
% ecualizador
y = filter(Bh, Ah, x4);

% Señal compensada, tras ser pasada por el sistema distorsionador
xc = filter(B,A,y);

% Representación de cada señal, x[n], y[n], xc[n]
figure(1)
    subplot(3,1,1)
    semilogx(x);
        title('Señal de entrada x[n]');
        xlabel('Frecuencia (Hz)'); ylabel('Amplitud') % En escala logarítmica con semilog. Amplitud en unidades lineales
    subplot(3,1,2)
    semilogx(y);
        title('Señal distorsionada y[n]'); 
        xlabel('Frecuencia (Hz)'); ylabel('Amplitud')
    subplot (3,1,3)
    semilogx(xc);
        title('Señal compensada x_{c}[n]');
        xlabel('Frecuencia (Hz)'); ylabel('Amplitud')

fprintf("Nota: para escuchar las señales, se puede usar la función 'soundsc', \n especificando la frecuencia de muestreo f_{s}\n")