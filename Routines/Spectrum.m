function [freq,amp] = Spectrum(x,fs);
% Donnes l'analyse FFT du signal x echantillonne a la frequence fs
% freq frequences
% amp amplitudes



% TESTS -------------------------------------------------------------------


test = false;

if test
    % Amplitudes a retrouver : continnu , fondamental, harmoniques
    a = [11:-1:0]
    % Phases bidon pour verifier invariance
    phi = 2*pi*rand(size(a));
    % Frequence fondamentale a retrouver
    f = 440;  % Hz
    % Freq echantillonnage
    fs = 44000; % hz
    % Nombre de cycles sur la fondamentale
    Nc = 100;
    % Signal temporel
    t = [0:1/fs:Nc/f];
    x = a(1)*ones(size(t));
    for p=2:size(a,2)
        x = x + a(p)*cos(2*pi*(p-1)*f*t + phi(p));
    end
end


% FFT ---------------------------------------------------------------------


% "Longueur" signal
L = size(x,2)-1;
y = fft(x);
% On conserve le seul premier cote de la FFT
P2 = abs(y/L);
amp = P2(round(1:L/2+1));
amp(2:end-1) = 2*amp(2:end-1);
% Frequences 
freq = fs*(0:(L/2))/L;


% SORTIES -----------------------------------------------------------------


if test
    figure(1);clf
    plot(t,x)
    xlabel t
    ylabel y

    figure(2);clf
    bar(freq,amp,1)
    xlabel frequency(Hz)
    ylabel intensity
    %sound(x,fs)
    disp(amp(1:size(a,2)));
    disp(freq(1:size(a,2)));
end