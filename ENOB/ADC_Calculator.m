%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for calculating Effective Number of Bit for an ADC
% Uses FFT test algorithm and calulates the signal to noise ratio
%
%   Dexter Elmendorf 10/02/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

%% Load data, calculate Fs
file = readmatrix("SARquantizer.dat");
x = transpose(file);
Ns = 300; % Number of samples
x = x(1:Ns);
N = length(x); % Full length of sample
fs = 100e3; % Sample Rate 
Ts = 1/fs; % Set the expected sample rate (find way to calculate from file)

time = linspace(0,N*Ts,N);

Vref = 3.3;
% voltage = voltage -(1024/2); % Remove the DC component

%% Calculate SNR using matlab functions
adc_SNR = snr(x,fs,3); % Already in dB
Enob_mat = (adc_SNR - 1.76)/6.02;
% Plot the signal to noise ratio
figure;
snr(x,fs, 3);
hold on;
tex3 = text(5, 20,"$ENOB =\mbox{}$" + Enob_mat + 'Bits','Interpreter','latex');
tex3.FontSize = 14;

f = linspace(0,0.5, N/2+1); % Normalized Frequency
fin = 5001;

%% Using the Delta Sigma Toolbox 
% create window:
%w = .5*(1 - cos(2*pi*(0:N-1)/N) ); % Hanning window
w = 0.54 - 0.46*cos(2*pi*(0:N-1)/N); % Hamming window
SPEC = fft(x.*w)/(N/4);

vfft = 20*log10(abs(SPEC(1:N/2+1))); % Only want 0 to fs/2
figure;
plot(f*fs, vfft);

grid on;
xlabel("frequency(Hz)")
ylabel("SNR(dB)")
title("Signal To Noise Ratio of ADC")

f1_bin = 9; % removing the DC component 
f2_bin = floor(N/2+1);
fin = floor((5001/fs)*N);
SPECsnr = SPEC(f1_bin:f2_bin);
snrCalc = calculateSNR(SPEC(f1_bin:f2_bin), fin , 9);

%% Calculate SNR
signalBins = [fin - 8: fin + 8];
signalBins = signalBins(signalBins > 0);
signalBins = signalBins(signalBins <= length(SPECsnr));
s = norm(SPECsnr(signalBins));
noiseBins = 1:length(SPECsnr);
noiseBins(signalBins) = [];
n = norm(SPECsnr(noiseBins));
snr = 20*log10(abs(s)/abs(n));


ENOB = (snrCalc - 1.76)/6.02;
tex1 = text(5000, 40,"$ENOB =\mbox{}$" + ENOB + 'Bits','Interpreter','latex');
tex2 = text(5000, 30,"$SNR =\mbox{}$" + snrCalc + 'dB','Interpreter','latex');
tex1.FontSize = 14;
tex2.FontSize = 14;

delta = Vref/(2^Enob_mat);

%% Plot the output voltage data
figure;
plot(time,x*delta);
title("Output of ADC")
xlabel("time(s)")
ylabel("Voltage(V)")
grid on
axis([0 2/5001 0 3.6]) ;