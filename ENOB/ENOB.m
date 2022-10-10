%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for calculating Effective Number of Bit for an ADC
% Uses FFT test algorithm and calulates the signal to noise ratio
%
%   Dexter Elmendorf 10/02/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

% Load data, calculate Fs
file = readmatrix("ADC_Out.dat");
voltage = transpose(file);
Ts = mean(diff(voltage(1,1:end)));
time = voltage(1,1:end);
voltage = voltage(2,1:end);

Fs = 1/Ts;

% Plot the output voltage data
figure;
plot(time,voltage);
title("Output of ADC")
xlabel("time(s)")
ylabel("Voltage(V)")
grid on

% Calculate SNR
adc_SNR = snr(voltage,Fs); % Already in dB

enob = (adc_SNR - 1.76)/6.02;

% Plot the signal to noise ratio
figure;
snr(voltage,Fs);