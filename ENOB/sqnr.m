#!/usr/bin/octave -qf

%
%script to produce a DFT plot and output key ratios
%   1000 time samples are stored in array x
%
%   copyright 2001  Eric Swanson
%
% Modified by gle on 7 Nov 2022
%

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$  INPUT REQUIRED FROM USER   $$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

% Debug mode or not

debug  = 0 ;

% Read from a file or not (1 = read the file)

readfile = 1 ;
filename = 'SARquantizer.dat' ;

% Resolution of ADC

Nbits = 10 ;

% Sampling frequency in Hz

fs = 100e3;

% Number of samples (always even)

Ns = 300 ;

% Input frequency

fin =  5001 ;

% Harmonics you wish to include in signal to distortion ratio calculation

harmonics = [2:3] ;

% Our true band of interest

fl = 3000 ;
fh = 13000 ;

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$    END USER INPUT SECTION   $$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

% Normalized frequency

fhat = fin / fs ;

% FFT Bin width in Hz

bin_wid = fs / Ns ;

% Find the bin where signal lives

signal_bin = floor( fin / bin_wid ) ;

% Fullscale output in bits

full_scale = 2 ** Nbits ;

% Amplitude of input signal

 A = full_scale / 2  ;

% *******************************************************
% Compute Hode window
% window and index defined in hodie_gle
% ********************************************************

N = Ns  ;
window = zeros(N, 1) ;
index = [0 : N-1]' ;
v = (2 * pi ) / N ;

%
%Hodie window cosine coefficients
%

a0 = 0.61640321314050;
a1 = 0.98537119272586;
a2 = 0.49603771622007;
a3 = 0.14992232793243;
a4 = 0.02458719103474;
a5 = 0.00176604651487;
a6 = 0.00003158118857;

%
%coefficients sum to N
%

for m=1:N;
    n1=m-.5;
    window(m,1)=a0-a1*cos(v*n1)+a2*cos(v*2*n1)-a3*cos(v*3*n1)+a4*cos(v*4*n1)-a5*cos(v*5*n1)+a6*cos(v*6*n1);
end;

% Plot window in time and frequency domains when in debug mode

if (debug == 1)
    figure() ;
    plot(index,window) ;
    axis([0 Ns 0 2.5]) ;
    title("Hodie Window Time Domain Plot", "fontsize", 16) ;
    xlabel("index value", "fontsize", 16) ;
    ylabel("window value", "fontsize", 16) ;
    grid on ;

    y = fft(window) ;
    ym = abs(y) ;
    ym = 20 * log10(ym) ;

    figure() ;
    plot(index, ym) ;
    axis([0 10 -300 100]) ;
    title("Hodie Window Frequency Domain Plot", "fontsize", 16) ;
    xlabel("index value", "fontsize", 16) ;
    ylabel("Magnitude (dB)", "fontsize", 16) ;
    grid on ;
endif

%
% If the readfile variable is set then we want to grab the x input vector from a the quantizer.dat file
%

if (readfile == 1)
    data = dlmread(filename, '\t' , 0, 0) ;
    x = data(1 : Ns ) ;
else
    vin = (full_scale / 2 ) + A * sin(2 * pi * fhat  .* index) ;
    x = round(vin) ;
endif

len_x = length(x) ;

%
% We need to apply Hode window to our data
%

x1 = window .* x ;

%
% Now we take the FFT
%
y = fft(x1) ;
%
% We are only interested out to fs/2

y(Ns / 2 + 2 : Ns) = [] ;

%
% Normalize (puts DC at 0 dB)
%

y1 = (abs(y) + 1e-25)  / ( Ns * (A/2)) ;
ydbfs = 20 * log10(y1) ;

%
% Compute ratios for  input
%

y2= y1 .* y1 ;

%
% Compute power in the "signal" bins
%

signal_bins = [signal_bin - 8 : signal_bin + 8] ;
len_signal_bins = length(signal_bins) ;
ysig = sum(y2(signal_bins)) ;

%
% Compute power at dc (first 9 bins)
%

dc_bins = [1 : 9] ;

len_dc_bins = length(dc_bins) ;
ydc = sum(y2(dc_bins)) ;

%
% Compute power in the distortion bins
%

harmonic_bins = [] ;
for i = 1 :  length(harmonics)
     sub_vector = harmonics(i) * signal_bin - 8  : harmonics(i) * signal_bin  + 8  ;
     harmonic_bins = [harmonic_bins  sub_vector] ;
end

ydist = sum(y2(harmonic_bins)) ;

%
% Compute power in noise bins
%
all_bins = [1: length(y2)] ;
len_all_bins = length(y2) ;

not_noise_bins = [dc_bins signal_bins harmonic_bins] ;

noise_bins = [] ;
for i = 1 : len_all_bins
    tmp = any(not_noise_bins == i) ;
    if (tmp == 0)
         noise_bins = [noise_bins  i]  ;
    endif
end


ynoise =sum(y2(noise_bins)) ;

%
%ratios
% Fudge factors account for hidden noise
%

len_noise_bins = length(noise_bins) ;
hidden_noise_factor_1 = len_all_bins / len_noise_bins ;
sn = 10 * log10(ysig / ( hidden_noise_factor_1 * ynoise)) ;
sd = 10 * log10(ysig / ydist) ;
hidden_noise_factor_2 = len_all_bins  / (len_all_bins - len_dc_bins) ;
snd = 10 * log10(ysig/(ydist + hidden_noise_factor_2 * ynoise)) ;

%
% Compute Effective Number of Bits (ENOB)
%

enob = (snd - 1.76) / 6.02 ;

%
% Compute ratios in our band of interest
%

fl_bin = floor(fl / bin_wid ) ;
fh_bin = floor(fh / bin_wid)  ;

noise_bins = [] ;

for i = fl_bin : fh_bin
    tmp = any(signal_bins == i) ;
    if (tmp == 0)
         noise_bins = [noise_bins  i]  ;
    endif
end

len_noise_bins = length(noise_bins) ;
hidden_noise_factor = 1 + (len_signal_bins / len_noise_bins) ;
ynoise =sum(y2(noise_bins)) ;
sqnr_dB = 10 * log10(ysig / (hidden_noise_factor * ynoise)) ;
enob_nb = (sqnr_dB - 1.76) / 6.02 ;

%
% ADC signal / quantization noise plot
%

freq = bin_wid * [0 : Ns / 2] ;

warning("off") ;

h = figure('name','SQNR Plot','numbertitle','off') ;
% figure() ;
plot(freq, ydbfs) ;
xmax = bin_wid * (Ns / 2) ;
ymax = -165 ;
axis([0 xmax  ymax 0]) ;
set (gca, "xminorgrid", "on");
set (gca, "yminorgrid", "on");
title("ADC Signal to Quantization Noise Plot", "fontsize", 16) ;
xlabel("Frequency (Hz)", "fontsize", 16) ;
ylabel("Magnitude (dB)", "fontsize", 16) ;

%
% Annotate the results on the plot!!!
%

str = sprintf("sn is %.2f dB\nsd is %.2f  dB\nsnd is %.2f\nENOB is %.2f bits\nNarrowband ENOB is %.2f bits", sn, sd, snd, enob, enob_nb) ;
xloc = fs / 6 ;
yloc = -25 ;
text(xloc, yloc, str, "fontsize", 16) ;

print(h,'-dpdf','-color', "enob.pdf") ;

grid on ;
hold off ;

