% ELEC2430 Communications Theory
% Matlab Assignment 1
% Amplitude Modulation System and Spectrum Analysis
% School of Electronic & Electrical Engineering 
% University of Leeds

% Last Update: 25 Mar 2018

clc;
clear;

% 1.Binary Input Data
N = 1000; % Number of bits, bps
data = rand(1,N)>0.5; % Binary Input Data
t = linspace(0,N*0.001,N*8); % Generate N*8 points between 0 and N*0.01

% 2.Sampling
fs = 8000; % Sampling frequency, samples/sec
k = 1;
for i = 1:N
    for j = 1:8            % Create the  % samplePerBit = fs/N = 8
        y(k) = data(i);    % bit stream
        k = k+1;
    end
end

% 3.DSB-SC Modulation (Upconverter)
fc = 2000; % Carrier frequency, Hz
c = cos(2*pi*fc*t); % Sampled signal/Carrier Signal
C = fft(c);
m = c.*y; % Modulated Signal
M = fft(m);

wc = 2*1000/fs;
[b,a] = butter(5,wc);
Drop = 10000;
snr = 1; % dB
for drop = 1:Drop
    for kk = 1:length(snr)
        y2 = awgn(m,snr(kk),'measured');
        Y2 = fft(y2);
        y3 = c.*y2;
        Y3 = fft(y3);
        Filtered_Sig = filtfilt(b,a,y3);
        FILTERED_SIG = fft(Filtered_Sig);
        for c3 = 1:N
            index = (c3-1)*8+1;
            y_decod(c3) = mean(Filtered_Sig(index:index+7))>0.25;
        end
        Erro = mod(data,double(y_decod));
        Number = sum(Erro);
        BER(kk) = Number/N;
    end
    Drop_BER(drop,:) = BER;
end
Average = mean(Drop_BER);