clear all;
close all;

%% Raw data
accelsID = fopen('accels2.txt','r');
accels = fscanf(accelsID, '%f');
%accels = accels(4:4:end);

figure(1);
subplot(2,1,1);
plot(accels);
title('Raw data: accels');

fs = 500;                       % sample frequency

n = pow2(nextpow2(length(accels)));    % Transform length
y = fft(accels,n);              % DFT
mag(1) = abs(y(1))/n;           % DC component
mag(n/2+1) = abs(y(n/2+1))/n;   % Nyquist frequency component
mag(2:n/2) = 2*abs(y(2:n/2))/n; % all other frequency components
freqs = linspace(0, 1, n/2+1);  % make x-axis as fraction of Nyquist freq

subplot(2,1,2);
plot(freqs, mag);               % plot the FFT magnitude plot
                             
% only need to plot from 0 to Nyquist frequency fs/2 (half the frequency 
% range), since second half is reflection of first.

xlabel('Frequency (Hz) as fraction of Nyquist Freq');
ylabel('Magnitude');
title('FFT of unfiltered accels');

%% MAF w/ 5 element FIFO
P = 5; 
fifo = [0 0 0 0 0]; 
accelsMAF = [];
for i=1:length(accels)
    fifo = [fifo(2:end) accels(i)];
    accelsMAF(i) = 1/(P)*sum(fifo); 
end

figure(2);
subplot(2,1,1);
plot(accelsMAF);
title('MAF Filtered accels');

%fft of MAF filtered
nMAF = pow2(nextpow2(length(accelsMAF)));
yMAF = fft(accelsMAF, nMAF);
magMAF(1) = abs(yMAF(1))/nMAF;                 % DC component
magMAF(nMAF/2+1) = abs(yMAF(nMAF/2+1))/nMAF;   % Nyquist frequency component
magMAF(2:nMAF/2) = 2*abs(yMAF(2:nMAF/2))/nMAF; % all other frequency components
freqsMAF = linspace(0, 1, nMAF/2+1);  % make x-axis as fraction of Nyquist freq

subplot(2,1,2);
plot(freqsMAF, magMAF); % plot the FFT magnitude plot

xlabel('Frequency (Hz) as fraction of Nyquist Freq');
ylabel('Magnitude');
title('FFT of MAF filtered accels');

%% FIR filter
b = fir1(10, 0.1); % 10th-order, 11-sample LPF with cutoff freq of 0.1 fN
figure(3);
freqz(b);

accelsFIR = conv(b, accels); 
figure(4);
subplot(2,1,1);
plot(accelsFIR);
title('FIR Filtered accels');

nFIR = pow2(nextpow2(length(accelsFIR)));
yFIR = fft(accelsFIR, nFIR);
magFIR(1) = abs(yFIR(1))/nFIR;           % DC component
magFIR(nFIR/2+1) = abs(yFIR(nFIR/2+1))/nFIR;   % Nyquist frequency component
magFIR(2:nFIR/2) = 2*abs(yFIR(2:nFIR/2))/nFIR; % all other frequency components
freqsFIR = linspace(0, 1, nFIR/2+1);  % make x-axis as fraction of Nyquist freq

subplot(2,1,2);
plot(freqsFIR, magFIR); % plot the FFT magnitude plot

xlabel('Frequency (Hz) as fraction of Nyquist Freq');
ylabel('Magnitude');
title('FFT of FIR filtered accels');