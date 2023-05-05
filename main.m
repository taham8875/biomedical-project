% For better understanding of the code and to see the output figure, please read the report

% Section 1 - Plot the original signal and remove the muscles noise
% load the signal
EKG1 = [struct2cell(load('ecg.mat')){:}];
fs = 500; % sampling frequency
time = (0:length(EKG1)-1)/fs; % time vector
plot(time, EKG1); 
set(gca, "fontsize", 24) % increase the fontsize of the plot
xlabel("Time (s)");
ylabel("Voltage (V)");
%%%%%%%%%%
% zoom in on the signal
t_start = 1.5; % start time for zoomed in plot
t_end = 1.52; % end time for zoomed in plot
idx_start = round(t_start*500); % index corresponding to start time
idx_end = round(t_end*500); % index corresponding to end time
plot(time(idx_start:idx_end), EKG1(idx_start:idx_end), "linewidth", 2);
set(gca, "linewidth", 2, "fontsize", 24);
xlabel("Time (s)");
ylabel("Voltage (V)");
%%%%%%%%%%
% remove the muscles noise
freq = (0:length(EKG1)-1)/length(EKG1)*500; % frequency vector
ecg_fft = fft(EKG1); % Fourier transform of signal
ecg_fft(freq < 0.5) = 0; % set frequencies below 0.5 Hz to zero
ecg1 = ifft(ecg_fft); % inverse Fourier transform of filtered signal
plot(time, EKG1, 'b', time, ecg1, 'r');
set(gca, "linewidth", 2, "fontsize", 24);
xlabel("Time (s)");
ylabel("Voltage (V)");
legend("Original signal", "Filtered signal");
%%%%%%%%
% Make sure that you filter has a real impulse response
filter_fft = ones(size(EKG1)); % initialize filter to all ones
filter_fft(freq < 0.5) = 0; % set frequencies below 0.5 Hz to zero
filter_ifft = ifft(ifftshift(filter_fft)); % inverse Fourier transform of filter
t = (-length(filter_ifft)/2:length(filter_ifft)/2-1)/fs; % create time vector
plot(t, real(fftshift(filter_ifft)), "linewidth", 2); % center filter in time domain before plotting
set(gca, "linewidth", 2, "fontsize", 24);
xlabel("Time (s)");
ylabel("Amplitude");
%%%%%%%%
%%%%%%%%
% Section 2 - remove the muscles noise with a notch filter
pkg load signal % I am using octave ¯\_(ツ)_/¯
fs = 500; % sampling frequency
f0 = 50; % notch frequency
bw = 1/(fs/2); % normalized notch bandwidth
Q = f0/bw; % quality factor
wo = f0/(fs/2); % normalized notch frequency
[b,a] = butter(2, [wo-bw/2, wo+bw/2], 'stop'); % design notch filter coefficients
ECG_signal_filt = filter(b, a, EKG1); % apply notch filter using filtfilt
plot(time, EKG1, 'b', time, ECG_signal_filt, 'r');
set(gca, "linewidth", 2, "fontsize", 24);
xlabel("Time (s)");
ylabel("Voltage (V)");
legend("Original signal", "Filtered signal");
%%%%%%%%
%%%%%%%%
% Section 3 - Increasing the signal-to-noise ratio (Low pass filter):
% Try different cut-off frequencies and see how the signal changes 
fs = 500; % sampling frequency in Hz
order = 2; % filter order
fc = [5 20 100]; % cut-off frequencies in Hz
colors = ['r', 'g', 'y'];
figure;
hold on;
plot(time, ECG_signal_filt, 'b', "linewidth", 2); % plot original signal
set(gca, "linewidth", 2, "fontsize", 24);
xlabel("Time (s)");
ylabel("Voltage (V)");
legend("Original signal");

% plot filtered signals
for i = 1:length(fc)
    [b,a] = butter(order, fc(i)/(fs/2), 'low');
    ecg_filt = filter(b, a, ECG_signal_filt);
    plot(time, ecg_filt, colors(i));
    set(gca, "linewidth", 2, "fontsize", 24);
end

title("Filtered Signals with Different Cut-off Frequencies");
legend("Original signal", "5", "10", "100");
hold off;
%%%%%%%%
%%%%%%%%
% Section 4 - Finding the heart rate using autocorrelation:
% As seen from the output figures in section 3, the lower the cut-off frequency, the better the signal looks.
% Therefore, we will use the signal filtered with a cut-off frequency of 3 Hz to find the heart rate.  
% The autocorrelation sequence of a periodic signal has the same cyclic characteristics as the signal itself.
% Thus, autocorrelation can help verify the presence of cycles and determine their durations.         
[b, a] = butter(5, 3/(fs/2), 'low');
ecg3 = filter(b, a, ECG_signal_filt);
% plot the filtered signal
plot(time, ecg3, 'b');
set(gca, "linewidth", 2, "fontsize", 24);
xlabel("Time (s)");
ylabel("Voltage (V)");
legend("Filtered signal");
% get the autocorrelation of the filtered signal
[autocor,lags] = xcorr(ecg3,'coeff');
% plot the autocorrelation
plot(lags/fs,autocor, "linewidth", 2)
xlabel('Lag')
ylabel('Autocorrelation')
% find the peaks of the autocorrelation
[peaks,locs] = findpeaks(autocor, "DoubleSided");
% find the average time between peaks
period = mean(diff(locs)) / fs
% find the heart rate by dividing 60 seconds by the average period
heart_rate = 60 / period
% plot the peaks on top of the autocorrelation
hold on
pks = plot(lags(locs)/fs,peaks,'or', "linewidth", 2);
set(gca, "linewidth", 2, "fontsize", 24);
title(sprintf('Heart rate: %.1f bpm', heart_rate));
hold off
legend(pks, 'peaks')