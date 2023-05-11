% For better understanding of the code and to see the output figure, please read the report
% Section 1 - Plot the original signal and remove the muscles noise
% load the signal
EKG1 = [struct2cell(load('ecg.mat')){:}];
fs = 500; % sampling frequency
time = (0:length(EKG1)-1)/fs; % time vector
figure
plot(time, EKG1);
set(gca, "fontsize", 24) % increase the fontsize of the plot
xlabel("Time (s)");
ylabel("Voltage (V)");
%%%%%%%%%%
% zoom in on the signal
t_start = 0.7; % start time for zoomed in plot
t_end = 1.5; % end time for zoomed in plot
idx_start = round(t_start*500); % index corresponding to start time
idx_end = round(t_end*500); % index corresponding to end time
figure
plot(time(idx_start:idx_end), EKG1(idx_start:idx_end), "linewidth", 2);
set(gca, "linewidth", 2, "fontsize", 24);
xlabel("Time (s)");
ylabel("Voltage (V)");
%%%%%%%%%%
% remove the muscles noise
freq = (0:length(EKG1)-1)/length(EKG1)*500; % frequency vector
ecg_fft = fft(EKG1); % Fourier transform of signal
ecg_fft(freq < 0.5) = 0; % set frequencies below 0.5 Hz to zero
ecg1 = real(ifft(ecg_fft)); % inverse Fourier transform of filtered signal
figure
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
figure
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
ecg2 = filter(b, a, ecg1); % apply notch filter using filtfilt
figure
plot(time, ecg1, 'b', time, ecg2, 'r');
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
plot(time, ecg2, 'b', "linewidth", 2); % plot original signal
set(gca, "linewidth", 2, "fontsize", 24);
xlabel("Time (s)");
ylabel("Voltage (V)");
legend("Original signal");

% plot filtered signals
for i = 1:length(fc)
    [b,a] = butter(order, fc(i)/(fs/2), 'low');
    ecg_filt = filter(b, a, ecg2);
    plot(time, ecg_filt, colors(i));
    set(gca, "linewidth", 2, "fontsize", 24);
end

title("Filtered Signals with Different Cut-off Frequencies");
legend("Original signal", "5", "10", "100");
hold off;
%%%%%%%%
%%%%%%%%
% Section 4 - Finding the heart rate using autocorrelation:
% The autocorrelation sequence of a periodic signal has the same cyclic characteristics as the signal itself.
% Thus, autocorrelation can help verify the presence of cycles and determine their durations.
[b, a] = butter(5, 40/(fs/2), 'low');
ecg3 = filter(b, a, ecg2);
% plot the filtered signal
figure
plot(time, ecg3, 'b');
set(gca, "linewidth", 2, "fontsize", 24);
xlabel("Time (s)");
ylabel("Voltage (V)");
legend("Filtered signal");
% get the autocorrelation of the filtered signal
[autocor,lags] = xcorr(ecg3,'coeff');
% plot the autocorrelation
figure
plot(lags/fs,autocor, "linewidth", 2)
xlabel('Lag')
ylabel('Autocorrelation')
% find the peaks of the autocorrelation
[peaks,locs] = findpeaks(autocor,'MinPeakheight',0.4, "DoubleSided");

% find the average time between peaks
period = mean(diff(locs)) / fs;
% find the heart rate by dividing 60 seconds by the average period
heart_rate = 60 / period
% plot the peaks on top of the autocorrelation
hold on
pks = plot(lags(locs)/fs,peaks,'or', "linewidth", 2);
set(gca, "linewidth", 2, "fontsize", 24);
title(sprintf('Heart rate: %.1f bpm', heart_rate));
hold off
legend(pks, 'peaks')
%%%%%%%%%%%%%%%%%%%%%%%%
% section 5 - Find QRS complex
function [qrs_peaks, ecg_filtered] = detect_qrs_complex(ecg_signal, fs, f1, f2, order, threshold_factor)
% Detect QRS complex in an ECG signal using the Pan-Tompkins algorithm
% Inputs:
%   ecg_signal: the ECG signal to process
%   fs: the sampling frequency of the ECG signal
%   f1: the lower cutoff frequency for the bandpass filter
%   f2: the higher cutoff frequency for the bandpass filter
%   order: the order of the Butterworth filter
%   threshold_factor: the factor to multiply the max amplitude of the filtered signal to set the threshold for QRS detection
% Output:
%   qrs_peaks: the indices of the R-peaks in the ECG signal

% Apply bandpass filter to the ECG signal
[b, a] = butter(order, [f1 f2]/(fs/2), 'bandpass');
ecg_filtered = filtfilt(b, a, ecg_signal);

% Set the threshold for QRS detection
threshold = threshold_factor * max(ecg_filtered);

% Find the R-peaks using a simple thresholding method
qrs_peaks = find(ecg_filtered > threshold);

end

f1 = 5;
f2 = 35;
order = 2;
threshold_factor = 0.5;

% Apply QRS detection to ecg2
[qrs_peaks_ecg2, ecg2_filtered]= detect_qrs_complex(ecg2, fs, f1, f2, order, threshold_factor);

figure
plot(time, ecg2, 'b', time, ecg2_filtered, 'r', 'linewidth', 2);
set(gca, "linewidth", 2, "fontsize", 24);
xlabel("Time (s)");
ylabel("Voltage (V)");
legend("Original signal", "Filtered signal");

% Apply QRS detection to EKG1
qrs_peaks_EKG1 = detect_qrs_complex(EKG1, fs, f1, f2, order, threshold_factor);

% Plot the original signal with the detected QRS complex for ecg2
figure;
subplot(2, 1, 1);
plot(time, ecg2, 'b', 'linewidth', 2);
set(gca, 'linewidth', 2, 'fontsize', 24);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Original Signal for ecg2');

subplot(2, 1, 2);
plot(time, ecg2, 'b', 'linewidth', 1);
hold on;
plot(time(qrs_peaks_ecg2), ecg2(qrs_peaks_ecg2), 'ro', 'linewidth', 2);
hold off;
set(gca, 'linewidth', 2, 'fontsize', 24);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('QRS Complex for ecg2');
legend('Original Signal', 'QRS Complex Peaks');

% Plot the original signal with the detected QRS complex for EKG1
figure;
subplot(2, 1, 1);
plot(time, EKG1, 'b', 'linewidth', 2);
set(gca, 'linewidth', 2, 'fontsize', 24);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Original Signal for EKG1');

subplot(2, 1, 2);
plot(time, EKG1, 'b', 'linewidth', 1);
hold on;
plot(time(qrs_peaks_EKG1), EKG1(qrs_peaks_EKG1), 'ro', 'linewidth', 2);
hold off;
set(gca, 'linewidth', 2, 'fontsize', 24);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('QRS Complex for EKG1');
legend('Original Signal', 'QRS Complex Peaks');

