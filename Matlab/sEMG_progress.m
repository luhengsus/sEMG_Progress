function result = sEMG_progress(sEMG_raw,fs_EMG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is basic data progressing data of emg signals             %
% input:EMG_raw = n*n matrix of raw sEMG signals                 %
%       fs      = sampling rate                                  %
%       output  = 带通滤波,全波整流,均方根平滑,线性包络                %
% Author: Heng (2024-05-08)                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize processed signal matrices
num_samples   = size(sEMG_raw, 1);
num_channels  = size(sEMG_raw, 2);
sEMGs         = nan(num_samples,num_channels);
EMG_Bandfilt  = nan(num_samples,num_channels);
EMG_RMS       = nan(num_samples,num_channels);
EMG_FWR       = nan(num_samples,num_channels);
EMG_Lowerfilt = nan(num_samples,num_channels);
%% Step 1:remove Mean DC offset 去除直流电
for i          = 1:num_channels
    sEMG_cut   =  sEMG_raw(:,i);
    sEMG_nonan =  sEMG_cut(~isnan(sEMG_cut));
    sEMGs(1:length(sEMG_nonan),i)  =  sEMG_nonan - mean(sEMG_nonan);
end
%% Step 2:Bandpass filter带通滤波，频域分析用滤波后的数据！！
% filtering from 20 to 400 Hz
% Define filter parameters
Fnyq    = fs_EMG/2;   %Nuquist Frequemcy
highcut = 20/Fnyq;     %Highpass cutoff Frequency
lowcut  = 400/Fnyq;    %lowpass  cutoff Frequency
order1  = 4;    % Filter order
% Design the Butterworth filter
[b,a] = butter(order1,[highcut,lowcut],'bandpass');%4th order butterworth BP filter
% Apply the filter to the EMG data
for i = 1:num_channels
    EMG_Bandfilt(:,i)  = filtfilt(b,a,sEMGs(:,i));%Apply BP bandpass filter to dadta
end
%% Step 3:RMS envelop(moving average) -mean power of signal 均方根平滑
% 优先使用这个！！！
% window size in ms
% 算 RER 用50-75ms
% onset is 3std
window_size      = 0.1 * fs_EMG;%window size in ms取决于动作速度
for i            = 1:num_channels
    EMG_RMS(:,i) = sqrt(movmean((EMG_Bandfilt(:,i).^2),window_size));
end
%% (not use)Step 3: FFT analysis  Perform FFT on the filtered EMG data
% for i = 1:size(EMG_filt,2)
%     frequencies = (0:length(EMG_filt(i))-1) * (fs / length(EMG_filt(i)));
%     fft_emg = fft(EMG_filt);
%     fft_amplitude(:,i) = abs(fft_emg(1:length(EMG_filt(i))/2));
% end
%% Step 4.1:full wave rectification 全波整流
EMG_FWR = abs(EMG_Bandfilt);
%% Step 4.2:lower-pass filter 低通滤波Linear Envelope 线性包络
%  Butterworth low-pass filter
order2  = 4;  % Filter order (adjust as needed)
cof_fs  = 20; % Cutoff frequency for the low-pass filter (Hz)
Wn      = cof_fs/Fnyq;
[d,c]   = butter(order2 , Wn , 'low'); %lowpass 4 order 6Hz
for i   = 1:num_channels
    EMG_Lowerfilt(:,i) = filtfilt(d,c,EMG_FWR(:,i));%Apply BP lowpass filter to data
end
%% output data
result = [EMG_Bandfilt,EMG_FWR,EMG_RMS,EMG_Lowerfilt];%带通滤波,全波整流,均方根平滑,低通滤波
end