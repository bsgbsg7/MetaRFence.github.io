clc;
clear all;
close all;


ref_df = readtable('ref1.csv', 'HeaderLines', 2, 'VariableNamingRule', 'preserve'); 
ref_grouped = findgroups(ref_df.EPC);

% Initialize arrays to store mean phase and amplit ude for each tag
mean_phase = zeros(max(ref_grouped), 1);
mean_amplitude = zeros(max(ref_grouped), 1);

for i = 1:max(ref_grouped)
    current_group_data = ref_df(ref_grouped == i, :);
    mean_phase(i) = mean(current_group_data.PhaseAngle);
    mean_amplitude(i) = mean(current_group_data.RSSI);
end

static_df = readtable('1.csv', 'HeaderLines', 2, 'VariableNamingRule', 'preserve'); 
moving_df = readtable('2.csv', 'HeaderLines', 2, 'VariableNamingRule', 'preserve'); 


static_grouped = findgroups(static_df.EPC);
moving_grouped = findgroups(moving_df.EPC);

epcs_static = unique(static_df.EPC);  
epcs_moving = unique(moving_df.EPC); 

static_time =  datetime(static_df.("// Timestamp"), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSX', 'TimeZone', 'local');
moving_time =  datetime(moving_df.("// Timestamp"), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSX', 'TimeZone', 'local');



time_intervals_static = diff(static_time);
desired_sampling_interval = seconds(1/35);

complex_signals_static = [];
for k = 1:length(epcs_static)
    idx_static = static_grouped == k;  
    rssi_phase_data_static = static_df(idx_static, {'RSSI', 'PhaseAngle'}); 
    rssi_phase_data_static = rssi_phase_data_static{1:end, :};
    b = sqrt(10.^(((rssi_phase_data_static(:,1)-mean_amplitude(k)) / 10)-3)) .* exp(1j * (rssi_phase_data_static(:,2)-mean_phase(k)));  
    b = angle(b);
    time_static_tag = static_time(idx_static);  
    time_interpolated = static_time(1):desired_sampling_interval:static_time(end); 
    signal_interpolated = interp1(time_static_tag, b.', time_interpolated, 'linear', 'extrap').';  
    complex_signals_static = [complex_signals_static,signal_interpolated];
end

complex_signals_moving = [];
for k = 1:length(epcs_moving)
    idx_moving = moving_grouped == k;
    rssi_phase_data_moving = moving_df(idx_moving, {'RSSI', 'PhaseAngle'});
    rssi_phase_data_moving = rssi_phase_data_moving{1:end, :};
    b = sqrt(10.^(((rssi_phase_data_moving(:,1)-mean_amplitude(k)) /10)-3)) .* exp(1j * (rssi_phase_data_moving(:,2)-mean_phase(k)));
    b = angle(b);
    time_moving_tag = moving_time(idx_moving);  
    time_interpolated = moving_time(1):desired_sampling_interval:moving_time(end); 
    signal_interpolated = interp1(time_moving_tag, b.', time_interpolated, 'linear', 'extrap').';  
    complex_signals_moving = [complex_signals_moving,signal_interpolated];
end


min_index = min(size(complex_signals_static,1),size(complex_signals_moving,1));
complex_signals_static = complex_signals_static(1:min_index,:);
complex_signals_moving = complex_signals_moving(1:min_index,:);

Y1_list = [];
Y2_list = [];

plot_size = 1;
figure;

for k = 1:length(epcs_moving)
    index = 1;
    for num = 1:1800*plot_size:length(complex_signals_moving)-1800*plot_size

        mean_phase_moving = complex_signals_moving(num:num+1800*plot_size,:);
        mean_phase_moving = mean_filter(mean_phase_moving(:,k),1800);

        mean_phase_static = complex_signals_static(num:num+1800*plot_size,:);
        mean_phase_static = mean_filter(mean_phase_static(:,k),1800);
       
%         figure;
%         subplot(2,1,1);
%         histogram(mean_phase_moving, 50, 'FaceColor', 'g');
%         xlabel('mean_phase_moving');
%         ylabel('频数');
%         subplot(2,1,2);
%         histogram(mean_phase_static, 50, 'FaceColor', 'g');
%         xlabel('mean_phase_static');
%         ylabel('频数');

        %mean_phase_static = mean_phase_static(:,k);
    
%         [moving_IMF, moving_u_hat,moving_CenFs_mvmd] = MVMD(mean_phase_moving, 2000, 0,6,0, 0, 1e-7); 
%         [static_IMF, static_u_hat,static_CenFs_mvmd] = MVMD(mean_phase_static, 2000, 0,6,0, 0, 1e-7); 
%         figure;
%         name='MVMD分解动态分量';[MVMD_IMF3]=Huatu(name,static_IMF,mean_phase_static,35);

    
        Fs = 35; 

        mean_phase_static = detrend(bandpass(mean_phase_static,[0.1,0.7],Fs));
        mean_phase_moving = detrend(bandpass(mean_phase_moving,[0.1,0.7],Fs));
        fft_num = 1;
        
        L1 = length(mean_phase_moving); 
        Y1 = fft((mean_phase_moving), L1 * fft_num); 
        frequencies_moving = (0:Fs/(fft_num*L1):Fs/2)  ;
        Y1_list(:, end+1) = abs(Y1(1:L1 * fft_num / 2 + 1));
        
        L2 = length(mean_phase_static); 
        Y2 = fft((mean_phase_static), L2 * fft_num); 
        frequencies_static = (0:Fs/(fft_num*L2):Fs/2)  ; 
        Y2_list(:, end+1) = abs(Y2(1:L2 * fft_num / 2 + 1));
        
        subplot(6, 3, 1 * (k - 1) + floor((k - 1) / 3) * 6 + index);
        plot(frequencies_moving*60, (abs(Y1(1:L1 * fft_num / 2 + 1)) * 2) / (fft_num * L1), 'red', 'DisplayName', 'mean\_phase\_moving');
%         xlim([0,0.7]);
        hold on;
        plot(frequencies_static*60, (abs(Y2(1:L2 * fft_num / 2 + 1)) * 2) / (fft_num * L2), '-black', 'DisplayName', 'mean\_phase\_static');
        xlabel('Breaths per minute');
%         xlim([0,0.7]);
        ylabel('Spectrum amplitude');
%         legend;
        grid on;
        
        subplot(6, 3, 1 * (k - 1) + floor((k - 1) / 3) * 6 + index + 3);
        plot(mean_phase_static);
        xlabel('Time');
        ylabel('Mean Phase Static');
        grid on;
        
        subplot(6, 3, 1 * (k - 1) + floor((k - 1) / 3) * 6 + index + 3 * 2);
        plot(mean_phase_moving);
        xlabel('Time');
        ylabel('Mean Phase Moving');
        grid on;

        index = index + 1;
    
        if index > 1
            break
        end
    
        % subplot(4,9,k+9*3);
        % plot(static_IMF(1,:), '-black', 'DisplayName', 'mean\_phase\_static');
        % hold on;
        % plot(moving_IMF(1,:), 'red', 'DisplayName', 'mean\_phase\_moving');
        % %legend;   
        % grid on;
    end
end

% figure;
% % subplot(2,1,1);
% plot(f2, mean(Y2_list,2), '-black', 'DisplayName', 'mean\_phase\_static');
% ylim([0,max(max(mean(Y2_list,2)),max(mean(Y1_list,2)))+0.2]);
% hold on;
% % legend;
% % subplot(2,1,2);
% plot(f1, mean(Y1_list,2), 'red', 'DisplayName', 'mean\_phase\_moving');
% legend;
% ylim([0,max(max(mean(Y2_list,2)),max(mean(Y1_list,2)))+0.2]);

% 
% figure;
% subplot(2, 1, 1);
% plot(f2, abs(Y2(1:L2*fft_num/2+1)).^2);
% title('Single-Sided Power Spectrum of mean\_phase\_static');
% xlabel('Frequency (Hz)');
% ylabel('Power');
% subplot(2, 1, 2);
% plot(f1, abs(Y1(1:L1*fft_num/2+1)).^2);
% title('Single-Sided Power Spectrum of mean\_phase\_moving');
% xlabel('Frequency (Hz)');
% ylabel('Power');
% 
% figure;
% subplot(2, 2, 1);
% plot(mean_magnitude_static);
% title('mean_magnitude_static');
% xlabel('Time');
% ylabel('Magnitude');
% ylim([0.025,0.035]);
% 
% subplot(2, 2, 2);
% plot(mean_phase_static);
% title('mean_phase_static');
% xlabel('Time');
% ylabel('Phase');
% ylim([-0.2,0.2]);
% 
% subplot(2, 2, 3);
% plot(mean_magnitude_moving);
% title('mean_magnitude_moving');
% xlabel('Time');
% ylabel('Magnitude');
% ylim([0.025,0.035]);
% 
% subplot(2, 2, 4);
% plot(mean_phase_moving);
% title('mean_phase_moving');
% xlabel('Time');
% ylabel('Phase');
% ylim([-0.2,0.2]);


% mean_phase_moving_detrended = detrend(mean_phase_moving);
% window_length = 64; 
% overlap_length = 32;
% sampling_frequency = 30; 
% nfft = 128; 
% [S, F, T] = spectrogram(mean_phase_moving_detrended, window_length, overlap_length, nfft, sampling_frequency);
% figure;
% imagesc(T, F, 10*log10(abs(S)));
% colorbar;
% axis xy;
% title('Short-time Fourier Transform of Mean Phase (Moving)');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');

