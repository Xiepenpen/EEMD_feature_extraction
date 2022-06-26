%This file is to do some testing witn new peak detction method
clc();
clear();
close all

% Files = dir(fullfile('C:\Users\Administrator\Desktop\毕业设计\数据\初始数据\dataverse_files\data\mat\*.mat'));
% LengthFiles = length(Files);
% MAE = zeros(1,LengthFiles);
% mae_AM_peaks = zeros(1,LengthFiles);
% mae_PM_peaks = zeros(1,LengthFiles);
% diff_peaks = zeros(1,LengthFiles);
% diff_peaks_ref = zeros(1,LengthFiles);
% for i=1:LengthFiles
%     name=Files(i).name;           %读取struct变量的格式
%     folder=Files(i).folder;
% s = load([folder,'\',name]);

%初始信号配置
fs = 300;
t = 0:1/fs:480;

%0329
s = load('C:\Users\Administrator\Desktop\毕业设计\数据\初始数据\dataverse_files\data\mat\0009_8min.mat');
ecg_y = s.signal.ecg.y.'; %ecg
ecg_y = mapminmax(ecg_y, -10, 10);
ppg_y = s.signal.pleth.y.'; %ppg
ppg_y = mapminmax(ppg_y, 0, 10);
%ecg_ref_y = s.reference.hr.ecg.y; 
co2_y = s.signal.co2.y.'; %co2 
%[m_peaks_co2,locs_co2] = findpeaks(co2_y,t,'minpeakdistance',4);
%ecg_ref_x = s.reference.hr.ecg.x; 
rs_ref_y = s.reference.rr.co2.y; %reference respiration y
rs_ref_x = s.reference.rr.co2.x; %reference respiration x
temp_ref_y = zeros(1,15);
count =0;
for k = 1:15
    count =  k*32;
    temp_cal = abs(rs_ref_x-count);
    [min_num,colum] = min(temp_cal);
    temp_ref_y(k) = rs_ref_y(colum); 
end


%rs_ref_x = s.reference.rr.co2.x; %reference respiration x
rs_ref_x = 32:32:480; %reference respiration x
rs_ref_y = temp_ref_y;
peaks = s.labels.ecg.peak.x;
peaks_ppg = s.labels.pleth.peak.x;

data = ecg_y;
mu = mean(data);
sig = std(data);
%通用
[m_peaks_ecg,locs_ecg] = findpeaks(data,t,'minpeakheight',mu+0.8*sig,'minpeakdistance',0.415,'MinPeakProminence',1);


%HRV & Ecg_AM
hrvy=diff(m_peaks_ecg)*(1/fs); %HRV
%E_AM = ecg_y(m_peaks)-mean(ecg_y);
E_AM = m_peaks_ecg;

if length(E_AM)>length(peaks)
    temp_add = zeros(1,abs(length(E_AM)-length(peaks)));
    temp_add = temp_add+peaks(length(peaks));
    peaks = [peaks(1:length(peaks)),temp_add];
elseif length(E_AM) < length(peaks)
    temp_add = zeros(1,abs(length(E_AM)-length(peaks)));
    temp_add = temp_add+mu;
    E_AM = [E_AM(1:length(E_AM)),temp_add];
end

mae_AM_peaks = sum(abs(E_AM - ecg_y(peaks)))/length(E_AM);

%end

% figure
% subplot(3,1,1)
% plot(t,ecg_y)
% hold on
% plot(locs_ecg,m_peaks_ecg,'o', 'MarkerSize', 3)
% hold on
% %plot(t(m_valley),ecg_y(m_valley),'o', 'MarkerSize', 3)
% title('Original Signal');
% subplot(3,1,2)
% %plot(ecg_y(peaks))
% plot(E_AM)
% hold on
% plot(ecg_y(peaks))
% title('E-AM');
% xlabel('Time(Second)') 
% ylabel('Amplitude') 
% 
% subplot(3,1,3)
% plot(hrvy)
% title('HRV');

%end
%PPG signal
data = ppg_y;
mu = mean(ppg_y);
sig = std(ppg_y);

[m_peaks_ppg,locs_ppg] = findpeaks(ppg_y,t,'minpeakheight',mu+0.5*sig,'minpeakdistance',0.26);

%PRV & PPG_AM
%m_peaks_ppg = m_peaks_ppg(1:2:length(m_peaks_ppg));
%locs_ppg = locs_ppg(1:2:length(locs_ppg));
P_AM = m_peaks_ppg;
prvy=diff(m_peaks_ppg)*(1/fs); %HRV

if length(P_AM)>length(peaks_ppg)
    temp_add = zeros(1,abs(length(P_AM)-length(peaks_ppg)));
    temp_add = temp_add+peaks_ppg(length(peaks_ppg));
    peaks_ppg = [peaks_ppg(1:length(peaks_ppg)),temp_add];
elseif length(P_AM) < length(peaks_ppg)
    temp_add = zeros(1,abs(length(P_AM)-length(peaks_ppg)));
    temp_add = temp_add+m_peaks_ppg(1:length(temp_add));
    P_AM = [P_AM(1:length(P_AM)),temp_add];
end

mae_PM_peaks = sum(abs(P_AM - ppg_y(peaks_ppg)))/length(P_AM);
%mae_PM_peaks = sum(abs(P_AM - ppg_y(peaks_ppg)))/length(P_AM);

% figure
% subplot(3,1,1)
% plot(t,ppg_y)
% hold on
% %plot(t(m_peaks_ppg),ppg_y(m_peaks_ppg),'o', 'MarkerSize', 3)
% plot(locs_ppg,m_peaks_ppg,'o', 'MarkerSize', 3)
% title('Original Signal');
% subplot(3,1,2)
% plot(P_AM)
% hold on
% plot(ppg_y(peaks_ppg))
% xlabel('Time(Second)') 
% ylabel('Amplitude') 
% title('P-AM');
% subplot(3,1,3)
% plot(prvy)
% title('PRV');
%%
% [m_peaks_co2,locs_co2] = findpeaks(co2_y,t,'minpeakdistance',2);
% [m_valley_co2,locs_co2_1] = findpeaks(-co2_y,t,'minpeakdistance',2);
% co2_ref = m_peaks_co2;

%% HRV(ECG) signal decomposition
%EEMD decomposion HRV signal
%get mix signal by IMF1 and IMF2
imf_HRV = pEEMDandFFT(hrvy,1,0.01,100);
imf2_HRV = imf_HRV(1,:);
imf3_HRV = imf_HRV(2,:);
imf4_HRV = imf_HRV(3,:);
imf_mix_HRV = imf3_HRV;%+imf4_HRV;
t_imf_mix = 0:1:length(imf_mix_HRV);

%t-aixs for 8min signal
complot_fs = 1/(480/length(imf_mix_HRV));
Ts_HRV = keep(1/keep(complot_fs,3),6);
t_imf_mix_1 = 0:Ts_HRV:480;
t_imf_mix_1 = t_imf_mix_1(1:length(imf_mix_HRV));


% %detect the peaks of imf_mix_HRV
imf_mix_HRV_zeromean = mapminmax(imf_mix_HRV, 0, 1);
[breath_num_HRV,save_peaks,locs_peaks,save_valleys,locs_valley,zero_point] = threePointDetect_can(imf_mix_HRV_zeromean,t_imf_mix_1,rs_ref_x);

figure
plot(t_imf_mix_1,imf_mix_HRV_zeromean)
hold on
plot(locs_peaks,save_peaks,'x')
hold on
plot(zero_point,mean(imf_mix_HRV_zeromean),'o')



%[temp_peaks,locs_peaks]= findpeaks(imf_mix_HRV_zeromean,t_imf_mix_1);
breath_num_HRV = comref(locs_peaks,t_imf_mix_1,rs_ref_x);

%HRV_MAE(i) = mae(breath_num_HRV,rs_ref_y);


%plot the figure
% figure
% subplot(4,1,1)
% plot(imf_mix_HRV)
% title('IMF2+IMF3')
% subplot(4,1,2)
% plot(t_imf_mix_1,imf_mix_HRV_zeromean)
% hold on
% %plot(locs_peaks,save_peaks,'o', 'MarkerSize', 3)
% plot(locs_peaks,temp_peaks,'o', 'MarkerSize', 3)
% title('IMF MIX Zeromean')
% subplot(4,1,3)
% plot(rs_ref_x,breath_num_HRV)
% subplot(4,1,4)
% plot(rs_ref_x,rs_ref_y)
%% PRV(PPG) signal decomposition

imf_PRV = pEEMDandFFT(prvy,1,0.01,100);
imf2_PRV = imf_PRV(1,:);
imf3_PRV = imf_PRV(2,:);
imf4_PRV = imf_PRV(3,:);
imf_mix_PRV = imf3_PRV;%+imf4_PRV;
t_imf_mix = 0:1:length(imf_mix_PRV);

%t-aixs for 8min signal
complot_fs = 1/(480/length(imf_mix_PRV));
Ts_PRV = keep(1/keep(complot_fs,3),6);
t_imf_mix_1 = 0:Ts_PRV:480;
t_imf_mix_1 = t_imf_mix_1(1:length(imf_mix_PRV));

imf_mix_PRV_zeromean = mapminmax(imf_mix_PRV, 0, 1);
%[breath_num_PRV,save_peaks,locs_peaks,mean_gap] = threePointDetect_can(imf_mix_PRV_zeromean,t_imf_mix_1,rs_ref_x);
[temp_peaks,locs_peaks]= findpeaks(imf_mix_PRV_zeromean,t_imf_mix_1);
breath_num_PRV = comref(locs_peaks,t_imf_mix_1,rs_ref_x);

PRV_MAE(i) = mae(breath_num_PRV,rs_ref_y);
% figure
% subplot(4,1,1)
% plot(imf_mix_PRV)
% title('IMF2+IMF3')
% subplot(4,1,2)
% plot(t_imf_mix_1,imf_mix_PRV_zeromean)
% hold on
% plot(locs_peaks,save_peaks,'o', 'MarkerSize', 3)
% title('IMF MIX Zeromean')
% subplot(4,1,3)
% plot(rs_ref_x,breath_num_PRV)
% subplot(4,1,4)
% plot(rs_ref_x,rs_ref_y)




%% ECG_AM signal decomposition
%EEMD decomposion AM signal
imf_AM = pEEMDandFFT(E_AM,1,0.01,100);


imf2_AM = imf_AM(1,:);
imf3_AM = imf_AM(2,:);
imf4_AM = imf_AM(3,:);
imf_mix_AM = imf3_AM;%+imf4_AM;
%imf_mix_AM = imf_mix_AM(1:814);
t_imf_mix = 0:1:length(imf_mix_AM);

complot_fs = 1/(480/length(imf_mix_AM));
Ts_AM = 1/keep(complot_fs,3);
t_imf_mix_1 = 0:Ts_AM:480;
t_imf_mix_1 = t_imf_mix_1(1:length(imf_mix_AM));


imf_mix_AM_zeromean_ecg = mapminmax(imf_mix_AM, 0, 1);
%imf_mix_HRV_zeromean(imf_mix_HRV_zeromean==0)=0.0001;
%[breath_num_AM,save_peaks,locs_peaks,mean_gap] = threePointDetect_can(imf_mix_AM_zeromean,t_imf_mix_1,rs_ref_x);
[temp_peaks,locs_peaks]= findpeaks(imf_mix_AM_zeromean_ecg,t_imf_mix_1);
breath_num_AM_ecg = comref(locs_peaks,t_imf_mix_1,rs_ref_x);

EAM_MAE(i) = mae(breath_num_AM_ecg,rs_ref_y);
% figure
% subplot(4,1,1)
% plot(imf_mix_AM)
% title('IMF2+IMF3')
% subplot(4,1,2)
% plot(t_imf_mix_1,imf_mix_AM_zeromean)
% hold on
% plot(locs_peaks,save_peaks,'o', 'MarkerSize', 3)
% %hold on
% %plot([1,length(imf_mix_AM_zeromean)],[mean(imf_mix_AM_zeromean),mean(imf_mix_AM_zeromean)])
% title('IMF MIX Zeromean')
% subplot(4,1,3)
% plot(rs_ref_x,breath_num_AM)
% subplot(4,1,4)
% plot(rs_ref_x,rs_ref_y)

%% PPG_AM signal decomposition

imf_AM_ppg = pEEMDandFFT(P_AM,1,0.01,100);


imf2_AM_ppg = imf_AM_ppg(1,:);
imf3_AM_ppg = imf_AM_ppg(2,:);
imf4_AM_ppg = imf_AM_ppg(3,:);
imf_mix_AM_ppg = imf3_AM_ppg;%+imf4_AM_ppg;
%imf_mix_AM = imf_mix_AM(1:814);
t_imf_mix = 0:1:length(imf_mix_AM_ppg);

complot_fs = 1/(480/length(imf_mix_AM_ppg));
Ts_AM = 1/keep(complot_fs,3);
t_imf_mix_1 = 0:Ts_AM:480;
t_imf_mix_1 = t_imf_mix_1(1:length(imf_mix_AM_ppg));

%detect the peaks of imf_mix_HRV
imf_mix_AM_zeromean_ppg = mapminmax(imf_mix_AM_ppg, 0, 1);
%imf_mix_HRV_zeromean(imf_mix_HRV_zeromean==0)=0.0001;
%[breath_num_AM_ppg,save_peaks_ppg,locs_peaks_ppg,mean_gap] = threePointDetect_can(imf_mix_AM_zeromean_ppg,t_imf_mix_1,rs_ref_x);
[temp_peaks,locs_peaks]= findpeaks(imf_mix_AM_zeromean_ppg,t_imf_mix_1);
breath_num_AM_ppg = comref(locs_peaks,t_imf_mix_1,rs_ref_x);


PAM_MAE(i) = mae(breath_num_AM_ppg,rs_ref_y);
% figure
% subplot(4,1,1)
% plot(imf_mix_AM_ppg)
% title('IMF2+IMF3')
% subplot(4,1,2)
% plot(t_imf_mix_1,imf_mix_AM_zeromean_ppg)
% hold on
% plot(locs_peaks_ppg,save_peaks_ppg,'o', 'MarkerSize', 3)
% %hold on
% %plot([1,length(imf_mix_AM_zeromean)],[mean(imf_mix_AM_zeromean),mean(imf_mix_AM_zeromean)])
% title('IMF MIX Zeromean')
% subplot(4,1,3)
% plot(rs_ref_x,breath_num_AM_ppg)
% subplot(4,1,4)
% plot(rs_ref_x,rs_ref_y)
%end
%% Evaluation of ECG
% 
% 
% SQI_AM = 0.25;
% SQI_HRV = 1-SQI_AM;
% 
% if length(imf_mix_AM_zeromean_ecg)>length(imf_mix_HRV_zeromean)
%       ECG=[(SQI_AM+0.2)*imf_mix_AM_zeromean_ecg(1:length(imf_mix_HRV_zeromean))+(SQI_HRV-0.2)*imf_mix_HRV_zeromean imf_mix_AM_zeromean_ecg(length(imf_mix_HRV_zeromean)+1:end)];
% elseif length(imf_mix_AM_zeromean_ecg)<length(imf_mix_HRV_zeromean)
%       ECG=[(SQI_HRV-0.2)*imf_mix_HRV_zeromean(1:length(imf_mix_AM_zeromean_ecg))+(SQI_AM+0.2)*imf_mix_AM_zeromean_ecg imf_mix_HRV_zeromean(length(imf_mix_AM_zeromean_ecg)+1:end)];
% else
%       ECG=(SQI_AM+0.2)*imf_mix_AM_zeromean_ecg+(SQI_HRV-0.2)*imf_mix_HRV_zeromean;
% end
% 
% ECG = mapminmax(ECG, 0, 1);
% 
% if length(imf_mix_AM_zeromean_ppg)>length(imf_mix_PRV_zeromean)
%       PPG=[SQI_AM*imf_mix_AM_zeromean_ppg(1:length(imf_mix_PRV_zeromean))+SQI_HRV*imf_mix_PRV_zeromean imf_mix_AM_zeromean_ppg(length(imf_mix_PRV_zeromean)+1:end)];
% elseif length(imf_mix_AM_zeromean_ppg)<length(imf_mix_PRV_zeromean)
%       PPG=[SQI_HRV*imf_mix_PRV_zeromean(1:length(imf_mix_AM_zeromean_ppg))+SQI_AM*imf_mix_AM_zeromean_ppg imf_mix_PRV_zeromean(length(imf_mix_AM_zeromean_ppg)+1:end)];
% else
%       PPG = SQI_AM*imf_mix_AM_zeromean_ppg+SQI_HRV*imf_mix_PRV_zeromean;
% end
% 
% PPG = mapminmax(PPG, 0, 1);
% 
% complot_fs = 1/(480/length(ECG));
% Ts_AM = 1/keep(complot_fs,3);
% t_imf_mix_1 = 0:Ts_AM:480;
% t_imf_mix_1 = t_imf_mix_1(1:length(ECG));
% [peaks_ecg_b,locs_ecg_b]= findpeaks(ECG,t_imf_mix_1);
% bn_ecg = comref(locs_ecg_b,t_imf_mix_1,rs_ref_x);
% 
% complot_fs = 1/(480/length(PPG));
% Ts_AM = 1/keep(complot_fs,3);
% t_imf_mix_1 = 0:Ts_AM:480;
% t_imf_mix_1 = t_imf_mix_1(1:length(PPG));
% [peaks_ppg_b,locs_ppg_b]= findpeaks(PPG,t_imf_mix_1);
% bn_ppg = comref(locs_ppg_b,t_imf_mix_1,rs_ref_x);
% 
% 
% 
% % ECG_breath_num = (SQI_AM+0.2)*breath_num_AM_ecg+(SQI_HRV-0.2)*breath_num_HRV;
% % PPG_breath_num = SQI_AM*breath_num_AM_ppg+SQI_HRV*breath_num_PRV;
% % 
% % MIX = 0.5*ECG_breath_num+0.5*PPG_breath_num;
% % figure
% % subplot(4,1,1)
% % plot(rs_ref_x,breath_num_AM)
% % subplot(4,1,2)
% % plot(rs_ref_x,breath_num_HRV)
% % subplot(4,1,3)
% % plot(rs_ref_x,test)
% % subplot(4,1,4)
% % plot(rs_ref_x,rs_ref_y)
% if length(PPG) > length(ECG)
%     co2_y = PPG;
% elseif length(PPG) < length(ECG)
%     co2_y = ECG;
% else
%     co2_y = ECG;
% end
% savepath = ['C:\Users\Administrator\Desktop\毕业设计\数据\初始数据\EEDM提取数据\480s数据','\',name];
% save(savepath,'imf_mix_HRV_zeromean','imf_mix_PRV_zeromean','imf_mix_AM_zeromean_ppg','imf_mix_AM_zeromean_ecg','ECG','PPG','rs_ref_y','co2_y')
%savepath = ['C:\Users\Administrator\Desktop\毕业设计\数据\初始数据\EEDM提取数据\480_分钟呼吸率','\',name];
%save(savepath,'breath_num_HRV','breath_num_PRV','breath_num_AM_ecg','breath_num_AM_ppg','ECG_breath_num','PPG_breath_num','rs_ref_y')





% PRV 优于 PPG——AM  HRV优于ECG——AM ECG优于PPG
% MAE_ECG_AM(i) = sum(abs(breath_num_AM_ecg - rs_ref_y))/length(breath_num_AM_ecg);
% MAE_ECG_HRV(i) = sum(abs(breath_num_HRV - rs_ref_y))/length(breath_num_HRV);
% MAE_PPG_AM(i) = sum(abs(breath_num_AM_ppg - rs_ref_y))/length(breath_num_AM_ppg);
% MAE_PPG_PRV(i) = sum(abs(breath_num_PRV - rs_ref_y))/length(breath_num_PRV);
% MAE_ECG(i) = sum(abs(bn_ecg - rs_ref_y))/length(bn_ecg);
% MAE_PPG(i) = sum(abs(bn_ppg - rs_ref_y))/length(bn_ppg);
% MAE_all(i) = sum(abs(MIX - rs_ref_y))/length(MIX);
%MAE = sum(abs(test - rs_ref_y))/length(test);
%E = MAE/abs(mean(rs_ref_y));

%end
%%
% index = zeros(1,length(mae_AM_peaks));
% for i = 1:length(mae_AM_peaks)
%    if mae_AM_peaks(i) >= 1
%        index(i)=i;
%    end
% end
% index(index==0)=[];
% 
% m1 = sort(MAE_ECG);
% m2 = sort(MAE_PPG);
% % m3 = sort(MAE_all);
% % m4 = sort(MAE_ECG_AM);
% % m5 = sort(MAE_ECG_HRV);
% % m6 = sort(MAE_PPG_AM);
% % m7 = sort(MAE_PPG_PRV);
% m1=[m1(11),m1(21),m1(32)]
% m2=[m2(11),m2(21),m2(32)]
% m3=[m3(11),m3(21),m3(32)]
% m4=[m4(11),m4(21),m4(32)]
% m5=[m5(11),m5(21),m5(32)]
% m6=[m6(11),m6(21),m6(32)]
% m7=[m7(11),m7(21),m7(32)]
%%
PAM_MAE = sort(PAM_MAE);
EAM_MAE = sort(EAM_MAE);
HRV_MAE = sort(HRV_MAE);
PRV_MAE = sort(PRV_MAE);


[PAM_MAE(11),PAM_MAE(21),PAM_MAE(32),mean(PAM_MAE)]
[EAM_MAE(11),EAM_MAE(21),EAM_MAE(32),mean(EAM_MAE)]
[HRV_MAE(11),HRV_MAE(21),HRV_MAE(32),mean(HRV_MAE)]
[PRV_MAE(11),PRV_MAE(21),PRV_MAE(32),mean(PRV_MAE)]


