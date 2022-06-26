%This part is to split the signal from 480s to 15 sections each section is
%32s.
Files = dir(fullfile('C:\Users\Administrator\Desktop\毕业设计\EEMD_extraction\EEMDextraction\EEMD_feature_extraction\extracted_data\480s数据op\*.mat'));
LengthFiles = length(Files);
ECG_MAE = zeros(1,LengthFiles);
PPG_MAE = zeros(1,LengthFiles);
HRV_MAE = zeros(1,LengthFiles);
PRV_MAE = zeros(1,LengthFiles);
ECGAM_MAE = zeros(1,LengthFiles);
PPGAM_MAE = zeros(1,LengthFiles);
PPG_op_MAE = zeros(1,LengthFiles);
ECG_op_MAE = zeros(1,LengthFiles);

ECG_breathnum = zeros(1,15);
PPG_breathnum = zeros(1,15);
HRV_breathnum = zeros(1,15);
PRV_breathnum = zeros(1,15);
ECGAM_breathnum = zeros(1,15);
PPGAM_breathnum = zeros(1,15);
ECG_op_breathnum = zeros(1,15);
PPG_op_breathnum =zeros(1,15);
count = 1;


for i=1:LengthFiles
    name=Files(i).name;           
    folder=Files(i).folder;
    s = load([folder,'\',name]);
    
    co2_ref = s.co2_y;
    ecg_am = s.imf_mix_AM_zeromean_ecg;
    ppg_am = s.imf_mix_AM_zeromean_ppg;
    hrv = s.imf_mix_HRV_zeromean;
    prv = s.imf_mix_PRV_zeromean;
    reference = s.rs_ref_y;
    ecg = s.ECG;
    ppg = s.PPG;
    
    
    lenP1 = round(length(ecg_am)/15);
    lenP2 = round(length(ppg_am)/15);
    lenP3 = round(length(hrv)/15);
    lenP4 = round(length(prv)/15);
    lenP5 = round(length(ecg)/15);
    lenP6 = round(length(ppg)/15);
    
    %Use of sqi weighting within the fragment
    if length(ecg_am)>length(hrv)
        ecg_op=[SQI(1,i)*ecg_am(1:length(hrv))+SQI(2,i)*hrv ecg_am(length(hrv)+1:end)];
    elseif length(ecg_am)<length(hrv)
        ecg_op=[SQI(2,i)*hrv(1:length(ecg_am))+SQI(1,i)*ecg_am hrv(length(ecg_am)+1:end)];
    else
        ecg_op = SQI(1,i)*ecg_am+SQI(2,i)*hrv;
    end
   
    if length(ppg_am)>length(prv)
        ppg_op=[SQI1(1,i)*ppg_am(1:length(prv))+SQI1(2,i)*prv ppg_am(length(prv)+1:end)];
    elseif length(ppg_am)<length(prv)
        ppg_op=[SQI1(2,i)*prv(1:length(ppg_am))+SQI1(1,i)*ppg_am prv(length(ppg_am)+1:end)];
    else
        ppg_op = SQI1(1,i)*ppg_am+SQI1(2,i)*prv;
    end
    
    if length(ppg_op) > length(ecg_op)
        co2_ref = ppg_op;
    elseif length(ppg_op) < length(ecg_op)
        co2_ref = ecg_op;
    else
        co2_ref = ppg_op;
    end
    
    lenP7 = round(length(ecg_op)/15);
    lenP8 = round(length(ppg_op)/15);
    lenP = round(length(co2_ref)/15);
    
    %spliting the data 
    for j = 1:15
        
        if j~=15 & j ~=1
            co2 = co2_ref(lenP*(j-1):lenP*j);
            ECG_AM = ecg_am(lenP1*(j-1):lenP1*j);
            PPG_AM = ppg_am(lenP2*(j-1):lenP2*j);
            HRV = hrv(lenP3*(j-1):lenP3*j);
            PRV = prv(lenP4*(j-1):lenP4*j);
            ECG = ecg(lenP5*(j-1):lenP5*j);
            PPG = ppg(lenP6*(j-1):lenP6*j);
            
            ECG_op = ecg_op(lenP7*(j-1):lenP7*j);
            PPG_op = ppg_op(lenP8*(j-1):lenP8*j);
        elseif j == 15
            co2 = co2_ref(lenP*(j-1):length(co2_ref));
            ECG_AM = ecg_am(lenP1*(j-1):length(ecg_am));
            PPG_AM = ppg_am(lenP2*(j-1):length(ppg_am));
            HRV = hrv(lenP3*(j-1):length(hrv));
            PRV = prv(lenP4*(j-1):length(prv));
            ECG = ecg(lenP5*(j-1):length(ecg));
            PPG = ppg(lenP6*(j-1):length(ppg));
            
            ECG_op = ecg_op(lenP7*(j-1):length(ecg_op));
            PPG_op = ppg_op(lenP8*(j-1):length(ppg_op));
        elseif j ==1
            co2 = co2_ref(1:lenP*j);
            ECG_AM = ecg_am(1:lenP1*j);
            PPG_AM = ppg_am(1:lenP2*j);
            HRV = hrv(1:lenP3*j);
            PRV = prv(1:lenP4*j);
            ECG = ecg(1:lenP5*j);
            PPG = ppg(1:lenP6*j);
            
            ECG_op = ecg_op(1:lenP7*j);
            PPG_op = ppg_op(1:lenP8*j);
            
        end
        reference_value = reference(j);
        
        

           
        
        ECG_op_breathnum(j) = (length(findpeaks(ECG_op))/32)*60;
        PPG_op_breathnum(j) = (length(findpeaks(PPG_op))/32)*60;
        ECG_breathnum(j) = (length(findpeaks(ECG))/32)*60;
        PPG_breathnum(j) = (length(findpeaks(PPG))/32)*60;
        HRV_breathnum(j) = (length(findpeaks(HRV))/32)*60;
        PRV_breathnum(j) = (length(findpeaks(PRV))/32)*60;
        ECGAM_breathnum(j) = (length(findpeaks(ECG_AM))/32)*60;
        PPGAM_breathnum(j) = (length(findpeaks(PPG_AM))/32)*60;
        
%         
%         q = find('.'==name);
%         imname = name(1: q-1);
%         mkdir(['C:\Users\Administrator\Desktop\毕业设计\数据\初始数据\EEDM提取数据\480s切片数据含EP\' imname])
%         save_name = strcat(num2str(i),'_',num2str(j),'.mat');
%         savepath = [['C:\Users\Administrator\Desktop\毕业设计\数据\初始数据\EEDM提取数据\480s切片数据含EP\' imname],'\',save_name];
%         save(savepath,'co2','ECG_AM','PPG_AM','HRV','PRV','reference_value','ECG_op','PPG_op')
        
    end
    
%     savepath = ['C:\Users\Administrator\Desktop\毕业设计\数据\初始数据\EEDM提取数据\480_分钟呼吸率','\',name];
%     save(savepath,'ECG_breathnum','PPG_breathnum','ECGAM_breathnum','PPGAM_breathnum','HRV_breathnum','PRV_breathnum','ECG_op_breathnum','PPG_op_breathnum','reference')
   

    ECG_MAE(i) = mae(ECG_breathnum,reference);
    PPG_MAE(i) = mae(PPG_breathnum,reference);
    HRV_MAE(i) = mae(HRV_breathnum,reference);
    PRV_MAE(i) = mae(PRV_breathnum,reference);
    ECGAM_MAE(i) = mae(ECGAM_breathnum,reference);
    PPGAM_MAE(i) = mae(PPGAM_breathnum,reference);
    ECG_op_MAE(i)= mae(ECG_op_breathnum,reference);
    PPG_op_MAE(i)= mae(PPG_op_breathnum,reference);
        
    

end

%%

        ECG_MAE = sort(ECG_MAE);
        PPG_MAE = sort(PPG_MAE);
        HRV_MAE = sort(HRV_MAE);
        PRV_MAE = sort(PRV_MAE);
        ECGAM_MAE = sort(ECGAM_MAE);
        PPGAM_MAE = sort(PPGAM_MAE);
        ECG_op_MAE = sort(ECG_op_MAE);
        PPG_op_MAE = sort(PPG_op_MAE);
        
        
        ECGa = [ECG_MAE(11),ECG_MAE(21),ECG_MAE(32),mean(ECG_MAE)];
        PPGa = [PPG_MAE(11),PPG_MAE(21),PPG_MAE(32),mean(PPG_MAE)];
        HRVa = [HRV_MAE(11),HRV_MAE(21),HRV_MAE(32),mean(HRV_MAE)];
        PRVa = [PRV_MAE(11),PRV_MAE(21),PRV_MAE(32),mean(PRV_MAE)];
        ECGAMa = [ECGAM_MAE(11),ECGAM_MAE(21),ECGAM_MAE(32),mean(ECGAM_MAE)];
        PPGAMa = [PPGAM_MAE(11),PPGAM_MAE(21),PPGAM_MAE(32),mean(PPGAM_MAE)];
        ECG_op_MAEa = [ECG_op_MAE(11),ECG_op_MAE(21),ECG_op_MAE(32),mean(ECG_op_MAE)];
        PPG_op_MAEa = [PPG_op_MAE(11),PPG_op_MAE(21),PPG_op_MAE(32),mean(PPG_op_MAE)];


%%
% Testing the AE and MAE value
% maey = zeros(1,42);
% for i = 1:42
%    
%    name=Files(i).name;
%    folder=Files(i).folder;
%    s = load([folder,'\',name]); 
%    reference = s.rs_ref_y;
%    
%    maey(i)=mae(refrr(i,:),reference); 
%     
%     
% end

