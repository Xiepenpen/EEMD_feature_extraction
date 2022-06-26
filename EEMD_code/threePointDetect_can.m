function [breath,save_peaks,locs_peak,save_valleys,locs_valley,zero_point] = three_point_detect(data,t,refx)
%data为输入信号
%t是归一化的横坐标 从0s开始到480s

for m = 1:length(refx)
    %按照refx指定时间切片
    index = 0;
    for m1=2:length(t)-1
        if t(m1-1)<refx(m) & t(m1+1)>refx(m)
            index = m1;
        elseif t(m1+1)<refx(m) & m>=round((length(refx)/40)*39)
            index = length(t);
        end
    end
    data_s = data(1:index);
    
%     %Use mypeaks function to detect the peaks
      mu = mean(data_s);
%     sig = std(data_s);
%     %peak探测的上下高度限制 e.x.不采集大于三倍本信号std的峰值
%     data_s(abs(data_s-mu) > 3*sig) = NaN;
%     %峰值探测窗口默认设置为60s
%     m_peaks = mypeaks(data_s,4,2,5).';
%     m_valley = mypeaks(-data_s,4,2,5).';
    %[m_peaks,locs_peak] = findpeaks(data,'minpeakheight',mu);
    %[m_valley,locs_valley] = findpeaks(-data,'minpeakheight',mu);
    [m_peaks,locs_peak] = findpeaks(data,t);
    [m_valley,locs_valley] = findpeaks(-data,t);
    save_peaks = m_peaks;
    save_valleys = m_valley;
    mean_gap = mean(diff(locs_peak));
    %calculate the zero point 
    t_data = t;
    t_peaks = locs_peak;%以480s为时间轴的峰值数据横坐标
    t_valley = locs_valley;

    zeropoint = 0;
    t_zeropoint = zeros(1,length(data_s));
    for i = 2:length(data_s)

        if (data_s(i-1) >=mu && data_s(i)<=mu)||(data_s(i-1) <=mu && data_s(i)>=mu)
            zeropoint = zeropoint+1;
            t_zeropoint(zeropoint) = t_data(i)+(t_data(i)-t_data(i-1))/2;
        end
    
    end

    %剔除数组中的0（多余位） 得到伪零点t坐标
    t_zeropoint(t_zeropoint==0)=[];
    zero_point = t_zeropoint;
    breath_num = 0;
    %breath_num = length(t_zeropoint)/2+1;
    for i = 2:2:length(t_zeropoint)
        if i==2
            max = t_zeropoint(i);
            min = 0;
        elseif i==length(t_zeropoint)
            max = t_zeropoint(i)+10;
            min = t_zeropoint(i-2);
        else
            max = t_zeropoint(i);
            min = t_zeropoint(i-2);
        end

        has_peak = 0;
        has_valley = 0;
        for j = 1:length(t_peaks)
           if t_peaks(j)<=max & t_peaks(j)>=min 
               has_peak = has_peak+1;
           end
        end

        for k = 1:length(t_peaks)
           if t_valley<=max & t_valley>=min
               has_valley = has_valley+1;
           end
        end

        if has_peak ~=0 | has_valley~=0
           breath_num = breath_num +1; 
        end
    
    end

    breath(m) = (breath_num/refx(m))*60;

end

end