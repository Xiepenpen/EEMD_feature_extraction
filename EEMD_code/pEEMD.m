function imf = pEEMD(data,FsOrT,Nstd,NE)
% plot EEMD decomposition chart
%%
if length(FsOrT) == 1  %if input is frequency value
    t = 1/FsOrT:1/FsOrT:length(data)/FsOrT;
else
    t = FsOrT;         %if input is time vector
end
imf=kEEMD(data,Nstd,NE);
%imf=eemd(data,Nstd,NE)';  
rows = size(imf,1);    %IMF layers number

figure('Name','EEMD Decomposition figure','Color','white');
subplot(rows+1,1,1);
plot(t,data,'k');grid on;
xlim([t(1) t(end)]);
ylabel('Original Signal');
title('EEMD Decomposition');

for i = 1:size(imf,1)
    subplot(rows+1,1,i+1);
    plot(t,imf(i,:),'k');
    xlim([t(1) t(end)]);
    ylabel(['IMF',num2str(i)]);
    if (i == size(imf,1))
        ylabel(['res']);
        xlabel('time');
    end
    grid on;
end
end

