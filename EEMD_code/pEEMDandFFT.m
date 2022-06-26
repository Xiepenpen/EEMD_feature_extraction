function imf = pEEMDandFFT(y,FsOrT,Nstd,NE)
% Plotting the EEMD decomposition of the signal against the spectrum of each IMF component
% Input.
% y is the signal to be decomposed
% FsOrT is the sampling frequency or the sampling time vector, if the sampling frequency, the variable is entered as a single value; if the time vector, the variable is a one-dimensional vector of the same length as y
% Nstd is the ratio of the standard deviation of the additional noise to the standard deviation of y
% NE is the number of times the signal is averaged
%%

if length(FsOrT) == 1
    t = 1/FsOrT:1/FsOrT:length(y)/FsOrT;
    Fs = FsOrT;
else
    t = FsOrT;
    Fs = 1/(t(2)-t(1));
end
imf = kEEMD(y,Nstd,NE);
figure
subplot(size(imf,1)+1,2,1);
plot(t,y,'k');grid on;
ylabel('Original data');
title('EEMD decomposition');
set(gca,'XTick',[]);
subplot(size(imf,1)+1,2,2);
pFFT(y,Fs);grid on;
title('Corresponding spectrum');
set(gca,'XTick',[]);
for i = 2:size(imf,1)+1
    subplot(size(imf,1)+1,2,i*2-1);
    plot(t,imf(i-1,:),'k');
    ylabel(['IMF',num2str(i-1)]);
    if (i ~= size(imf,1)+1)
        set(gca,'XTick',[]);
    end
    if (i == size(imf,1)+1)
        ylabel(['res']);
        xlabel('time/s');
    end
    grid on;
    subplot(size(imf,1)+1,2,i*2);
    pFFT(imf(i-1,:),Fs);
    if (i ~= size(imf,1)+1)
        set(gca,'XTick',[]);
    end
    if (i == size(imf,1)+1)
        xlabel('frequency/Hz');
    end
    grid on;
end
end
