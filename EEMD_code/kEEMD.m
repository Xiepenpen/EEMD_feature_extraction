function imf = kEEMD(data,Nstd,NE)
%% Get the current version number
vtemp = version;
vv = {vtemp(strfind(vtemp,'R')+1:strfind(vtemp,'R')+4),vtemp(strfind(vtemp,'R')+5)};

%% Demo usage of setting random seeds
% Random seeds fix the white noise added, setting random seeds ensures consistent results from one decomposition to the next
rng(123)
%%
imfTemp = eemd(data,Nstd,NE)';  % Note that imf(1,:) here is the original signal
imf = imfTemp(2:end,:);
end