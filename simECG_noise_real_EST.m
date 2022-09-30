function [noiseAll_15] = simECG_noise_real_EST(ecgLength, ecgParameters)
% noiseAll_15 = simECG_noise_real_EST() returns one real Exercise Stress Test noise signal adapted to the ECG length.
%
% Copyright (c), Cristina Perez, University of Zaragoza, 09/2022

% load('DATA_noises_MN_real_L3');
load('DATA_noises_EST_real_L3');


sigNum = randi([1 25])  % Select randomly a stress test noise from 25 possible
% sigNum = 8
noise = EST_real_noise(sigNum).noise; %3 leads
noisePeak = EST_real_noise(sigNum).peak;
noiseFs = EST_real_noise(sigNum).fs;

%Exercise
if noisePeak/noiseFs > ecgParameters.peak
    noise_e = noise(:,noisePeak-ecgParameters.peak*1e3:noisePeak -1);
else
    diff = round(ecgParameters.peak*1e3-noisePeak);
    noiseIni = repmat(noise(:,1:60*noiseFs),1,ceil(diff/(60*noiseFs))); %samples from the begining
    add_e = noiseIni(:,1:diff); %samples from the begining
    noise_e = [add_e noise(:,1:noisePeak)];
end

%Recovery
if (size(noise,2) - noisePeak)/noiseFs > (ecgLength/1000 - ecgParameters.peak)
    noise_r = noise(:,noisePeak+1:noisePeak+round(ecgLength - ecgParameters.peak*1e3));
else
    diff = (ecgLength - ecgParameters.peak*1e3) - length(noise(1,noisePeak+1:end));
    noiseEnd = repmat(noise(:,size(noise,2)-60*noiseFs+1:size(noise,2)),1,ceil(diff/(60*noiseFs))); %samples belong to the end
    add_r = noiseEnd(:,1:diff); %samples belong to the end
    noise_r = [noise(:,noisePeak+1:end) add_r];
end

noiseAll = [noise_e noise_r].*1e-3; %in mVolts



% Transform to the 15 leads
%1)Obtain augmented unipolar limb leads
noiseAll_8 = leadcalc(noiseAll,'stan');% V1,V2,V3,V4,V5,V6,I,II
noiseAll_12 = leadcalc(noiseAll_8,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III

noiseAll_15=vertcat(noiseAll_8(7,:),noiseAll_8(8,:),noiseAll_12(12,:),...
    -noiseAll_12(9,:),-noiseAll_12(7,:),noiseAll_12(11,:),noiseAll_8(1:6,:),noiseAll);

end