function [noiseAll_15] = simECG_noise_real_EST(ecgLength, ecgParameters, noiseRMS, signal)
% noiseAll_15 = simECG_noise_real_EST() returns one real Exercise Stress Test noise signal adapted to the ECG length.
%
% Copyright (c), Cristina Perez, University of Zaragoza, 09/2022

load('DATA_MNnoises_EST_real');

fs = ecgParameters.fs;
sigNum = randi([1 25]);  % Select randomly a stress test noise from 25 possible
% sigNum = 8
noise = EST_real_noise(sigNum).MNnoise; %3 leads
noisePeak = EST_real_noise(sigNum).peak;
noiseFs = EST_real_noise(sigNum).fs;

%Exercise
if noisePeak/noiseFs > ecgParameters.peak
    noise_e = noise(:,noisePeak-ecgParameters.peak*1e3:noisePeak);
else
    diff = round(ecgParameters.peak*1e3-noisePeak);
    noiseIni = repmat(noise(:,1:60*noiseFs),1,ceil(diff/(60*noiseFs))); %samples from the begining
    add_e = noiseIni(:,1:diff); %samples from the begining
    noise_e = [add_e noise(:,1:noisePeak)];
end

%Recovery
if (size(noise,2) - noisePeak)/noiseFs > (ecgLength/1000 - ecgParameters.peak)
    noise_r = noise(:,noisePeak:noisePeak+round(ecgLength - ecgParameters.peak*1e3));
else
    diff = (ecgLength - ecgParameters.peak*1e3) - length(noise(1,noisePeak+1:end));
    noiseEnd = repmat(noise(:,size(noise,2)-60*noiseFs+1:size(noise,2)),1,ceil(diff/(60*noiseFs))); %samples belong to the end
    add_r = noiseEnd(:,1:diff); %samples belong to the end
    noise_r = [noise(:,noisePeak:end) add_r];
end

noiseAll = [noise_e noise_r].*1e-3; %in mVolts


if noiseRMS > 1 %info in dB
    DefnoiseRMS = noiseRMS;
    Rpos = round(cumsum(ecgParameters.RR).*fs); %in samples
    noiseRMS = zeros(1,size(signal,1));ppQRS = zeros(1,size(signal,1));
    for l = 1:size(signal,1)
        qrsAll = zeros(100,length(Rpos));
        for ii = 2:length(Rpos)-1
            qrsAll(:,ii) = signal(l,Rpos(ii)-50:Rpos(ii)+50-1);
        end
        Aqrs = mean(qrsAll,2);
        ppQRS(l) = max(Aqrs)-min(Aqrs);
        noiseRMS(l) = ((max(Aqrs)-min(Aqrs))/(10^(DefnoiseRMS/20)));%in mV
    end     
end

%--> Adapt real Noise
%             factors = noiseRMS./mean(abs(realNoise(1:60*fa,:)),1);
factors = noiseRMS'./(mean(abs(noise(:,noisePeak-30*noiseFs+1:noisePeak+30*noiseFs,:)),2).*1e-3); %in mVolts
realNoise = noiseAll.*factors;
realNoise = realNoise(:,1:ecgLength);


% Transform to the 15 leads
%1)Obtain augmented unipolar limb leads
noiseAll_12 = leadcalc(realNoise,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III
noiseAll_xyz = leadcalc(realNoise,'hcuzst');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III
noiseAll_15=vertcat(realNoise(7:9,:),-noiseAll_12(9,:),noiseAll_12(7,:),noiseAll_12(11,:),...
    realNoise(1:6,:),noiseAll_xyz);

end