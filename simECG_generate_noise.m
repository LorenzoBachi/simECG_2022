function multileadNoise = simECG_generate_noise(ecgLength, noiseType, noiseRMS, ESTParameters)
%
% multileadNoise = simECG_generate_noise() returns physiological noises obtained
% from the MIT-BIH Noise Stress Test Database. The signals of the
% the MIT–BIH Noise Stress Test database is used as a noise source for
% orthogonal leads X and Y, whereas noise for lead Z is obtained as the sum
% of squares of these two noises of leads X and Y. Then, an inverse of
% Dowers’s transformation matrix is applied to generate noise in the remaining
% 12 ECG leads. Ultimately, the noise is rescaled to the desired RMS value.
%
% noiseType: value from 0 to 4
% 0 - no noise added (noise RMS = 0 mV)
% 1 - motion artefacts
% 2 - electrode movement
% 3 - baseline wander
% 4 - mixture of noises
% 5 -
% 6 - Exercise stress test noise (from R. Bailón) %CPerez 03/2022 

disp('Generating noise ...');

noiseLength = 1805556; % Noise length
if ecgLength < noiseLength
    noiseStart = randi([1 (noiseLength-ecgLength)]); % take starting point randomly
else
    noiseStart = randi([1 noiseLength/2]); %1805556 - cut half of length noise segments
    cycles = ceil(ecgLength/(noiseLength/2));
    noiseTemp = [];
end

switch noiseType
    case 0     % no noise added
        multileadNoise = zeros(15, ecgLength);
    case 1     % motion artefacts
        data = load('DATA_noises_real');
        if ecgLength < noiseLength
            multileadNoise = data.motion_artefacts(:, noiseStart+1:noiseStart+ecgLength);
        else
            multileadNoise = data.motion_artefacts(:, noiseStart+1:noiseStart+(noiseLength/2));
            for i = 1:cycles
                noiseTemp = [noiseTemp multileadNoise];
            end
            multileadNoise = noiseTemp(:,1:ecgLength);
        end
        
    case 2     % electrode movement
        data = load('DATA_noises_real');
        if ecgLength < noiseLength
            multileadNoise = data.electrode_movement(:, noiseStart+1:noiseStart+ecgLength);
        else
            multileadNoise = data.electrode_movement(:, noiseStart+1:noiseStart+(noiseLength/2));
            for i = 1:cycles
                noiseTemp = [noiseTemp multileadNoise];
            end
            multileadNoise = noiseTemp(:,1:ecgLength);
        end
        
    case 3     % baseline wander
        data = load('DATA_noises_real');
        if ecgLength < noiseLength
            multileadNoise = data.baseline_wander(:, noiseStart+1:noiseStart+ecgLength);
        else
            multileadNoise = data.baseline_wander(:, noiseStart+1:noiseStart+(noiseLength/2));
            for i = 1:cycles
                noiseTemp = [noiseTemp multileadNoise];
            end
            multileadNoise = noiseTemp(:,1:ecgLength);
        end
        
    case 4     % mixture of noises
        data = load('DATA_noises_real');
        if ecgLength < noiseLength
            multileadNoise = data.mixture_of_noises(:, noiseStart+1:noiseStart+ecgLength);
        else
            multileadNoise = data.mixture_of_noises(:, noiseStart+1:noiseStart+(noiseLength/2));
            for i = 1:cycles
                noiseTemp = [noiseTemp multileadNoise];
            end
            multileadNoise = noiseTemp(:,1:ecgLength);
        end
        
    case 5     % baseline wander and motion artefacts because electrode motion has resudual ECG 10/02/2021 Alba
        data = load('DATA_noises_real');
        if ecgLength < noiseLength
            multileadNoise = data.baseline_wander(:, noiseStart+1:noiseStart+ecgLength) + data.motion_artefacts(:, noiseStart+1:noiseStart+ecgLength);
        else
            multileadNoise = data.baseline_wander(:, noiseStart+1:noiseStart+(noiseLength/2)) + data.motion_artefacts(:, noiseStart+1:noiseStart+(noiseLength/2));
            for i = 1:cycles
                noiseTemp = [noiseTemp multileadNoise];
            end
            multileadNoise = noiseTemp(:,1:ecgLength);
        end
        
    case 6     % Exercise stress test noise %CPerez 04/2022
        load('DATA_noises_EST_real');
%         sigNum = randi([1 25])  % Select randomly a stress test noise from 25 possible
        sigNum = 17
        noise = simECG_construct_real_EST_noise(EST_real_noise(sigNum).noise); %15 leads
        noisePeak = EST_real_noise(sigNum).peak;
        noiseFs = EST_real_noise(sigNum).fs;
        clear DATA_noises_ST_real
        
        %Exercise
        if noisePeak/noiseFs > ESTParameters.peak
            noise_e = noise(:,noisePeak-ESTParameters.peak*1e3:noisePeak -1);     
        else
            diff = round(ESTParameters.peak*1e3-noisePeak);
            noiseIni = repmat(noise(:,1:60*noiseFs),1,ceil(diff/(60*noiseFs))); %samples from the begining
            add_e = noiseIni(:,1:diff); %samples from the begining
            noise_e = [add_e noise(:,1:noisePeak)];
        end

        %Recovery
        if (size(noise,2) - noisePeak)/noiseFs > (ecgLength/1000 - ESTParameters.peak)
            noise_r = noise(:,noisePeak+1:noisePeak+round(ecgLength - ESTParameters.peak*1e3));     
        else
            diff = (ecgLength - ESTParameters.peak*1e3) - length(noise(1,noisePeak+1:end));
            noiseEnd = repmat(noise(:,size(noise,2)-60*noiseFs+1:size(noise,2)),1,ceil(diff/(60*noiseFs))); %samples belong to the end
            add_r = noiseEnd(:,1:diff); %samples belong to the end
            noise_r = [noise(:,noisePeak+1:end) add_r];
        end
        
        %Exercise-Stress-Test noise
        multileadNoise=[noise_e noise_r].*1e-3; %in mVolts        
end

if noiseType > 0
    % Adjust to desired noise RMS value
    for i = 1:15
        multileadNoise(i,:) = noiseRMS*(multileadNoise(i,:)/std(multileadNoise(i,:)));
    end
end
end