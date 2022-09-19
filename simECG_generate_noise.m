function [multileadNoise, poles] = simECG_generate_noise(ecgLength, noiseType, noiseRMS, ecgParameters)
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
% 5 - bw + ma because em has residual ECG;
% 6 - Simulated Muscular Noise %CPerez 07/2022 
% 7 - Real Exercise stress test noise (from R. Bailón)
% 8 - Motion artifacts

disp('Generating noise ...');

poles=[];

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
        
    case 6     %Simulated Muscular Noise
        [multileadNoise, poles] = simECG_muscular_noise(ecgLength,ecgParameters, noiseRMS);%in mVolts
        
    case 7     %Real Exercise stress test noise (from R. Bailón)
        multileadNoise = simECG_noise_real_EST(ecgLength,ecgParameters);
        
    case 8     %Motion artifacts
        multileadNoise = simECG_generate_motion_artifact(ecgParameters.MA_Prob, ecgParameters.MA_sigBernoGauss, ecgParameters.MA_Flag, noiseRMS, ecgLength);
               
end

if noiseType > 0 && (noiseType ~=6 || noiseType ~=8)
    % Adjust to desired noise RMS value
    for i = 1:15
        multileadNoise(i,:) = noiseRMS*(multileadNoise(i,:)/std(multileadNoise(i,:)));
    end
end
end
