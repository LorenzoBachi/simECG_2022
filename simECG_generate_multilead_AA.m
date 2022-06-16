function multileadAA = simECG_generate_multilead_AA(targets_beats, QRSindex, fibFreqz, realAAon, ecgLength, AFburden, ecgParameters)
% multileadAA = simECG_gen_multilead_AA() returns multilead (15 lead) atrial
% activity. 
%
% Generated leads:
% multileadAA(1,:) - I      multileadAA(7,:) - V1    multileadAA(13,:) - X     
% multileadAA(2,:) - II     multileadAA(8,:) - V2    multileadAA(14,:) - Y  
% multileadAA(3,:) - III    multileadAA(9,:) - V3    multileadAA(15,:) - Z  
% multileadAA(4,:) - aVR    multileadAA(10,:) - V4 
% multileadAA(5,:) - aVL    multileadAA(11,:) - V5 
% multileadAA(6,:) - aVF    multileadAA(12,:) - V6 

%Last update: CPerez 05/2022

% Decide which rhythm type to generate
if AFburden == 0
    rhythmType = 0; % SR
elseif AFburden == 1
    rhythmType = 1; % AF
elseif AFburden > 0 && AFburden < 1
    rhythmType = 2; % PAF
else
    error('AF burden must be a value between 0 and 1')
end

disp('Generating atrial activity ...');

multileadAA = zeros(15, ecgLength);
% Generate atrial activity
switch rhythmType
    case 0 % Entire signal is SR
       if realAAon == 0
           Nrr = length(targets_beats);
           P_waves = simECG_generate_multilead_P_waves(Nrr); 
           
           % Replace P waves during ectopic beats
           numABPs = length(find(targets_beats == 3));
           if numABPs > 0
               P_waves_APBs = simECG_generate_multilead_P_waves(numABPs);
               iAPB = 1;
               for n = 1:Nrr 
                   if targets_beats(n) == 3
                       P_waves(:,n,:) = P_waves_APBs(:,iAPB,:);
                       iAPB = iAPB + 1;
                   end
               end
           end
           
           % Insert P waves% CPerez 03/2022
           Fs = 1000; 
           RR = diff(QRSindex./Fs); %seconds
           QR = 50;
           for p_num = 2:Nrr
               bias = 0;
               if ~isfield(ecgParameters,'peak') || (QRSindex(p_num) <= ecgParameters.peak*1e3) %Exercise or normal
                   if RR(p_num-1) < 0.520 %Change point
                       PQ = round(152 + 358*(RR(p_num-1)-0.520));%-34.53 + 361.21*RR
                   else
                       PQ = randi([140 160],1); %constant 152
                   end
               end
               
               if isfield(ecgParameters,'peak') && QRSindex(p_num) > ecgParameters.peak*1e3 %Recovery
                   if RR(p_num-1) < 0.430 %Change point
                       PQ = round(-65 + 470*RR(p_num-1)); %-65.26 + 470.24*RR
                   else
                       PQ = randi([130 150],1); %constant 139
                   end
               end
               
               lenP = length(QRSindex(p_num)-(PQ + QR -1):QRSindex(p_num)-QR);
               bias = lenP - size(P_waves,3);
               if bias < 0
                   multileadAA(:, QRSindex(p_num)-(PQ + QR -1 -bias):QRSindex(p_num)- QR) = P_waves(:,p_num,:);
               else
                   multileadAA(:, QRSindex(p_num)-(PQ + QR -1):QRSindex(p_num)- (QR + bias)) = P_waves(:,p_num,:);
               end
               %                else
               %                    multileadAA(:, QRSindex(p_num)-249:QRSindex(p_num)-100) = P_waves(:,p_num,:);
               %                end
           end
       end   
       
       
       
    case 1 % Entire signal is AF
        if realAAon == 0
            multileadAA = simECG_generate_multilead_f_waves(fibFreqz, ecgLength);
        else
            sigNum = randi([1 20]);
            load('DATA_f_waves_real')
            f_waves = DATAfWaves(sigNum).f_waves;
            
            clear DATA_f_waves_real
            
            if ecgLength > length(f_waves)   % Concatenate if RR series is shorter than the desired length
                nCycles = ceil(ecgLength/length(f_waves));
                multileadAA = [];
                for i = 1:nCycles
                    multileadAA = [multileadAA f_waves];
                end
                multileadAA = multileadAA(:,1:ecgLength);  
            else
                fStart = randi([1 (length(f_waves)-ecgLength)]); % Randomly select start position
                multileadAA = f_waves(:,fStart+1:fStart+ecgLength);
            end   
        end
                
    case 2 % PAF
        Nrr = length(targets_beats);
        if realAAon == 0  % Synthetic atrial activity is prefered
           P_waves = simECG_generate_multilead_P_waves(Nrr); 
           % Replace P waves during ectopic beats
           numABPs = length(find(targets_beats == 3));
           if numABPs > 0
               P_waves_APBs = simECG_generate_multilead_P_waves(numABPs);
               iAPB = 1;
               for n = 1:Nrr 
                   if targets_beats(n) == 3
                       P_waves(:,n,:) = P_waves_APBs(:,iAPB,:);
                       iAPB = iAPB + 1;
                   end
               end
           end
           % Generate f waves
           f_waves = simECG_generate_multilead_f_waves(fibFreqz, ecgLength);
           % Insert P waves
           for p_num = 2:Nrr
               multileadAA(:, QRSindex(p_num)-249:QRSindex(p_num)-100) = P_waves(:,p_num,:);  
           end
           % Insert f waves
           for p_num = 2:Nrr
               if targets_beats(p_num) == 2
                   praIndex = QRSindex(p_num);
                   while ( targets_beats(p_num) == 2 ) && ( p_num < Nrr )
                      p_num = p_num + 1;
                   end 
                   pabIndex = QRSindex(p_num);
                   multileadAA(:, praIndex:pabIndex) = f_waves(:,praIndex:pabIndex);   
               end
           end
               
        else % Real atrial activity is prefered
            sigNum = randi([1 20]);
            load('DATA_f_waves_real')
            f_waves = DATAfWaves(sigNum).f_waves;
            
            clear DATA_f_waves_real
            
            if ecgLength > length(f_waves)   % Concatenate if RR series is shorter than the desired length
                nCycles = ceil(ecgLength/length(f_waves));
                multileadfWaves = [];
                for i = 1:nCycles
                    multileadfWaves = [multileadfWaves f_waves];
                end
                multileadfWaves = multileadfWaves(:,1:ecgLength);  
            else
                fStart = randi([1 (length(f_waves)-ecgLength)]); % Randomly select start position
                multileadfWaves = f_waves(:,fStart+1:fStart+ecgLength);
            end 
            for p_num = 2:Nrr
                if targets_beats(p_num) == 2
                    praIndex = QRSindex(p_num);
                    while targets_beats(p_num) == 2 && p_num < Nrr
                        p_num = p_num + 1;
                    end 
                    if p_num>length(QRSindex)
                        temp=1;
                    end
                    pabIndex = QRSindex(p_num);
                    multileadAA(:, praIndex:pabIndex) = multileadfWaves(:,praIndex:pabIndex); 
                end
            end
            
        end
             
end

