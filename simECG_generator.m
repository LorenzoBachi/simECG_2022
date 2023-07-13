function [simECGdata, initialParameters, annotations] = simECG_generator(sigLength, realRRon, realVAon, realAAon, noiseType, noiseRMS, onlyRR, arrhythmiaParameters, simECGdata)
% [] = simECG_generator() returns a 15-by-N matrix of 15-lead ECG. Standard
% leads I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6 and Frank leads
% X, Y, Z are generated, with a sampling frequency of 1000 Hz. Please note
% that simulated RR intervals, simulated ventricular activity and simulated
% atrial activity allow for the maximum control over simulation outcomes.
% This functions also returns annotations of the simulated record in the
% MIT-BIH Arrhythmia Database annotation style.
%
% Generated leads:
% multileadVA(1,:) - I      multileadVA(7,:) - V1    multileadVA(13,:) - X     
% multileadVA(2,:) - II     multileadVA(8,:) - V2    multileadVA(14,:) - Y  
% multileadVA(3,:) - III    multileadVA(9,:) - V3    multileadVA(15,:) - Z  
% multileadVA(4,:) - aVR    multileadVA(10,:) - V4 
% multileadVA(5,:) - aVL    multileadVA(11,:) - V5 
% multileadVA(6,:) - aVF    multileadVA(12,:) - V6 
%
% Input arguments:
% sigLength - desired ECG length of the simulated record, in seconds.
% realRRon - 0 for simulated RR intervals, 1 for real RR intervals.
% realVAon - 0 for simulated ventricular activity, 1 for real ventricular
% activity.
% realAAon - 0 for simulated atrial activity, 1 for real atrial activity.
% noiseType - an array containing the desired types of noise (eg [3, 6,
% 8]).
% noiseRms - the desired root mean square level of each desired noise type.
% onlyRR 1 - 0 for multilead ECG simulation, 1 for RR series-only
% simulation.
% arrhythmiaParameters - struct of arrhythmia simulation parameters defined
% in the main script.
% simECGdata - struct of ECG simulation parameters defined in the main
% script.
%
% Output arguments:
% simECGdata - struct of ECG simulation parameters defined in the main
% script (updated).
% initialParameters - struct of user-defined parameters as well as f waves
% fibrillatory frequency.
% annotations - struct of simulated ECG record annotations in the MIT-BIH
% Arrhythmia Database style.
%
% Known problems: AV node model used for generating RR intervals during AF
% is relatively slow.
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

disp('    ECG generator: simulation starting ...');

if simECGdata.MA_Prob > 1, simECGdata.MA_Prob = 1 / simECGdata.MA_Prob; end

switch onlyRR
    case 1 % only RR intervals are generated
        % Generate initial parameters (fibrillatory frequency)
        fibFreqz = simECG_fibrillation_frequency();
        % Generate RR intervals
        [rr,annotations,targets_beats,simECGdata,hrArray,state_history] = simECG_global_rr_intervals(sigLength, fibFreqz, realRRon, arrhythmiaParameters, simECGdata);
        rr(cumsum(rr)>sigLength)= [];
        rrLength = numel(rr);
        simECGdata.rr = rr;
        simECGdata.multileadECG = [];
        simECGdata.multileadVA = [];
        simECGdata.multileadAA = [];
        simECGdata.multileadNoise = [];
        simECGdata.QRSindex = [];
        simECGdata.targets_beats = targets_beats;
        simECGdata.ecgLength = sigLength;
        simECGdata.Fr = [];
        simECGdata.poles = [];
        simECGdata.state_history = state_history;
        simECGdata.hrArray = hrArray;
        
        initialParameters.fibFreqz = fibFreqz;
        initialParameters.rrLength = rrLength;
        initialParameters.realRRon = realRRon;
        initialParameters.realVAon = realVAon;
        initialParameters.realAAon = realAAon;
        initialParameters.noiseType = noiseType;
        initialParameters.noiseRMS = noiseRMS;
        
    case 0 % multilead ECG is generated
        % Check for errors:
        if (realVAon == 0) && (realAAon == 1)
            msg = ('Selection of synthetic ventricular activity and real atrial activity is not allowed');
            error('MyComponent:incorrectType', msg);
        end

        % Generate initial parameters (fibrillatory frequency)
        fibFreqz = simECG_fibrillation_frequency();   
        % Generate RR intervals
        [rrIn,annotations,targets_beats,simECGdata,hrArray,state_history] = simECG_global_rr_intervals(sigLength,fibFreqz, realRRon, arrhythmiaParameters, simECGdata);
        % Generate multilead ventricular activity
        [QRSindex, TendIndex,rr, multileadVA, ecgLength, simECGdata] = simECG_generate_multilead_VA(targets_beats, rrIn, realVAon, simECGdata,state_history); %_QTC adde by Alba 19/03
        % Generate multilead atrial activity
        multileadAA = simECG_generate_multilead_AA(targets_beats, QRSindex, fibFreqz, realAAon, ecgLength, arrhythmiaParameters.B_af, simECGdata);
        % Generate multilead noise
        for ii = 1:numel(noiseType)
            [multileadNoise_All(:,:,ii)] = simECG_generate_noise(ecgLength, noiseType(ii), noiseRMS(ii), simECGdata, multileadVA);
        end
        multileadNoise = sum(multileadNoise_All,3);
        % Generate multilead noise
        multileadECG = multileadVA + multileadAA + multileadNoise;

        simECGdata.rr = rr;
        simECGdata.multileadECG = multileadECG;
        simECGdata.multileadVA = multileadVA;
        simECGdata.multileadAA = multileadAA;
        simECGdata.multileadNoise = multileadNoise;
        simECGdata.multileadNoise_All = multileadNoise_All;
        simECGdata.QRSindex = QRSindex;
        simECGdata.TendIndex = TendIndex;
        simECGdata.targets_beats = targets_beats;
        simECGdata.state_history = state_history;
        simECGdata.ecgLength = ecgLength';
        simECGdata.Fr = simECGdata.Fr';
        simECGdata.state_history = state_history;
        simECGdata.hrArray = hrArray;
        
        initialParameters.fibFreqz = fibFreqz;
        initialParameters.rrLength = length(rr);
        initialParameters.realRRon = realRRon;
        initialParameters.realVAon = realVAon;
        initialParameters.realAAon = realAAon;
        initialParameters.noiseType = noiseType;
        initialParameters.noiseRMS = noiseRMS;
        
end

end