function [simECGdata, initialParameters, annotations, ecgParameters] = simECG_generator(sigLength, realRRon, realVAon, realAAon, noiseType, noiseRMS, onlyRR, arrhythmiaParameters, ecgParameters)
% [] = simECG_generator() returns a 15-by-N matrix containing 15 lead
% ECGs. Three types of ECG signals can be generated: SR (AF burden set to 0, 
% AF (AF burden set to 1) or PAF (AF burden any value from the interval 
% (0, 1)). Standard leads I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6
% and Frank leads X, Y, Z are generated(sampling frequence 1000 Hz).
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
% rrLength indicates the length of the desired ECG signal (in RR intervals)
%
% realRRon 1 indicates that real RR intervals are used, 0 - synthetic
% realVAon 1 indicates that real ventricular activity is used, 0 - synthetic
% realAAon 1 indicates that real atrial activity is used, 0 - synthetic
%
% onlyRR 1 - only RR intervals are generated, 0 - multilead ECG is generated
%
% noiseType: a number from 0 to 4
% 0 - no noise added (noise RMS = 0 mV)
% 1 - motion artefacts
% 2 - electrode movement artefacts
% 3 - baseline wander
% 4 - mixture of type 1, type 2 and type 3 noises
%
% noiseLevel - noise level in milivolts, i.e. 0.02 corresponds to 20 uV
%
% arrhythmiaParameters is a structure containing the following fields:
% - AFburden is a value between 0 and 1. 0: the entire signal is SR, 
% 1: the entire signal is AF.
% - stayInAF denotes the probability to stay in AF
%
% Output arguments:
% simECGdata returns generated data (multilead ECG, multilead ventricular
% activity, multilead atrial activity, QRS index, etc.). initialParameters 
% returns initial parameter values used to generated ECG signals.   
%
% Known problems:
% AV node model used for generating RR intervals during AF is relatively slow.
%
% Synthetic P wave amplitude is nearly 1.5 lower in several leads than that
% observed in reality (at least for healthy patients). Parameters for 
% simulating Type 2 P waves are taken from the paper by Havmoller et al. 
% Age-related changes in P wave morphology in healthy subjects. 
% BMC Cardiovascular Disorders, 7(1), 22, 2007.
%
% Interpolated TQ intervals (using a cubic spline interpolation) sometimes 
% do not look realistic.

disp('    ECG generator: simulation starting ...');

switch onlyRR
    case 1 % only RR intervals are generated
        % Generate initial parameters (fibrillatory frequency)
        fibFreqz = simECG_fibrillation_frequency();
        % Generate RR intervals
        [rr,annTime,annType,annRhythm,targets_beats,pafBoundaries, pafEpisLength, ecgParameters,state_history] = simECG_global_rr_intervals(sigLength, fibFreqz, realRRon, arrhythmiaParameters, ecgParameters);
        rr(cumsum(rr)>sigLength)= [];
        rrLength = numel(rr);
        simECGdata.rr = rr;
        simECGdata.multileadECG = [];
        simECGdata.multileadVA = [];
        simECGdata.multileadAA = [];
        simECGdata.multileadNoise = [];
        simECGdata.QRSindex = [];
        simECGdata.targets_beats = targets_beats;
        simECGdata.pafBoundaries = pafBoundaries;
        simECGdata.pafEpisLength = pafEpisLength;
        simECGdata.ecgLength = [];

        ecgParameters.fibFreqz = fibFreqz;
        ecgParameters.rrLength = rrLength;
        ecgParameters.realRRon = realRRon;
        ecgParameters.realVAon = realVAon;
        ecgParameters.realVAon = realAAon;
        ecgParameters.noiseType = noiseType;
        ecgParameters.noiseRMS = noiseRMS;
        
        annotations.annTime = annTime;
        annotations.annType = annType;
        annotations.annRhythm = annRhythm;
        
    case 0 % multilead ECG is generated
        % Check for errors:
        if (realVAon == 0) && (realAAon == 1)
            msg = ('Selection of synthetic ventricular activity and real atrial activity is not allowed');
            error('MyComponent:incorrectType', msg);
        end

        % Generate initial parameters (fibrillatory frequency)
        fibFreqz = simECG_fibrillation_frequency();   
        % Generate RR intervals
        [rrIn,annTime,annType,annRhythm,targets_beats,pafBoundaries, pafEpisLength, ecgParameters,state_history] = simECG_global_rr_intervals(sigLength,fibFreqz, realRRon, arrhythmiaParameters, ecgParameters);
        rrLength = numel(rrIn);
        % Generate multilead ventricular activity
        [QRSindex, TendIndex,rr, multileadVA, ecgLength] = simECG_generate_multilead_VA(rrLength, targets_beats, rrIn, realVAon, ecgParameters,state_history); %_QTC adde by Alba 19/03
        % Generate multilead atrial activity
        multileadAA = simECG_generate_multilead_AA(targets_beats, QRSindex, fibFreqz, realAAon, ecgLength, arrhythmiaParameters.AFburden, ecgParameters);
        % Generate multilead noise
        multileadNoise = simECG_generate_noise(ecgLength, noiseType, noiseRMS, ecgParameters);
        % Generate multilead noise
        multileadECG = multileadVA + multileadAA + multileadNoise;

        simECGdata.rr = rr;
        simECGdata.multileadECG = multileadECG;
        simECGdata.multileadVA = multileadVA;
        simECGdata.multileadAA = multileadAA;
        simECGdata.multileadNoise = multileadNoise;
        simECGdata.QRSindex = QRSindex;
        simECGdata.TendIndex = TendIndex;
        simECGdata.targets_beats = targets_beats;
        simECGdata.pafBoundaries = pafBoundaries;
        simECGdata.pafEpisLength = pafEpisLength;
        simECGdata.state_history = state_history;
        simECGdata.ecgLength = ecgLength';
        simECGdata.Fr = ecgParameters.Fr';
               

        initialParameters.fibFreqz = fibFreqz;
        initialParameters.rrLength = rrLength;
        initialParameters.realRRon = realRRon;
        initialParameters.realVAon = realVAon;
        initialParameters.realVAon = realAAon;
        initialParameters.noiseType = noiseType;
        initialParameters.noiseRMS = noiseRMS;
        
        annotations.annTime = annTime;
        annotations.annType = annType;
        annotations.annRhythm = annRhythm;
        
end

end