function [rr, annTime, annBeats, annRhythm, targets_beats, pafBoundaries, pafEpisodeLength, ecgParameters] = simECG_global_rr_intervals(sigLength, fibFreqz, realRRon, arrhythmiaParameters, ecgParameters)
% [] = simECG_global_rr_intervals() generates the global RR series, which
% may include sinus rhythm, atrial fibrillation or other arrhythmias.

% The dominant length of PAF episodes, and the percentage of the total time AF is
% present (AF burden) is determined by a first-order two-state Markov chain.
%
% An atrioventricular node model in which the ventricles are assumed to be
% activated by the arriving atrial impulses according to a Poisson process
% is used to generate realistic RR series during an episode of AF

% Ventricular rhythm during SR is simulated using RR interval generator
% proposed by in which both the impact of parasympathetic stimulation
% (respiratory sinus arrhythmia) and baroreflex regulation (Mayer waves) is
% modeled by a bimodal power spectrum

% In case of real data the entire MIT-BIH Normal Sinus Rhythm database,
% consisting of 18 long-term ECG recordings, was taken as a basis for the set
% of 18 RR interval series of sinus rhythm.
%
% The Long Term Atrial Fibrillation database was used to compose a set
% of AF rhythm. In total, 69 RR interval series were extracted out of
% 84 long-term ECG recordings. The remaining records were excluded since
% they did not meet the criteria of ?5000 RR intervals of AF.

% With respect to the previous contributions to the simulator, in this
% version the RR series is generated from four separate pacing
% processes: the sinus node pacing, the atrial ectopic pacing, the
% ventricular ectopic pacing and the atrial fibrillation pacing. The four
% stimulation processes are combined in a single RR series. The pacing
% processes may behave in different manners depending on the dominating
% rhythm.

% Additionally, the simulator outputs beat and rhythm annotations in the
% same format of the MIT-BIH Arrhythmia Database.

% The desired ECG signal length in number of RR intervals 
% rrLength = round(sigLength/0.2);
rrLength = sigLength;

% signal length in ms
sigLengthMs = sigLength * 1000;
ecgParameters.Fr = [];

% fetching simulation parameters
AFburden = arrhythmiaParameters.AFburden;
stayInAF = arrhythmiaParameters.stayInAF;
% beat codes for annotation production
beatCodes = ['N','a','V','+']; % please note: if new beat types are added, '+' should always be last
rhythmCodes = {'(N','(AFIB','(SVTA'};

% Choice of rhythm to generate
if AFburden == 0 && ecgParameters.ESTflag == 1 %Cris 03/2022
    rhythmType = 3; % SR and EST
elseif AFburden == 0
    rhythmType = 0; % SR
elseif AFburden == 1
    rhythmType = 1; % AF
elseif AFburden > 0 && AFburden < 1
    rhythmType = 2; % PAF
else
    error('AF burden must be a value between 0 and 1')
end

disp('Generating RR intervals ...');

switch rhythmType  % 0 - sinus rhythm, 1 - AF, 2 - PAF

    case 0 % The entire rhythm is SR
        if realRRon == 1 % Use real RR series
            rr_sr = simECG_get_real_RR_intervals(0, rrLength); % sinus rhythm rr series (real)
            hrMean = (1./diff(rr_sr))*60 ;
        else % Generate sinus activity with the McSharry spectral model
            [rr_sr,hrMean,ecgParameters] = simECG_generate_sinus_rhythm(rrLength,ecgParameters); % sinus rhythm rr series
%             rr_sr = rr_sr(1:rrLength);
        end
        rhythm_states = ones(1,sigLength);
        rr_af = Inf;
        
    case 1 % The entire rhythm is AF
        if realRRon == 1 % Use real RR series
            rr_af = simECG_get_real_RR_intervals(1, rrLength); %atrial fibrillation rr series
        else % Generate atrial fibrillation activity with the AV model of V. Corino
            rr_af = simECG_generate_AF_intervals(fibFreqz,rrLength); %atrial fibrillation rr series
            rr_af = rr_af(1:rrLength);
        end
        rhythm_states = 2*ones(1,sigLength);
        rr_sr = Inf;
        
    case 2 % PAF
        % Generate alternating SR and AF episodes
        goToNS = 1-stayInAF;
        goToAF = (goToNS*AFburden)/(1-AFburden);
        stayInNS = 1 - goToAF;
        
        trans = [stayInNS,1-stayInNS;1-stayInAF, stayInAF];
        [~,rhythm_states] = hmmgenerate(sigLength,trans,[0; 1]);
        
        % Generate sinus node pacing activity
        srLength = rrLength;
        if realRRon == 1 % Use real RR series
            rr_sr = simECG_get_real_RR_intervals(0, srLength); % sinus rhythm rr series (real)
            hrMean = (1./diff(rr_sr))*60 ;
        else % Use simulated RR series
            [rr_sr,hrMean,ecgParameters] = simECG_generate_sinus_rhythm(srLength,ecgParameters); % sinus rhythm rr series
            sigLength = ceil(ecgParameters.Duration);
%             rr_sr = rr_sr(1:srLength);
        end
        
        % Generate atrial fibrillation pacing activity
        afLength = rrLength;
        if realRRon == 1 % Use real RR series
            rr_af = simECG_get_real_RR_intervals(1, afLength);  %atrial fibrillation rr series
        else
            rr_af = simECG_generate_AF_intervals(fibFreqz,afLength); %atrial fibrillation rr series
            rr_af = rr_af(1:afLength);
        end
        
        % Make sure that the average RR value during SR is larger than that in AF
        if mean(rr_sr) < mean(rr_af)
            rr_sr = rr_sr + (mean(rr_af) - mean(rr_sr));
        end
        
        
    case 3 % The entire rhythm is SR but for EST signal %CPerez 06/2021
        srLength = [];
        if realRRon == 1 % Use real RR series
            [rr_sr, ecgParameters.peak, fa ecgParameters.ecgnr] = simECG_get_real_RR_intervals(2,rrLength); % If opt == 1 - AF, If opt == 2 - SR for stress test
             
            L = 300;
            if  realRRon %CPerez 04/2022. To take into account the influence of the "QT memory" before starting the EST
                rep = ceil(L/rr_sr(1));
                rr_sr = [repmat(rr_sr(1),rep,1); rr_sr];
            end
             
            dHR = 60./rr_sr; %instantaneous HR
            timebeats = [0; cumsum(rr_sr(2:end))]; %In seconds
            [pos, dHRu] = resamp(dHR,timebeats.*fa,fa,4);%4Hz resample
            pos = pos./fa;
            [b,a] = butter(2,0.03); %Fc = 0.03Hz
            hrMeanu = filtfilt(b,a,dHRu);
            hrMean = interp1(pos, hrMeanu, timebeats);
            ecgParameters.Duration = sum(rr_sr); %in seconds
            
            %Recalculate the exercise peak
            posPeak = sort(find(round(rr_sr,3)==min(round(rr_sr,3)))); %peak position
            if length(posPeak)>1
                posPeak = posPeak(round(length(posPeak)/2));
            end
            ecgParameters.peak = sum(rr_sr(1:posPeak));
            
            
            ecgParameters.Frini = [];
            ecgParameters.Frpeak = [];
            ecgParameters.Frend = [];
            
            sigLength = ceil(ecgParameters.Duration);
            sigLengthMs = sigLength * 1000;
            
        else % Use synthetic RR series CPerez 04/2022
            [rr_sr, hrMean, ecgParameters] = simECG_generate_sinus_rhythm(rrLength, ecgParameters);
            sigLength = fix(ecgParameters.Duration); %in sec
            sigLengthMs = sigLength * 1000; %in msec
        end
        rhythm_states = ones(1,sigLength);
        rr_af = Inf;
end

%% Atrial and ventricular ectopic pacing generation
[ap,vp] = simECG_generate_ectopic_pacing(sigLength,arrhythmiaParameters,hrMean);
% sinus pacing
sp = round(cumsum(rr_sr)*1000);
% AF pacing
afp = round(cumsum(rr_af)*1000);
% setting last element to Inf allows fast beat search by "min"
sp(end+1) = Inf; afp(end + 1) = Inf; ap(end + 1) = Inf; vp(end + 1) = Inf;

%% Construct global RR series
rr = Inf;
% global rr counter
k = 1;
% rhythm counters
c = [1, 1, 1, 1]; % c(1): sinus rhythm, c(2): af, c(3): apbs, c(4): vpbs
% annotation counter
kann = 1;
% time counter
t = 0;
% rhythm check - used to annotate rhythm changes
rc = 0;
% flags used to determine wheter an atrial or ventricular beat was inserted
ab_flag = 0; vb_flag = 0;
% rr series and targets initialization
targets_beats=0;
% convert rhythm states from seconds to ms
rhythm_states_s = zeros(1,sigLengthMs);
for n = 1:sigLength
    rhythm_states_s(1,1+(n-1)*1000:n*1000) = rhythm_states(1,n)*ones(1,1000);
end

while t<=sigLengthMs
    % discard beats occurred before t + 200 ms
    while sp(c(1))<t
        c(1) = c(1) + 1;
    end
    while afp(c(2))<t
        c(2) = c(2) + 1;
    end
    while ap(c(3))<t
        c(3) = c(3) + 1;
    end
    while vp(c(4))<t
        c(4) = c(4) + 1;
    end
    [next,idx] = min([sp(c(1)),ap(c(3)),vp(c(4))]);
    if next>sigLengthMs
        break;
    end
    % check whether next beat is in SR of AF
    if (rhythm_states_s(next)==1)
        % sinus rhythm
        if idx==2
            % insert APB
            % check if the beat is part of a run
            if ab_flag > 0
                % check if the run has not been annotated yet
                if rc~=3
                    % rhythm annotation - SVTA
                    % shuffle last apb annotation after the SVTA one
                    annTime(kann,1) = annTime(kann-1,1);
                    annBeats(kann,1) = annBeats(kann-1,1);
                    annRhythm{kann,1} = annRhythm{kann-1,1};
                    % actual SVTA annotation
                    annTime(kann-1,1) = t - round(1.5*(rr(k-1)*1000));
                    annBeats(kann-1,1) = beatCodes(end);
                    annRhythm{kann-1,1}= rhythmCodes(3);
                    kann=kann+1;
                    rc=3;
                end
                % current rr interval of the run
                rrtemp = ( ap(c(3)) - ap(c(3)-1) ) / 1000;
                % sinus rhythm rr series is delayed
                sp(c(1):end) = sp(c(1):end) + ( ap(c(3)) - ap(c(3)-1) );
            else
                % prematurity factor
                f1 = max(min(( 0.725 + (0.125.*randn) ),0.9),0.325);
                % rr interval of the apb
                rrtemp = rr_sr(c(1)) * f1; % premature intervals from the MIT BIH database
                % delay factor
                f2 = max(min(( 1.15 + (0.1.*randn) ),1.6),0.5);
                % apb compensatory pause
                rr_sr(c(1)) = rr_sr(c(1)) * f2;
                % sinus rhythm rr series is delayed
                sp(c(1):end) = sp(c(1):end) + t + round( (rrtemp + rr_sr(c(1)))* 1000 ) - sp(c(1));
            end
            % insert premature atrial beat
            rr(k) = rrtemp;
            % update targets array
            targets_beats(k) = 3;
            % update time counter
            t = t + round(rr(k)*1000);
            % update global rr counter
            k = k + 1;
            % ap series is adjusted
            ap(c(3):end) = ap(c(3):end) + t-ap(c(3));
            % beat annotation
            annTime(kann,1) = round(t);
            annBeats(kann,1) = beatCodes(2);
            annRhythm{kann,1} = [];
            kann=kann+1;
            % this was an atrial beat
            ab_flag = 1;
            % update atrial beats counter
            c(3) = c(3) + 1;
        else
            if idx==3
                % insert VPB
                % rr interval of the vpb
                rrtemp = rr_sr(c(1)) * ( 0.5 + 0.2*(rand) );
                % vpb compensatory pause
                rr(k) = rr_sr(c(1)) * ( 1.4 + 0.2*(rand) );
                % insert premature ventricular beat
                rr(k) = rrtemp;
                % update targets array
                targets_beats(k) = 4;
                % update time counter
                t = t + round(rr(k)*1000);
                % vp series is adjusted
                vp(c(4):end) = vp(c(4):end) + t-vp(c(4));
                % update global rr counter
                k = k + 1;
                % beat annotation
                annTime(kann,1) = round(t);
                annBeats(kann,1) = beatCodes(4);
                annRhythm{kann,1} = [];
                kann=kann+1;
                % this was a ventricular beat
                vb_flag = 1;
                % update ventricular beats counter
                c(4) = c(4) + 1;
            else
                % insert normal beat
                rr(k) = rr_sr(c(1));
                % update time counter
                t = t + round(rr(k)*1000);
                % sp series is adjusted
                sp(c(1):end) = sp(c(1):end) + t-sp(c(1));
                % update sinus rhythm counter
                c(1) = c(1) + 1;
                % check if rhythm annotation should change
                if rc~=1
                    % rhythm annotation - N
                    annTime(kann,1) = t + round((rr(k)*1000)/2);
                    annBeats(kann,1) = beatCodes(end);
                    annRhythm{kann,1}= rhythmCodes(1);
                    kann=kann+1;
                    % rhythm check
                    rc = 1;
                end
                % beat annotation
                annTime(kann,1) = round(t);
                annBeats(kann,1) = beatCodes(1);
                annRhythm{kann,1} = [];
                kann=kann+1;
                % update targets array
                targets_beats(k) = 1;
                % update global rr counter
                k = k + 1;
                % this was not an atrial or ventricular beat
                ab_flag = 0; vb_flag = 0;
            end
        end
    else
        % atrial fibrillation
        [next,idx] = min([afp(c(2)),vp(c(4))]); % no APBs in AF
        if next>sigLengthMs
            break;
        end
        if idx==2
            % insert VPB
            % rr interval of the vpb
            rrtemp = rr_sr(c(1)) * ( 0.5 + 0.2*(rand) );
            % insert premature ventricular beat
            rr(k) = rrtemp;
            % update targets array
            targets_beats(k) = 4;
            % update time counter
            t = t + round(rr(k)*1000);
            % vp series is adjusted
            vp(c(4):end) = vp(c(4):end) + t-vp(c(4));
            % update global rr counter
            k = k + 1;
            % beat annotation
            annTime(kann,1) = round(t);
            annBeats(kann,1) = beatCodes(4);
            annRhythm{kann,1} = [];
            kann=kann+1;
            % this was a ventricular beat
            vb_flag = 1;
            % update ventricular beats counter
            c(4) = c(4) + 1;
        else
            % insert af beat from Corino's AV model
            rr(k) = rr_af(c(2));
            % update time counter
            t = t + round(rr(k)*1000);
            % afp series is adjusted
            afp(c(2):end) = afp(c(2):end) + t-afp(c(2));
            % update af counter
            c(2) = c(2) + 1;
            % check if rhythm annotation should change
            if rc~=2
                % rhythm annotation - N
                annTime(kann,1) = t + round((rr(k)*1000)/2);
                annBeats(kann,1) = beatCodes(end);
                annRhythm{kann,1}= rhythmCodes(2);
                kann=kann+1;
                % rhythm check
                rc = 2;
            end
            % beat annotation
            annTime(kann,1) = round(t);
            annBeats(kann,1) = beatCodes(2);
            annRhythm{kann,1} = [];
            kann=kann+1;
            % update targets array
            targets_beats(k) = 2;
            % update global rr counter
            k = k+1;
            % this was not an atrial or ventricular beat
            ab_flag = 0; vb_flag = 0;
        end
    end
end

% Find the number of AF episodes and the length of each episode
k = 1;
pafEpisodeLength = [];
for p = 1:length(rhythm_states_s)-1
    if rhythm_states_s(p) == 2
        if rhythm_states_s(p+1) == 2
            k = k + 1;
            if p == (length(rhythm_states_s)-1)
                pafEpisodeLength = [pafEpisodeLength k];
            end
        else
            pafEpisodeLength = [pafEpisodeLength k];
            k = 1;
        end
    end
end
% Find boundaries of each PAF episode
if prod(rhythm_states_s == 1)==1 % The entire signal is SR
    pafBoundaries(1,1) = 0;
    pafBoundaries(1,2) = 0;
elseif prod(rhythm_states_s == 2)==1 % The entire signal is AF
    pafBoundaries(1,1) = 1;
    pafBoundaries(1,2) = length(rr);
else % The signal with PAF
    diffTar = diff(rhythm_states_s);
    j = 1;
    k = 1;
    flag = 1;
    for i = 1:length(diffTar)
        if diffTar(i) == 1
            pafBoundaries(j,1) = i + 1;
            j = j + 1;
            flag = 1;
        end
        if diffTar(i) == -1
            pafBoundaries(k,2) = i;
            k = k + 1;
            flag = 2;
        end
        
        if i == length(diffTar)
            if flag == 1
                pafBoundaries(k,2) = i+1;
            end
        end
    end
end
temp = 1; % debug

end