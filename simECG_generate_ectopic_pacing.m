function [aa,va] = simECG_generate_ectopic_pacing(sigLength,arrhythmiaParameters,hrMean)
% [] = simECG_generate ectopic activity() generates the pacing activity
% from ectopic foci in the atria and the ventricles.

% please insert additional comments here

% initializing atrial and ventricular ectopic event series
aa = Inf;
va = Inf;
% number of ectopic events
naa=0;
nva=0;
% refractory period
rp = 200;

%% Atrial ectopic pacing
% atrial ectopic pacing event rate - from APBs per hour to APBs per s
lambdaA = arrhythmiaParameters.APBph / ( 60 * 60 );
% atrial ectopic event times, in [seconds]
aee = exprnd(1/lambdaA,sigLength,1);
% atrial ectopic events in the selected time frame
aee = cumsum(aee(cumsum(aee)<sigLength));
% atrial ectopic event times, in [ms] and, thus, samples
aee = round(aee.*1000);
% number of events
N = length(aee);
% atrial run rate
ARRate = arrhythmiaParameters.ARRate;
% atrial run creation
% fetching porbability distribution of atrial run length
ARDist = arrhythmiaParameters.ARDist;
ARDist =ARDist./sum(ARDist);
cumADist = cumsum(ARDist);
% range of possible ( AR BPM ) / ( SR BPM ) ratio
arBpmRange = [1.5, 2.2];
%maximum index for activity
lengthMs = sigLength*1000;
for k=1:N
    % new atrial ectopic event
    naa = naa + 1;
    % timestap of the event
    aa(naa) = aee(k);
    % chance of an atrial run
    chance = rand;
    if chance<ARRate
        % an atrial run begins.
        % pick a duration: sampling cumulative distribution
        x = find(rand<=cumADist);
        % duration is the first nonzero index of rand<=cumDist +1 since the
        % first element of ARDist refers to couplets of APBs
        d = x(1)+1;
        % pick a BPM frequency
        bpm = (rand*(arBpmRange(2)-arBpmRange(1)) + arBpmRange(1)) * hrMean;
        % turn BPM in RR interval, in [ms]
        rrMean = (60/bpm)*1000;
        % rr interval standard devation
        rrStd = rrMean/20;
        % counter for run insertion
        l=2;
        while l<d
            % rr distance from previous atrial beat
            rri = round( rrMean + round(rrStd*randn) );
            % position of the next atrial beat in run
            temp = aa(naa) + rri ;
            %check if the run is not over the max duration of the simulated
            %record
            if aa(naa)<lengthMs
                % new atrial event
                naa=naa+1;
                % timestap of the event
                aa(naa) = temp;
                l=l+1;
            else
                %quit loop
                break;
            end
        end
    end
end
% move ectopic events if they are too close to each other (refractory
% period)
for k=2:naa
    if (aa(k)-aa(k-1)) < rp
        aa(k)=aa(k-1) + rp - (aa(k)-aa(k-1)) + round(rp*rand*0.25);
    end
end


%% Ventricular ectopic pacing
% atrial ectopic pacing event rate
lambdaV = arrhythmiaParameters.VPBph / ( 60 * 60 );
% ventricular ectopic event times, in [seconds]
vee = exprnd(1/lambdaV,sigLength,1);
% ventricular ectopic events in the selected time frame
vee = cumsum(vee(cumsum(vee)<sigLength));
% ventricular ectopic event times, in [ms] and, thus, samples
vee = round(vee.*1000);
% number of events
N = length(vee);
for k=1:N
    % new ventricular ectopic event
    nva = nva + 1;
    % timestap of the event
    va(nva) = vee(k);
end
% move ectopic events if they are too close to each other (refractory
% period)
for k=2:nva
    if (va(k)-va(k-1)) < rp
        va(k)=va(k-1) + rp - (aa(k)-aa(k-1)) + round(rp*rand*0.25);
    end
end




end

