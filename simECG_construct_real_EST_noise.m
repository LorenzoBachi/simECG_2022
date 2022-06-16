function noise_15=simECG_construct_real_EST_noise(noise_9)

%Function
% construct_EST_noise
% 
%Purpose
% obtain 15 leads of a noise signal extracted of an Exercise Stress Test recording.
% (R.Bailón)
%
%Input
% noise_9: 9 leads of stress test noise (V1, RV4, V3-V6, I, II, III), 
%whose sample frequency is 1000 Hz.
%
%Output:
% noise_15: 15 leads of stress test noise (I, II, III, aVR, aVL, aVF, V1-V6, XYZ)
%
% Copyright (c), Cristina Perez, University of Zaragoza, 06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1)Obtain augmented unipolar limb leads
noise_12 = leadcalc(noise_9,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III

%2) Obtain synthesized VCG (leads X, Y and Z)
noise_3 = leadcalc(noise_9,'hcuzst'); %XYZ

%3) Final 15-lead noise
noise_15=vertcat(noise_12(8,:),noise_12(10,:),noise_12(12,:),...
    -noise_12(9,:),noise_12(7,:),noise_12(11,:),noise_12(1:6,:),noise_3);
end

