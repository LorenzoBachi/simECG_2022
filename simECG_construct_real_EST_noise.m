function [noise_15] = simECG_construct_real_EST_noise(noise_9)
% [] = simECG_construct_real_EST_noise() returns a 15 lead matrix of noise
% extracted of an exercise stress test recording (R.Bailón).
% 
% Input arguments:
% noise_9 - 9 leads of stress test noise (V1, RV4, V3-V6, I, II, III), 
% sampled at 1000 Hz.
% 
% Output arguments:
% noise_15 - 15 leads of stress test noise (I, II, III, aVR, aVL, aVF,
% V1-V6, XYZ).
% 
% Author: Cristina Perez, University of Zaragoza
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

%1)Obtain augmented unipolar limb leads
noise_12 = leadcalc(noise_9,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III

%2) Obtain synthesized VCG (leads X, Y and Z)
noise_3 = leadcalc(noise_9,'hcuzst'); %XYZ

%3) Final 15-lead noise
noise_15=vertcat(noise_12(8,:),noise_12(10,:),noise_12(12,:),...
    -noise_12(9,:),noise_12(7,:),noise_12(11,:),noise_12(1:6,:),noise_3);
end

