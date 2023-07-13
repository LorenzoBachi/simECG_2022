function fibFreqz = simECG_fibrillation_frequency()
% fibFreqz = simPAF_gen_initial_parameters returns initial parameters (at
% this point only a dominant fibrillatory frequency). The dominant
% frequency of f-waves is PAF-specific, from the interval [3 7] Hz. See
% a paper by Petrutiu et al. Abrupt changes in fibrillatory wave
% characteristics of paroxysmal atrial fibrillation in humans. Europace, 9,
% 466-470. 2007. 
%
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html 

fibFreqz = randi([30 70])/10;
end