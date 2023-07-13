function [states] = simECG_state_gen(l2h, h2l)
% [] = simECG_state_gen() generates a markov chain with two states with 30
% samples, which are then upsampled by a factor of 1000 to match with a 30s
% thumb ECG with sampling frequency of 1000 Hz. Thresholding is done to
% remove the ringing that otherwise would have occured.
% 
% Input arguments:
% l2h - probability of transition to noisy segment.
% h2l - probability of transition to clean segment.
% 
% Output arguments:
% states - the states of the Markov chain.
% 
% Author: Hesam Halvaei, Lund University
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html
%%
N = 30;
P = ([1-l2h, l2h; h2l, 1-h2l]); % transition matrix
mc = dtmc(P);
states = simulate(mc, N);
states = interp(states, 1000);

states(states > 1.5) = 2;
states(states < 1.5) = 1;

states(states == 1) = 0;
states(states == 2) = 1;
end

