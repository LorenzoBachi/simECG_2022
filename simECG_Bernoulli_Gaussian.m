function [bernogauss] = simECG_Bernoulli_Gaussian(n, prob, sigma)
% [] = simECG_Bernoulli_Gaussian() returns the realization of a Bernoulli-Gaussian
% process.
% 
% Input arguments:
% n - length of signal, in samples.
% prob - probability of success in the Bernoulli Gaussian Process.
% sigma - noise level, standard deviation of guassian noise.
% 
% Output arguments:
% bernogauss - Bernoulli-Gaussian process realization.
% 
% Author: Hesam Halvaei, Lund University
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

N = ones(n, 1);
p = N.*prob;
berno = binornd(N, p);
gauss= normrnd(0, sigma, [n, 1]);
bernogauss = gauss.*berno;
end
