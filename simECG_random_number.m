function number = simECG_random_number(rangeLow, rangeHigh)
% number = simECG_random_number() returns a random number from the interval
% [rangeLow, rangeHigh].
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

number = (rangeHigh-rangeLow)*rand(1,1) + rangeLow;
number = round(number*100)/100;

end