function number = simECG_random_number(rangeLow, rangeHigh)
%
% number = simECG_random_number() returns a random number from the interval
% [rangeLow, rangeHigh].
%
number = (rangeHigh-rangeLow)*rand(1,1) + rangeLow;
number = round(number*100)/100;

end