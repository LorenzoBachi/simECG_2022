function pqrst = simECG_correct_baseline(pqrstIn)
% [] = simECG_correct_baseline() returns the baseline-corrected PQRST
% complex.
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

pqrst(1,:) = pqrstIn;
x = [1 length(pqrst)];
xi = 1:1:length(pqrst);
y = [pqrst(1,1) pqrst(1,end)];
baseline = interp1(x,y,xi,'linear');
pqrst = pqrst - baseline;
end