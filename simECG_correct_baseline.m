function pqrst = simECG_correct_baseline(pqrstIn)
pqrst(1,:) = pqrstIn;
x = [1 length(pqrst)];
xi = 1:1:length(pqrst);
y = [pqrst(1,1) pqrst(1,end)];
baseline = interp1(x,y,xi,'linear');
pqrst = pqrst - baseline;
end