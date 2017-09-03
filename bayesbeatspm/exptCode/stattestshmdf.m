st = [1 14 46 52 1; 1 42 68 81 1];
en = [13 45 51 59 59; 41 67 80 92 92];

bval = {'bar', 'sec'};
dval = {'long', 'short'};

for d = 1:2
    bardata = dlmread(['bar_', dval{d}, '.txt']);
    secdata = dlmread(['sec_', dval{d} '.txt']);
    bardata = bardata(:,[1 4]);
    secdata = secdata(:,[1 4]);
    for t = 1:5
        [shVal(d,t) spVal(d,t)] = ttest(bardata(st(d,t):en(d,t),1),secdata(st(d,t):en(d,t),1));
        [mhVal(d,t) mpVal(d,t)] = ttest(bardata(st(d,t):en(d,t),2),secdata(st(d,t):en(d,t),2));
    end
end
