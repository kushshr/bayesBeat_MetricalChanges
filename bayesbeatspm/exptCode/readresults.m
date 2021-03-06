function res = readresults(fpath)
% colHead1 = 'File,Tala,MedianTempo,EstMeter,EstRhythm,EstTempo,sfMeas,sPrec,sRecall,';
% colHead2 = 'sCMLt,sAMLt,sInfoGain,bfMeas,bPrec,bRecall,bCMLt,bAMLt,bInfoGain';
fp = fopen(fpath, 'rt');
header = fgetl(fp);
A = textscan(fp, '%s%s%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
fclose(fp);
for k = 1:length(A{1})
    res(k).fname = A{1}{k};
    res(k).talaName = A{2}{k};
    res(k).bpmGT = A{3}(k);
    res(k).meter = A{4}{k};
    res(k).rhythm = A{5}{k};
    res(k).bpm = A{6}(k);
    res(k).sfMeas = A{7}(k)/100;
    res(k).sPrecision = A{8}(k)/100;
    res(k).sRecall = A{9}(k)/100;
    res(k).sCmlt = A{10}(k);
    res(k).sAmlt = A{11}(k);
    res(k).sInfoGain = A{12}(k);
    res(k).bfMeas = A{13}(k)/100;
    res(k).bPrecision = A{14}(k)/100;
    res(k).bRecall = A{15}(k)/100;
    res(k).bCmlt = A{16}(k);
    res(k).bAmlt = A{17}(k);
    res(k).bInfoGain = A{18}(k);
    %res(k).fID = (str2double(res(k).fname(4:5)) - 19)*30 + str2double(res(k).fname(1:2));
    res(k).talaID = (str2double(res(k).fname(4:5)) - 19); % CAREFUL!!
end
