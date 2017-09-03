clear
close all
clc
basepath = '/media/Code/UPFWork/PhD/BayesResultsFull';
dataset = 'HMDs';
% dataset = 'BallroomDataset';
exptType = 'Tracking';
methods = {'PF_bar_AMPF_PriorTrans_NoHop_1500pp', ...
    'PF_section_AMPF_PriorTrans_NoHop_1500pp'};
methName = {'PF_bar', 'PF_sec'};
nTrials = [3 3];

% talas = {'adi', 'rupaka', 'mChapu', 'kChapu'};
talas = {'teen', 'ek', 'jhap', 'rupak'};
nPatts = [1];
perfMeas = {'sfMeas', 'bfMeas'};

% First read all the results
for m = 1:length(methods)
    results(m).name = methods{m};
    for p = 1:length(nPatts)
        results(m).patts(p).R = nPatts(p);
        for t = 1:length(talas)
            results(m).patts(p).tala(t).name = talas{t};
            load(fullfile(basepath,dataset,exptType,methods{m},...
                    talas{t}, ['nPatts_' num2str(nPatts(p))], 'Parameters.mat'));
            % results(m).patts(p).tala(t).medRunTime = median(Params.runTime{end,end});
            clear Params
            for nt = 1:nTrials(m)
                results(m).patts(p).tala(t).trial(nt).id = nt;
                fpath = fullfile(basepath,dataset,exptType,methods{m},...
                    talas{t}, ['nPatts_' num2str(nPatts(p))],...
                    [methods{m} '_' talas{t} '_nPatts_' num2str(nPatts(p)) '_'...
                    num2str(nt) '_allResults.txt']);
                res = readresults(fpath);
                results(m).patts(p).tala(t).trial(nt).res = res';
                clear res;
            end
        end
        % Now pool all tala results together
        allTalaInd = length(talas)+1;
        results(m).patts(p).tala(allTalaInd).name = 'All';
        for nt = 1:nTrials(m)
            results(m).patts(p).tala(allTalaInd).trial(nt).id = nt;
            results(m).patts(p).tala(allTalaInd).trial(nt).res = [];
            for t = 1:length(talas)
                results(m).patts(p).tala(allTalaInd).trial(nt).res = ...
                    [results(m).patts(p).tala(allTalaInd).trial(nt).res; ...
                    results(m).patts(p).tala(t).trial(nt).res];  
            end
        end
    end
end
%% And now for something completely different!
talas = [talas 'All'];  % First add "All" as a tala
% Intra trial significance testing
% Write to a report
fpresult = fopen(fullfile(basepath, dataset, exptType, [dataset '_' exptType '_statTestingResults_Trials.txt']),'wt');
fpreport = fopen(fullfile(basepath, dataset, exptType, [dataset '_' exptType '_statTestingReport_Trials.txt']),'wt');
for m = 1:length(methods)
    for p = 1:length(nPatts)
        for t = 1:length(talas)
            for k = 1:length(perfMeas)
                results(m).patts(p).tala(t).sigTestTrial(k).measure = perfMeas{k};
                results(m).patts(p).tala(t).sigTestTrial(k).pMat = ...
                    NaN(length(results(m).patts(p).tala(t).trial));
                results(m).patts(p).tala(t).sigTestTrial(k).hMat = ...
                    NaN(length(results(m).patts(p).tala(t).trial));
                accumRes = [];
                for a = 1:length(results(m).patts(p).tala(t).trial)
                    for b = 1:length(results(m).patts(p).tala(t).trial)
                        aVec = [results(m).patts(p).tala(t).trial(a).res.(perfMeas{k})];
                        bVec = [results(m).patts(p).tala(t).trial(b).res.(perfMeas{k})];
                        % [pVal hVal] = ranksum(aVec,bVec);
                        % [pVal hVal] = signrank(aVec,bVec);
                        [hVal pVal] = ttest(aVec,bVec);
                        results(m).patts(p).tala(t).sigTestTrial(k).pMat(a,b) = pVal;
                        results(m).patts(p).tala(t).sigTestTrial(k).hMat(a,b) = hVal;
                        fprintf(fpreport, 'Stat test: Method-%s,Patts-%d,Tala-%s,Measure-%s. Between trials: %d and %d, p-value = %.4f\n',...
                                methName{m}, nPatts(p), talas{t}, perfMeas{k}, a, b, pVal);
                        if ~isnan(hVal) && hVal
                            fprintf(fpresult, 'Significantly different: Method-%s,Patts-%d,Tala-%s,Measure-%s. p-value = %.4f for test between Trials %d and %d\n',...
                                methName{m}, nPatts(p), talas{t}, perfMeas{k}, pVal, a, b)
                            fprintf('Significantly different: Method-%s,Patts-%d,Tala-%s,Measure-%s. p-value = %.4f for test between Trials %d and %d\n',...
                                methName{m}, nPatts(p), talas{t}, perfMeas{k}, pVal, a, b)
                        end
                        clear pVal hVal aVec bVec
                    end
                    accumRes = [accumRes [results(m).patts(p).tala(t).trial(a).res.(perfMeas{k})]'];
                end
                results(m).patts(p).tala(t).meanRes.(perfMeas{k}) = nanmean(accumRes,2);
            end
        end
    end
end
fclose(fpreport);
fclose(fpresult);
% Inter method significance test, keeping patterns fixed
fpresult = fopen(fullfile(basepath, dataset, exptType, [dataset '_' exptType '_statTestingResults_Methods.txt']),'wt');
fpreport = fopen(fullfile(basepath, dataset, exptType, [dataset '_' exptType '_statTestingReport_Methods.txt']),'wt');
for p = 1:length(nPatts)
    for t = 1:length(talas)
        for k = 1:length(perfMeas)
            sigTestMethod(p).perfMeasure(k).name = perfMeas{k};
            sigTestMethod(p).perfMeasure(k).tala(t).pMat = NaN(length(methods));
            sigTestMethod(p).perfMeasure(k).tala(t).hMat = NaN(length(methods));
            for m1 = 1:length(methods)
                for m2 = 1:length(methods)
                    aVec = [results(m1).patts(p).tala(t).meanRes.(perfMeas{k})];
                    bVec = [results(m2).patts(p).tala(t).meanRes.(perfMeas{k})];
                    % [pVal hVal] = ranksum(aVec,bVec);
                    % [pVal hVal] = signrank(aVec,bVec);
                    [hVal pVal] = ttest(aVec,bVec);
                    sigTestMethod(m).perfMeasure(k).tala(t).pMat(m1,m2) = pVal;
                    sigTestMethod(m).perfMeasure(k).tala(t).hMat(m1,m2) = hVal;
                    fprintf(fpreport, 'Stat test across methods: Pattern-%d,Tala-%s,Measure-%s. Between methods: %s and %s, p-value = %.4f\n',...
                                nPatts(p), talas{t}, perfMeas{k}, methName{m1}, methName{m2}, pVal);
                    if ~isnan(hVal) && hVal
                        fprintf(fpresult, 'Significantly different: Pattern-%d,Tala-%s,Measure-%s. p-value = %.4f for test between Methods %s and %s\n',...
                                nPatts(p), talas{t}, perfMeas{k}, pVal, methName{m1}, methName{m2})
                        fprintf('Significantly different: Pattern-%d,Tala-%s,Measure-%s. p-value = %.4f for test between Methods %s and %s\n',...
                                nPatts(p), talas{t}, perfMeas{k}, pVal, methName{m1}, methName{m2})
                    end
                    clear pVal hVal aVec bVec
                end
            end
        end
    end
end
fclose(fpreport);
fclose(fpresult);