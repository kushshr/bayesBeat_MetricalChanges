clear
close all
clc
% 1. R = 1, combine *all* results for a method from all datasets
addpath('/media/Code/UPFWork/PhD/githubCode/CommonPool');
addpath('/media/Code/UPFWork/PhD/CommonPoolCodeGeneral/export_fig/');
basepath = '/media/Code/UPFWork/PhD/';
nExp = 3;
amlRatios = [0.25 0.5 1 2 4];
numPatts = [1];
taskpaths = {'BayesResultsThesis', 'BayesResultsThesisServer'};
measures = {'bFmeas', 'bAMLt', 'bInfoGain', 'sFmeas', 'tempoCML', 'tempoAML'};
% Some specific things here
%%%%%%%% Set of figures-1
datasets = {'CMCMDa_small', 'HMDs', 'HMDl'};
fExt = {'.txt', '.txt', '.p7.txt'};
talaLists = {{'adi', 'rupaka', 'mChapu', 'kChapu'}, ...
    {'teen', 'ek', 'jhap', 'rupak'}, ...
    {'teen', 'ek', 'jhap', 'rupak'}};
methods = {{'Inference', 'HMM_prior_allHop_bar'}, ...
    {'Inference' 'AMPF_prior_allHop_bar'}, ...
    {'Tracking' 'HMM_prior_allHop_bar'}, ...
    {'Tracking' 'AMPF_prior_allHop_bar'}, ...
    {'Tracking' 'AMPF_mix_allHop_bar'}, ...
    {'Tracking' 'AMPF_prior_allHop_section'}, ...
    {'tempoInfTracking' 'AMPF_prior_allHop_bar'}, ...
    {'tempoInfTracking' 'AMPF_prior_allHop_section'}, ...
    {'tempoSamaInfTracking' 'AMPF_prior_allHop_bar'}, ...
    {'tempoSamaInfTracking' 'AMPF_prior_allHop_section'}};
methodNames = {'InfHMMo', 'InfAMPFo', 'TrackHMMo', 'TrackAMPFo', ...
    'TrackAMPFm', 'TrackAMPFs', 'tempoInfTrackAMPFo', 'tempoInfTrackAMPFs', ...
    'tempoSamaInfTrackAMPFo', 'tempoSamaInfTrackAMPFs'}; 
figStub = '-statTest';
fileStub = 'statTest-all-IndianMusic';
%%%%%% Set of figures-2
% datasets = {'CMCMDa_small', 'HMDs'};
% fExt = {'.txt', '.txt'};
% talaLists = {{'adi', 'rupaka', 'mChapu', 'kChapu'}, ...
%     {'teen', 'ek', 'jhap', 'rupak'}};
% methods = {{'Tracking' 'AMPF_prior_allHop_bar'}, ...
%     {'Tracking' 'AMPF_acc_allHop_bar'}, ...
%     {'Tracking' 'AMPF_prior_peakHop_bar'}, ...
%     {'Tracking' 'AMPF_prior_obsHop_bar'}};
% methodNames = {'TrackAMPFo', 'TrackAMPFe', 'TrackAMPFp', 'TrackAMPFg'}; 
% figStub = '-statTest-infExt';
% fileStub = 'statTest-infExt-IndianMusic';
for m = 1:length(methods)
    % Accumulate values for each dataset together
    resAll = [];
    resNamesAll = {};
    for d = 1:length(datasets)
        resFull = [];
        resNames = {};
        for t = 1:length(talaLists{d})
            if strcmp(methods{m}{1}, 'Inference')
                fstub = fullfile(datasets{d}, methods{m}{1}, methods{m}{2}, 'nPatts_1');
            else
                fstub = fullfile(datasets{d}, methods{m}{1}, methods{m}{2}, talaLists{d}{t}, 'nPatts_1');
            end
            if exist(fullfile(basepath, taskpaths{1}, fstub),'dir')
                fbasepath = fullfile(basepath, taskpaths{1}, fstub);
            elseif exist(fullfile(basepath, taskpaths{2}, fstub),'dir')
                fbasepath = fullfile(basepath, taskpaths{2}, fstub);
            else
                error('Data not found!!!');
            end
            for ex = 1:nExp
                fname = [fbasepath filesep methods{m}{2} '_' talaLists{d}{t} ...
                    '_nPatts_1_' num2str(ex) '_allResults' fExt{d}];
                tbl = readtable(fname,'delimiter',',');
                tempoCML = double(abs(tbl.EstTempo - tbl.MedianTempo)./(tbl.MedianTempo) < 0.05);
                for jj = 1:length(tbl.EstTempo)
                    tempoAML(jj,1) = double((sum((abs(tbl.EstTempo(jj) - tbl.MedianTempo(jj).*amlRatios)./(tbl.MedianTempo(jj).*amlRatios)) < 0.05)) > 0);
                end
                res(:,1,ex) = tbl.bfMeas/100;
                res(:,2,ex) = tbl.bAMLt/100;
                res(:,3,ex) = tbl.bInfoGain;
                res(:,4,ex) = tbl.sfMeas/100;
                res(:,5,ex) = tempoCML;
                res(:,6,ex) = tempoAML;
                fnames = tbl.File;
                clear tempoCML tempoAML tbl
            end
            resFull = [resFull; res];
            resNames = [resNames; fnames];
            clear res fnames
        end
        resAll = [resAll; resFull];
        resNamesAll = [resNamesAll; resNames];
        clear resFull resNames
    end
    resOp(m).fullRes = resAll;
    resOp(m).meanRes = nanmean(resAll,3);
    resOp(m).stdRes = std(resAll,1,3);
    resOp(m).fnames = resNamesAll;
    resOp(m).measures = measures;
    resOp(m).methodName = methodNames{m};
    resOp(m).meanmeanRes = nanmean(resOp(m).meanRes);
    clear resAll resNamesAll
    m
end
fp = fopen([fileStub '.txt'], 'wt');
fprintf(fp, 'Algorithm && \\fmeas & \\amlt & \\infoGain && $\\fmeas_s$ && \\multicolumn{2}{c}{Tempo} \\tabularnewline \n');
fprintf(fp, '&& & & Bits && && CML & AML \\tabularnewline \\midrule \n');
for m = 1:length(methods)
    fprintf(fp, '%s && %.3f & %.3f & %.2f && %.3f && %.2f & %.2f \\tabularnewline \n', ...
       resOp(m).methodName, resOp(m).meanmeanRes); 
end
fclose(fp); 
for mm = 1:length(measures)
    statTest(mm).meas = measures{mm};
    for k1 = 1:length(methods)
        for k2 = 1:length(methods)
            [statTest(mm).hVal(k1,k2), statTest(mm).pVal(k1,k2)] = ttest(resOp(k1).meanRes(:,mm),resOp(k2).meanRes(:,mm));
            if ~isnan(statTest(mm).hVal(k1,k2)) && statTest(mm).hVal(k1,k2) && mm <= 4
                fprintf('Significantly different:- Measure: %s, Methods: %s and %s with p value = %.6f\n', ...
                    statTest(mm).meas, methodNames{k1}, methodNames{k2}, statTest(mm).pVal(k1,k2))
            end
        end
    end
    statTest(mm).methodName = methodNames;
end
%% Now to generate plots for each case
txoff = -0.1;
tyoff = -0.045;
for mm = 1:length(measures)
    hFig = figure('Position', [100,100,60*length(methods),50*length(methods)]);
    imagesc(statTest(mm).hVal)
    colormap([256 256 256; 220 220 220]/256);
    hold on, 
    for k = 1:(length(methods)-1)
        line([0.5+k,0.5+k],[0.5,length(methods)+0.5],'LineStyle',':','Color',[0 0 0]);
        line([0.5,length(methods)+0.5],[0.5+k,0.5+k],'LineStyle',':','Color',[0 0 0]);
    end
    xtxtPos = (1:length(methods)) + txoff;
    ytxtPos = (1:length(methods)) + tyoff;
    txtCoord = [];
    txtText = {};
    for k1 = 1:length(xtxtPos)
        for k2 = 1:length(ytxtPos)
            txtCoord = [txtCoord; [xtxtPos(k1) ytxtPos(k2)]];
            if ~isnan(statTest(mm).hVal(k1,k2))
                txtText = [txtText; num2str(statTest(mm).hVal(k1,k2))];
            else
                txtText = [txtText; '-'];
            end
        end
    end
    text(txtCoord(:,1),txtCoord(:,2),txtText,'FontSize',5,'FontName','Verdana');
    set(gca,'xaxisLocation','top')
    sparams.dimFig = [5,5];
    sparams.bkColor = 'None';
    set(gca, 'Ticklength', [0 0])
    saveThesisFigure(hFig, [statTest(mm).meas figStub], sparams, {'fig', '-png'});
    close all
    % set(gca,'YTickLabel',methodNames)
    % fSize = 12; 
    % fType = 'Times New Roman';
    % hText = getTextHandles(hFig); % All text handles
    % set(hText,'FontSize',fSize,'FontName',fType); % All text converted
end