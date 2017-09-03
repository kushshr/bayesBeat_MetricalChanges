clear
close all
clc
addpath('../../CommonPoolCodeGeneral/Davies_beat_error_histogram/');
basepath = '/media/Code/UPFWork/PhD/';
nExp = 3;
nFold = 2;
amlRatios = [0.25 0.5 1 2 4];
numPatts = [1 2];
dataset = 'CMCMDa_small';
% dataset = 'BallroomDataset';
% dataset = 'HMDs';
% dataset = 'HMDf';
% dataset = 'HMDl';
% dataset = 'CMCMDa_v2';
% dataset = 'TurkishMakam';
% dataset = 'CretanLeapingDances';
% task = 'Inference'; taskStub = 'inference'; taskID = 'Inf';
task = 'Tracking'; taskStub = 'tracking'; taskID = 'Tr';
% task = 'tempoInfTracking'; taskStub = 'tempo-informed tracking'; taskID = 'Itr';
% task = 'tempoSamaInfTracking'; taskStub = 'tempo-sama-informed tracking'; taskID = 'Istr';
% task = 'samaInfTracking'; taskStub = 'sama-informed tracking'; taskID = 'Isr';
bpath = [fullfile(basepath, 'BayesResultsNpartSPM', dataset, task) filesep];
dpath = [fullfile(basepath, 'Data', dataset) filesep];
annpath = [dpath 'annotations/beats/'];
meterpath = [dpath 'annotations/meter/'];
exptName = 'AMPF_prior_allHop_bar';
exppath = [bpath exptName filesep];
if strcmp(dataset, 'CMCMDa_small') ||strcmp(dataset, 'CMCMDa_v2')
    talaID = 10:13;
    talaName = {'adi', 'rupaka', 'mChapu', 'kChapu'};
    tableOpName = {'\Gls{adi}', '\Gls{rupaka}', '\Gls{mishra chapu}', '\Gls{khanda chapu}'};
    if strcmp(dataset, 'CMCMDa_small')
        dStub = '\acrshort{CMDs}';
        dName = 'CMDs';
    else
        dStub = '\acrshort{CMDf}';
        dName = 'CMDf';
    end
elseif strcmp(dataset, 'HMDf') || strcmp(dataset, 'HMDl') || strcmp(dataset, 'HMDs')
    talaID = 20:23;
    talaName = {'teen', 'ek', 'jhap', 'rupak'};
    tableOpName = {'\Gls{teental}', '\Gls{ektal}', '\Gls{jhaptal}', '\Gls{rupak}'};
    dStub = ['\acrshort{HMD' dataset(4) '}'];
    dName = ['HMD' dataset(4)];
elseif strcmp(dataset, 'BallroomDataset')
    talaID = 11:18;
    talaName = {'ChaChaCha', 'Jive' , 'Quickstep', 'Rumba' , 'Samba' , 'Tango', 'VienneseWaltz', 'Waltz'};
    tableOpName = {'Cha-Cha-Cha', 'Jive' , 'Quickstep', 'Rumba' , 'Samba' , 'Tango', 'Viennese waltz', 'Waltz'};
    dStub = 'Ballroom';
    dName = dStub;
elseif strcmp(dataset, 'TurkishMakam')
    talaID = 30:32;
    talaName = {'aksak', 'curcuna', 'duyek'};
    tableOpName = talaName;
    dStub = 'Turkish';
    dName = dStub;
end
prms = be_params;
for t = 1:length(talaName)
    fpp = fopen([dpath 'filelist_' talaName{t} '.txt']);
    filelist = textscan(fpp,'%s');
    filelist = filelist{1};
    fclose(fpp);
    for p = 1:length(numPatts)
        if strcmp(task,'Inference')
            talabasepath = [exppath 'nPatts_' num2str(numPatts(p)) filesep];
        else
            talabasepath = [exppath talaName{t} filesep 'nPatts_' num2str(numPatts(p)) filesep];
        end
        for ex = 1:nExp
            % Process each file now
            for k = 1:length(filelist)
                % Read the ground truth first
                annFile = [filelist{k} '.beats'];
                talaIndicator(k) = find(ismember(talaID, str2num(annFile(4:5))));
                fileTalaName{k} = talaName{talaIndicator(k)};
                fileID(k) = str2num(annFile(4:8));
                fileNum(k) = str2num(annFile(1:2));
                ann = load([annpath annFile],'-ascii');
                % If sama informed tracking, then consider only annotations starting first sama
                if strcmp(task, 'tempoSamaInfTracking')  
                    firstInd = find(ann(:,2) == 1,1);
                    ann = ann(firstInd:end,:);
                end
                samas = ann(ann(:,2)==1,1);
                beats = ann(:,1);
                % Meter GT
                fpm = fopen([meterpath filelist{k} '.meter'],'rt');
                temp = textscan(fpm,'%s');
                res(k).meterGT = temp{1}{1};
                fclose(fpm);
                % Tempo GT
                temp = sort(diff(beats));
                lenn = length(beats);
                res(k).bpmGT = 60./median(temp(round(lenn/10):round(0.9*lenn)));
                % Get the output file name next
                gotit = 0;
                for f = 1:nFold
                    sim_id = 1000*ex+f;
                    oFileName = [talabasepath num2str(sim_id) filesep filelist{k} '.beats.txt'];
                    if exist(oFileName, 'file')
                        oFiles{k} = oFileName;
                        gotit = 1;
                        break;
                    end
                end
%                 annout = dlmread(oFiles{k});
%                 annout(:,2) = round(10*(annout(:,2) - floor(annout(:,2))));
                if gotit
                    fp1 = fopen(oFiles{k},'rt');
                    anntemp = textscan(fp1,'%s %s\n');
                    fclose(fp1);
                    annout(:,1) = str2double(anntemp{1});
                    annout(:,2) = str2double(anntemp{2});
                    opSamas = annout((annout(:,2) == 1),1);
                    opBeats = annout(:,1);
                    % Read tempo
                    % res(k).bpm = dlmread(fullfile(filepath, [fname{k} '.bpm']));
                    % Read output meter
                    fpm = fopen([talabasepath num2str(sim_id) filesep filelist{k} '.meter.txt'],'rt');
                    temp = textscan(fpm,'%s');
                    res(k).meter = temp{1}{1};
                    fclose(fpm);
                    res(k).meterCorr = double(strcmp(res(k).meter, res(k).meterGT));
                    % Read output rhythm/s
                    fpr = fopen([talabasepath num2str(sim_id) filesep filelist{k} '.rhythm.txt'],'rt');
                    temp = textscan(fpr,'%s');
                    res(k).rhythm = [temp{1}{:}];
                    fclose(fpr);
                    % Tempo Estimated
                    temp2 = sort(diff(opBeats));
                    lenn2 = length(opBeats);
                    res(k).bpm = 60./median(temp2(round(lenn2/10):round(0.9*lenn2)));
                    res(k).tempoCML = double(((abs(res(k).bpmGT - res(k).bpm)/res(k).bpmGT) < 0.05));
                    res(k).tempoAML = double((sum((abs(res(k).bpm - res(k).bpmGT.*amlRatios) ...
                        ./(res(k).bpmGT.*amlRatios)) < 0.05)) > 0);                    
                    % Sama metrics
                    res(k).sama.pScore = be_pScore(samas,opSamas);
                    [res(k).sama.fMeas res(k).sama.precision res(k).sama.recall res(k).sama.Ameas]...
                        = be_fMeasure(samas,opSamas);
                    if length(opSamas) == 1   % hack for long cycles
                        if sum(abs(samas-opSamas) < prms.fMeasure.thresh)
                            res(k).sama.infoGain = log2(prms.informationGain.numBins);
                            res(k).sama.cmlC = 50;
                            res(k).sama.cmlT = 50;
                            res(k).sama.amlC = 50;
                            res(k).sama.amlT = 50;
                            disp('Only one sama found!!')
                        else
                            res(k).sama.infoGain = 0;
                            res(k).sama.cmlC = 0;
                            res(k).sama.cmlT = 0;
                            res(k).sama.amlC = 0;
                            res(k).sama.amlT = 0;
                        end
                    else
                        res(k).sama.infoGain = be_informationGain(samas,opSamas);
                        [res(k).sama.cmlC res(k).sama.cmlT res(k).sama.amlC res(k).sama.amlT]...
                            = be_continuityBased(samas,opSamas);
                    end
                    
                    % Beat metrics
                    res(k).beat.pScore = be_pScore(beats,opBeats);
                    [res(k).beat.fMeas res(k).beat.precision res(k).beat.recall res(k).beat.Ameas]...
                        = be_fMeasure(beats,opBeats);
                    res(k).beat.infoGain = be_informationGain(beats,opBeats);
                    [res(k).beat.cmlC res(k).beat.cmlT res(k).beat.amlC res(k).beat.amlT]...
                        = be_continuityBased(beats,opBeats);
                    % fprintf('Exp-%d: Processing file... %s\n',ex, oFiles{k});
                    clear annout ann
                else
                    res(k).bpm = nan;
                    res(k).meter = nan;
                    res(k).rhythm = nan;
                    res(k).tempoCML = nan;
                    res(k).tempoAML = nan;
                    res(k).meterCorr = nan;
                    % Sama metrics
                    res(k).sama.pScore = nan;
                    res(k).sama.fMeas = nan; 
                    res(k).sama.precision = nan; 
                    res(k).sama.recall = nan; 
                    res(k).sama.Ameas = nan; 
                    res(k).sama.infoGain = nan;
                    res(k).sama.cmlC = nan;
                    res(k).sama.cmlT = nan;
                    res(k).sama.amlC = nan;
                    res(k).sama.amlT = nan;
                    % Beat metrics
                    res(k).beat.pScore = nan;
                    res(k).beat.fMeas = nan;
                    res(k).beat.precision = nan;
                    res(k).beat.recall = nan;
                    res(k).beat.Ameas = nan;
                    res(k).beat.infoGain = nan;
                    res(k).beat.cmlC = nan; 
                    res(k).beat.cmlT = nan; 
                    res(k).beat.amlC = nan; 
                    res(k).beat.amlT = nan;
                    fprintf('Exp-%d: Did not find file while processing file... %s\n',ex, filelist{k});
                end                
            end
            % First sort the results on fileID
            [fileID sortInd] = sort(fileID,'ascend');
            res = res(sortInd);
            talaIndicator = talaIndicator(sortInd);
            fileTalaName = fileTalaName(sortInd);
            fileNum = fileNum(sortInd);
            fname = filelist(sortInd);
            % Overall results for each tala
            sm = [res.sama]; s.pScore = [sm.pScore]; s.fMeas = [sm.fMeas];
            s.precision = [sm.precision]; s.recall = [sm.recall]; s.infoGain = [sm.infoGain];
            s.cmlC = [sm.cmlC]; s.cmlT = [sm.cmlT]; s.amlC = [sm.amlC]; s.amlT = [sm.amlT];
            bt = [res.beat]; b.pScore = [bt.pScore]; b.fMeas = [bt.fMeas];    
            b.precision = [bt.precision]; b.recall = [bt.recall]; b.infoGain = [bt.infoGain];
            b.cmlC = [bt.cmlC]; b.cmlT = [bt.cmlT]; b.amlC = [bt.amlC]; b.amlT = [bt.amlT];
            % Write overall results to a file
            colHead1 = 'File,Tala,MedianTempo,EstMeter,EstRhythm,EstTempo,sfMeas,sPrec,sRecall,';
            colHead2 = 'sCMLt,sAMLt,sInfoGain,bfMeas,bPrec,bRecall,bCMLt,bAMLt,bInfoGain';
            fp = fopen([talabasepath exptName '_' talaName{t} '_nPatts_' num2str(numPatts(p)) '_' num2str(ex) '_allResults.txt'], 'wt');
            fprintf(fp, '%s\n', [colHead1 colHead2]);
            for k = 1:length(filelist)
                fprintf(fp, '%s,%s,%.2f,%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n',...
                    fname{k},fileTalaName{k},res(k).bpmGT,res(k).meter,res(k).rhythm,...
                    res(k).bpm,s.fMeas(k),s.precision(k),s.recall(k),s.cmlT(k),s.amlT(k),s.infoGain(k),...
                    b.fMeas(k),b.precision(k),b.recall(k),b.cmlT(k),b.amlT(k),b.infoGain(k));
            end
            fclose(fp);
            sMat = [];
            bMat = [];
            samaResults(t).patt(p).expt(ex).allres = res;
            samaResults(t).patt(p).expt(ex).op = s;
            samaResults(t).patt(p).expt(ex).tempoCML = [res.tempoCML];
            samaResults(t).patt(p).expt(ex).tempoAML = [res.tempoAML];
            samaResults(t).patt(p).expt(ex).meterCorr = [res.meterCorr];
            samaResults(t).patt(p).expt(ex).fname = fname;
            beatResults(t).patt(p).expt(ex).op = b;
            beatResults(t).patt(p).expt(ex).tempoCML = [res.tempoCML];
            beatResults(t).patt(p).expt(ex).tempoAML = [res.tempoAML];
            beatResults(t).patt(p).expt(ex).meterCorr = [res.meterCorr];
            beatResults(t).patt(p).expt(ex).fname = fname;
            clear res s b oFiles fname fileTalaName fileID talaIndicator fileNum
        end
    end
    t
end
%% Now to generate the mean column for the excel sheet
opFullCol = [];
opFullTabCol = [];
opFullTabColTempo = [];
for m = 1:2
    if m == 2
        backupSamaResults = samaResults;
        samaResults = beatResults;
    end
    for p = 1:length(numPatts)
        fMeasAccFull = [];
        precAccFull = [];
        recallAccFull = [];
        cmltAccFull = [];
        amltAccFull = [];
        infoGainAccFull = [];
        tempoCMLAccfull = [];
        tempoAMLAccfull = [];
        meterCorrAccfull = [];
        for t = 1:length(talaName)
            for ex = 1:nExp
                % f-Meas
                fMeasVal = samaResults(t).patt(p).expt(ex).op.fMeas';
                fMeas(t+1,ex) = nanmean(fMeasVal);
                fMeasAcc(:,ex) = fMeasVal;
                % 
                precVal = samaResults(t).patt(p).expt(ex).op.precision';
                prec(t+1,ex) = nanmean(precVal);
                precAcc(:,ex) = precVal;
                % 
                recallVal = samaResults(t).patt(p).expt(ex).op.recall';
                recall(t+1,ex) = nanmean(recallVal);
                recallAcc(:,ex) = recallVal;
                % 
                cmltVal = samaResults(t).patt(p).expt(ex).op.cmlT';
                cmlt(t+1,ex) = nanmean(cmltVal);
                cmltAcc(:,ex) = cmltVal;
                %
                amltVal = samaResults(t).patt(p).expt(ex).op.amlT';
                amlt(t+1,ex) = nanmean(amltVal);
                amltAcc(:,ex) = amltVal;
                %
                infoGainVal = samaResults(t).patt(p).expt(ex).op.infoGain';
                infoGain(t+1,ex) = nanmean(infoGainVal);
                infoGainAcc(:,ex) = infoGainVal;
                %
                tempoCMLVal = samaResults(t).patt(p).expt(ex).tempoCML(:);
                tempocml(t+1,ex) = nanmean(tempoCMLVal);
                tempocmlacc(:,ex) = tempoCMLVal;
                %
                tempoAMLVal = samaResults(t).patt(p).expt(ex).tempoAML(:);
                tempoaml(t+1,ex) = nanmean(tempoAMLVal);
                tempoamlacc(:,ex) = tempoAMLVal;
                %
                meterCorrVal = samaResults(t).patt(p).expt(ex).meterCorr(:);
                metercorr(t+1,ex) = nanmean(meterCorrVal);
                metercorracc(:,ex) = meterCorrVal;
            end
            fMeasAccFull = [fMeasAccFull; fMeasAcc];
            precAccFull = [precAccFull; precAcc];
            recallAccFull = [recallAccFull; recallAcc];
            cmltAccFull = [cmltAccFull; cmltAcc];
            amltAccFull = [amltAccFull; amltAcc];
            infoGainAccFull = [infoGainAccFull; infoGainAcc];
            tempoCMLAccfull = [tempoCMLAccfull; tempocmlacc];
            tempoAMLAccfull = [tempoAMLAccfull; tempoamlacc];
            meterCorrAccfull = [meterCorrAccfull; metercorracc];
            clear fMeasAcc precAcc recallAcc cmltAcc amltAcc infoGainAcc tempocmlacc tempoamlacc metercorracc
        end
        fMeas(1,:) = nanmean(fMeasAccFull);
        prec(1,:) = nanmean(precAccFull);
        recall(1,:) = nanmean(recallAccFull);
        cmlt(1,:) = nanmean(cmltAccFull);
        amlt(1,:) = nanmean(amltAccFull);
        infoGain(1,:) = nanmean(infoGainAccFull);
        tempocml(1,:) = nanmean(tempoCMLAccfull);
        tempoaml(1,:) = nanmean(tempoAMLAccfull);
        metercorr(1,:) = nanmean(meterCorrAccfull);
        zz = zeros(1,nExp);
        opCol(:,p) = nanmean([fMeas; zz; prec; zz; recall; zz; cmlt; zz; amlt; zz; infoGain; zz; zz],2);       
        tabFmeas = [fMeas(2:end,:); fMeas(1,:)]./100;
        tabamlt = [amlt(2:end,:); amlt(1,:)]./100;
        tabinfoGain = [infoGain(2:end,:); infoGain(1,:)];
        tabtempoCML = [tempocml(2:end,:); tempocml(1,:)];
        tabtempoAML = [tempoaml(2:end,:); tempoaml(1,:)];
        tabmeterCorr = [metercorr(2:end,:); metercorr(1,:)];
        opTabCol(:,p) = nanmean([tabFmeas; tabamlt; tabinfoGain],2);
        opTabColTempo(:,p) = nanmean([tabtempoCML; tabtempoAML; tabmeterCorr],2);  % Overwritten, but does not matter
    end
    opFullCol = [opFullCol; opCol];
    tempMat1 = [];
    tempMat2 = [];
    nRows = length(talaName)+1;
    for kk = 1:3
        tempMat1 = [tempMat1 opTabCol((kk-1)*nRows+1:kk*nRows,:)];
        tempMat2 = [tempMat2 opTabColTempo((kk-1)*nRows+1:kk*nRows,:)];
    end
    opFullTabCol = [opFullTabCol; tempMat1];
    opFullTabColTempo = [opFullTabColTempo; tempMat2];
    if m == 2
        samaResults = backupSamaResults;
    end
    clear opCol opTabCol opTabColTempo tempMat1 tempMat2 backupSamaResults 
end
exptStub = [dName '_' task '_' exptName];
for p = 1:length(numPatts)
    dlmwrite([bpath exptName filesep 'opCol_' exptStub '_nPatt_' num2str(numPatts(p)) '.txt'], opFullCol(:,p), 'precision', '%.2f');
end
save([bpath exptName filesep exptStub '.mat'], 'samaResults', 'beatResults', 'exptStub');
% Now to create latex tables: Two tables are done: One with all talas, one
% showing just the average over all talas.
nMeas = 3;
nTalas = length(talaName);
nPatts = length(numPatts);
exptNameL = exptName;
exptNameL(exptNameL == '_') = '-';
% Mean table first
fp = fopen([bpath exptName filesep exptStub '_MEAN.tex'],'wt');
fprintf(fp,'\\begin{table}\n\\centering\n\\begin{tabular}{@{}lcccccccc@{}} \\toprule\n');
fprintf(fp,'Measure & \\multicolumn{2}{c}{\\fmeas} && \\multicolumn{2}{c}{\\amlt} && ');
fprintf(fp,'\\multicolumn{2}{c}{\\infoGain}\\tabularnewline \\addlinespace[2pt]\n');
fprintf(fp,'$\\nrhythmPatts$ ');
ss = repmat([sprintf('& %d ',numPatts) '&'],1,nMeas);
ss = ss(1:end-2);
fprintf(fp,'%s\\tabularnewline \\midrule\n\\Gls{sama} Tracking ',ss);
ss = [sprintf('& %.3f ',opFullTabCol(nTalas+1,  1:nPatts)), '&' , ...
    sprintf('& %.3f ',opFullTabCol(nTalas+1,  nPatts+1:nPatts*2)), '&', ...
    sprintf('& %.2f ',opFullTabCol(nTalas+1,  nPatts*2+1:3*nPatts))];
ss = ss(1:end-1);
fprintf(fp, '%s\\tabularnewline \\addlinespace[3pt] \n Beat Tracking ', ss);
ss = [sprintf('& %.3f ',opFullTabCol(2*(nTalas+1),  1:nPatts)), '&' , ...
    sprintf('& %.3f ',opFullTabCol(2*(nTalas+1),  nPatts+1:nPatts*2)), '&', ...
    sprintf('& %.2f ',opFullTabCol(2*(nTalas+1),  nPatts*2+1:3*nPatts))];
ss = ss(1:end-1);
fprintf(fp, '%s\\tabularnewline \\bottomrule\n\\end{tabular}\n\\caption[', ss);
fprintf(fp, '%s]{%s}\\label{tab:mtmeanres:%s}\n\\end{table}\n', ...
    [taskStub ' with ' exptNameL ' on ' dStub ' dataset'], ...
    ['Results of ' taskStub ' with ' exptNameL ' on ' dStub ' dataset'], ...
    [taskID dName exptNameL(1:3)]);
fclose(fp);
% All values table next
fp = fopen([bpath exptName filesep exptStub '_ALL.tex'],'wt');
fprintf(fp,'\\begin{table}\n\\centering\n\\begin{tabular}{@{}llcccccccc@{}} \\toprule\n');
fprintf(fp,' & Measure & \\multicolumn{2}{c}{\\fmeas} && \\multicolumn{2}{c}{\\amlt} && ');
fprintf(fp,'\\multicolumn{2}{c}{\\infoGain}\\tabularnewline \\addlinespace[2pt]\n');
fprintf(fp,' & $\\nrhythmPatts$ ');
ss = repmat([sprintf('& %d ',numPatts) '&'],1,nMeas);
ss = ss(1:end-2);
fprintf(fp,'%s\\tabularnewline \\midrule',ss);
tableNames = [tableOpName {'Mean'}];
strLabel = {['\multirow{' num2str(nTalas+1) '}{*}{\rotatebox[origin=c]{90}{\Gls{sama}}}'], ...
    ['\multirow{' num2str(nTalas+1) '}{*}{\rotatebox[origin=c]{90}{Beat}}']};
strEndLabel = {'\midrule', '\bottomrule'};
for k = 1:2
    for t = 1:nTalas
        if t == 1
            str1 = strLabel{k};
        else
            str1 = [];
        end
        fprintf(fp, '\n%s & %s ', str1, tableNames{t});
        ss = [sprintf('& %.3f ',opFullTabCol((k-1)*(nTalas+1)+t,  1:nPatts)), '&' , ...
            sprintf('& %.3f ',opFullTabCol((k-1)*(nTalas+1)+t,  nPatts+1:nPatts*2)), '&', ...
            sprintf('& %.2f ',opFullTabCol((k-1)*(nTalas+1)+t,  nPatts*2+1:3*nPatts))];
        ss = ss(1:end-1);
        fprintf(fp, '%s\\tabularnewline', ss);
    end
    fprintf(fp, ' \\addlinespace[2pt]\n & \\textbf{%s} ', tableNames{end});
    ss = [sprintf('& \\textbf{%.3f} ',opFullTabCol(k*(nTalas+1),  1:nPatts)), '&' , ...
        sprintf('& \\textbf{%.3f} ',opFullTabCol(k*(nTalas+1),  nPatts+1:nPatts*2)), '&', ...
        sprintf('& \\textbf{%.2f} ',opFullTabCol(k*(nTalas+1),  nPatts*2+1:3*nPatts))];
    ss = ss(1:end-1);
    fprintf(fp, '%s\\tabularnewline %s', ss, strEndLabel{k});
end
fprintf(fp,'\n\\end{tabular}\n\\caption[');
fprintf(fp, '%s]{%s}\\label{tab:mtres:%s}\n\\end{table}\n', ...
    [taskStub ' with ' exptNameL ' on ' dStub ' dataset'], ...
    ['Results of ' taskStub ' with ' exptNameL ' on ' dStub ' dataset'], ...
    [taskID dName exptNameL(1:3)]);
fclose(fp);
% A mean highlighted table next
fp = fopen([bpath exptName filesep exptStub '_ALL_highlight.tex'],'wt');
fprintf(fp,'\\begin{table}\n\\setlength{\\tabcolsep}{1.5\\tabcolsep}\n\\centering\n\\begin{tabular}{@{}LLCCCCCCCC@{}} \\toprule\n');   
fprintf(fp,' & Measure & \\multicolumn{2}{c}{\\fmeas} && \\multicolumn{2}{c}{\\amlt} && ');
fprintf(fp,'\\multicolumn{2}{c}{\\infoGain}\\tabularnewline \\addlinespace[2pt]\n');
fprintf(fp,' & $\\nrhythmPatts$ ');
ss = repmat([sprintf('& %d ',numPatts) '&'],1,nMeas);
ss = ss(1:end-2);
fprintf(fp,'%s\\tabularnewline \\midrule',ss);
tableNames = [tableOpName {'Mean'}];
strLabel = {['\multirow{-' num2str(nTalas+1) '}{*}{\rotatebox[origin=c]{90}{\Gls{sama}}}'], ...
    ['\multirow{-' num2str(nTalas+1) '}{*}{\rotatebox[origin=c]{90}{Beat}}']};
strEndLabel = {'\midrule', '\bottomrule'};
for k = 1:2
    for t = 1:nTalas
        fprintf(fp, '\n & %s ', tableNames{t});
        ss = [sprintf('& %.3f ',opFullTabCol((k-1)*(nTalas+1)+t,  1:nPatts)), '&' , ...
            sprintf('& %.3f ',opFullTabCol((k-1)*(nTalas+1)+t,  nPatts+1:nPatts*2)), '&', ...
            sprintf('& %.2f ',opFullTabCol((k-1)*(nTalas+1)+t,  nPatts*2+1:3*nPatts))];
        ss = ss(1:end-1);
        fprintf(fp, '%s\\tabularnewline', ss);
    end
    fprintf(fp, ' \\addlinespace[2pt]\n \\rowcolor{tabgray}[0.1pt][0.1pt] %s & \\textbf{%s} ', strLabel{k}, tableNames{end});
    ss = [sprintf('& %.3f ',opFullTabCol(k*(nTalas+1),  1:nPatts)), '&' , ...
        sprintf('& %.3f ',opFullTabCol(k*(nTalas+1),  nPatts+1:nPatts*2)), '&', ...
        sprintf('& %.2f ',opFullTabCol(k*(nTalas+1),  nPatts*2+1:3*nPatts))];
    ss = ss(1:end-1);
    fprintf(fp, '%s\\tabularnewline %s', ss, strEndLabel{k});
end
fprintf(fp,'\n\\end{tabular}\n\\caption[');
fprintf(fp, '%s]{%s}\\label{tab:mtres:%s}\n\\end{table}\n', ...
    [taskStub ' with ' exptNameL ' on ' dStub ' dataset'], ...
    ['Results of ' taskStub ' with ' exptNameL ' on ' dStub ' dataset'], ...
    [taskID dName exptNameL(1:3)]);
fclose(fp);
%% All values table of tempo and tala recognition
fp = fopen([bpath exptName filesep exptStub '_ALL_tempometer.tex'],'wt');
fprintf(fp,'\\begin{table}\n\\centering\n\\begin{tabular}{@{}lcccccccc@{}} \\toprule\n');
fprintf(fp,'Measure & \\multicolumn{2}{c}{tempoCML} && \\multicolumn{2}{c}{tempoAML} && ');
fprintf(fp,'\\multicolumn{2}{c}{talaID}\\tabularnewline \\addlinespace[2pt]\n');
fprintf(fp,'$\\nrhythmPatts$ ');
ss = repmat([sprintf('& %d ',numPatts) '&'],1,nMeas);
ss = ss(1:end-2);
fprintf(fp,'%s\\tabularnewline \\midrule',ss);
tableNames = [tableOpName {'Mean'}];
strLabel = {['\multirow{' num2str(nTalas+1) '}{*}{\rotatebox[origin=c]{90}{\Gls{sama}}}'], ...
    ['\multirow{' num2str(nTalas+1) '}{*}{\rotatebox[origin=c]{90}{Beat}}']};
strEndLabel = {'\midrule', '\bottomrule'};
for t = 1:nTalas
    fprintf(fp, '\n%s ', tableNames{t});
    ss = [sprintf('& %.3f ',opFullTabColTempo(t, 1:nPatts)), '&' , ...
        sprintf('& %.3f ',opFullTabColTempo(t, nPatts+1:nPatts*2)), '&', ...
        sprintf('& %.2f ',opFullTabColTempo(t, nPatts*2+1:3*nPatts))];
    ss = ss(1:end-1);
    fprintf(fp, '%s\\tabularnewline', ss);
end
fprintf(fp, ' \\addlinespace[2pt]\n\\textbf{%s} ', tableNames{end});
ss = [sprintf('& \\textbf{%.3f} ',opFullTabColTempo(nTalas+1, 1:nPatts)), '&' , ...
    sprintf('& \\textbf{%.3f} ',opFullTabColTempo(nTalas+1, nPatts+1:nPatts*2)), '&', ...
    sprintf('& \\textbf{%.2f} ',opFullTabColTempo(nTalas+1, nPatts*2+1:3*nPatts))];
ss = ss(1:end-1);
fprintf(fp, '%s\\tabularnewline \\bottomrule', ss);
fprintf(fp,'\n\\end{tabular}\n\\caption[');
fprintf(fp, '%s]{%s}\\label{tab:mttempores:%s}\n\\end{table}\n', ...
    [taskStub ' tempo and tala with ' exptNameL ' on ' dStub ' dataset'], ...
    ['Results of ' taskStub ' tempo and tala with ' exptNameL ' on ' dStub ' dataset'], ...
    [taskID dName exptNameL(1:3)]);
fclose(fp);
