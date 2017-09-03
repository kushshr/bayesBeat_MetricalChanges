% clear
close all
clc
if size(whos('serverFlag'),1) && serverFlag
    % Running on server
    disp('Running on server...')
    addpath('../../bayesbeatSPM');
    base_path = '/homedtic/amurthy/UPFWork_Server/PhD';
else
    clear; 
    serverFlag = 0;        % Default is that we are running on local machine
    addpath('../../bayesbeatspm');
    %base_path = '/media/Code/UPFWork/PhD';
    base_path = 'D:\trainSet';
end
numExp = 1;
folds = 1;
procLongFiles = 1;      % Set it to 1 if you wish to process longer files
timeFormat = 'HH:MM:SS.FFF dd-mmm-yyyy';
frame_length = 0.02;
numPatts = [1 2];
% The order of variables is this: 
% Variable order in configList: dataset, system, patt_trans_opt, peakInfMode, pattern_size
% dataset: 'CMCMDa_small', 'CMCMDa_v2', 'HMDf', 'HMDl', 'HMDs', 'BallroomDataset', 'Cretan'
% system: Tracking (0), tempoInfTracking (1), tempoSamaInfTracking (2), Inference (3), samaInfTracking (4)
% patt_trans_opt: noTrans (0), prior (1), mix (2), acc (3)
% peakInfMode: allHop (0) peakHop (1) fixHop (2) obsHop (3)
% pattern_size: 'bar' , 'section', 'beat'
% infAlgo: 'HMM' or 'PF' (default is PF)
% Testing number of Particles with SPM
configList = {{'Metrical', 0, 1, 0, 'bar'}, ...
    };


% 
% configList = {{'CMCMDa_small', 0, 1, 0, 'section'}, ...
%    {'CMCMDa_small', 0, 1, 0, 'bar'}};
% SPM informed tracking expts
% configList = {{'CMCMDa_small', 4, 1, 0, 'section'}, ...
%      {'HMDs', 4, 1, 0, 'section'}, ...
%      {'HMDl', 4, 1, 0, 'section'}, ...
%      {'BallroomDataset', 1, 1, 0, 'section'}, ...
%      {'BallroomDataset', 2, 1, 0, 'section'}, ...
%      {'BallroomDataset', 4, 1, 0, 'section'}, ...
%      {'TurkishMakam', 1, 1, 0, 'section'}, ...
%      {'TurkishMakam', 2, 1, 0, 'section'}, ...
%      {'TurkishMakam', 4, 1, 0, 'section'}};
% 15/07/2016: Turkish music expts
% configList = {{'TurkishMakam', 0, 1, 0, 'bar'}, ...
%     {'TurkishMakam', 0, 1, 0, 'section'}, ...
%     {'TurkishMakam', 1, 1, 0, 'bar'}, ...
%     {'TurkishMakam', 2, 1, 0, 'bar'}, ...
%     {'TurkishMakam', 4, 1, 0, 'bar'}};
% 15/07/2016: Journal paper expts
% configList = {{'CMCMDa_small', 4, 1, 0, 'bar'}, ...
%     {'HMDs', 4, 1, 0, 'bar'}, ...
%     {'HMDl', 4, 1, 0, 'bar'}, ...
%     {'BallroomDataset', 4, 1, 0, 'bar'}, ...
%     {'BallroomDataset', 1, 1, 0, 'bar'}, ...
%     {'BallroomDataset', 2, 1, 0, 'bar'}};
% configList = {{'BallroomDataset', 0, 1, 0, 'section'}, ...
%     };

% 28/02/2016: local rerun HMDs
% configList = {{'HMDs', 1, 1, 0, 'section'}, ...
%     {'HMDs', 2, 1, 0, 'section'}, ...
%     {'HMDs', 0, 1, 0, 'section'}, ...
%     {'HMDf', 1, 1, 0, 'section'}, ...
%     {'HMDf', 2, 1, 0, 'section'}};
% 
% % 27/02/2016: local rerun
% configList = {{'BallroomDataset', 0, 1, 0, 'bar'}, ...
%     {'BallroomDataset', 0, 1, 0, 'bar', 'HMM'}, ...
%     {'BallroomDataset', 0, 2, 0, 'bar'}, ...
%     {'BallroomDataset', 0, 1, 3, 'bar'}, ...
%     {'BallroomDataset', 0, 1, 1, 'bar'}, ...
%     {'BallroomDataset', 0, 3, 0, 'bar'}};

%{'CMCMDa_small', 0, 1, 0, 'bar', 'HMM'}, ...
%{'BallroomDataset', 0, 1, 0, 'bar', 'HMM'}, ...
%{'HMDs', 0, 1, 0, 'bar', 'HMM'}, ...
%{'HMDf', 0, 1, 0, 'bar', 'HMM'}, ...
              
for cf = 1:length(configList)
    dataset = configList{cf}{1};
    def_tempoInfMode = configList{cf}{2};
    def_patt_trans_opt = configList{cf}{3};
    def_peakInfMode = configList{cf}{4};
    def_pattern_size = configList{cf}{5};
    if length(configList{cf}) > 5 && strcmp(configList{cf}{6},'HMM')
        allParams = HMM_config(base_path, dataset,def_tempoInfMode, def_patt_trans_opt, ...
            def_peakInfMode, def_pattern_size);
    else 
        allParams = PF_config(base_path, dataset, def_tempoInfMode, def_patt_trans_opt, ...
            def_peakInfMode, def_pattern_size);
    end
    % Now start with the simulation
    nTalas = length(allParams.meter_names);
    talaIDs = allParams.meter_names;
    for r = 1:length(numPatts)
        for t = 1:nTalas
            % Copy params
            Params = allParams;  %  Copy all params and then modify whatever is needed
            Params.R = numPatts(r);
            Params.simStartTime = datestr(clock, timeFormat);
            Params.meters = allParams.meters(t,:);
            Params.meter_names = allParams.meter_names(t);
            Params.sections = allParams.sections(t);
            Params.sectionLens = allParams.sectionLens(t);
            Params.section_names = allParams.section_names(t);
            Params.min_tempo = allParams.min_tempo(t);
            Params.max_tempo = allParams.max_tempo(t);
            % Some param values need to change for this iteration
            if strcmp(Params.pattern_size,'section')
                secLens = Params.sectionLens{:} ./ Params.meters(2);
                Params.M = Params.Minit * 4 * max(secLens);
                clear secLens
            else
                % Get the longest cycle and scale it to that
                Params.M = Params.Minit * 4 * Params.meters(1) / Params.meters(2);
            end
            % Set a name to store the results
            if Params.inferenceMethod(1:2) == 'HM'
                disp('An exact inference using HMM chosen');
            elseif Params.inferenceMethod(1:2) == 'PF'
                fprintf('Approximate inference using a Particle Filter: %s\n', Params.store_name);
                if strcmp(Params.pattern_size,'section')
                    % Params.nParticles = 1500*Params.R*length(Params.sections{:}); % 1500 particles/pattern/section
                    Params.nParticles = 1500*Params.R; % 1500 particles in total
                else
                    Params.nParticles = 500; %*Params.R*length(Params.sections{:}); % Params.nParticles/nTalas*numPatts(r);  % 1500 particles/pattern/section
                end
            end
            % Now for many trials with diff folds
            for ex = 1:numExp
                for fld = 1:folds
                    % specify a simulation id
                    sim_id = 1000*ex+fld;
                    Params.results_path = fullfile(allParams.results_path, Params.dataset, ...
                        Params.system, Params.store_name, Params.meter_names{1}, ...
                        ['nPatts_' num2str(Params.R)], num2str(sim_id));
                    if ~isdir(Params.results_path)
                        mkdir(Params.results_path);
                    end
                    Params.train_set = ['train_' num2str(fld) '_' Params.meter_names{1}];
                    % Path to lab files
                    Params.trainLab = fullfile(Params.base_path, 'Data', Params.dataset, ...
                        ['train_' num2str(fld) '_' Params.meter_names{1} '.lab']);
                    Params.testLab = fullfile(Params.base_path, 'Data', Params.dataset, ...
                        ['test_' num2str(fld) '_' Params.meter_names{1} '.lab']);
                    Params.testLab = 'D:\trainSet\train_1.lab';
                    Params.trainLab = 'D:\trainSet\train_1.lab';
                    % CLUSTERING THE DATASET
                    data_save_path = Params.results_path;
                    Clustering = RhythmCluster(Params.trainLab, Params.feat_type, frame_length,...
                        data_save_path, Params.pattern_size);
                    % cluster the dataset according to the meter of each file
                    % Params.clusterIdFln = Clustering.make_cluster_assignment_file('meter');
                    Clustering.make_feats_per_patt(Params.whole_note_div);
                    Params.clusterIdFln = Clustering.do_clustering(Params.pattern_size, Params.R, ...
                        'meters', Params.meters, 'meter_names', Params.meter_names,...
                        'sections', Params.sections, 'section_names', Params.section_names,...
                        'save_pattern_fig', Params.fig_store_flag, ...
                        'plotting_path', Params.results_path);
                    % TRAINING THE MODEL
                    % create beat tracker object
                    BT = BeatTracker(Params, sim_id);
                    % set up training data
                    BT.init_train_data();
                    % set up test_data
                    BT.init_test_data();
                    % initialize probabilistic model
                    BT.init_model();
                    % train model
                    BT.train_model();
                    % TEST THE MODEL
                    % do beat tracking
                    dur = zeros(1,length(BT.test_data.file_list));
                    timePerFile = NaN(1,length(BT.test_data.file_list));
                    for k = 1:length(BT.test_data.file_list)
                        % [samps Fs] = wavread(BT.test_data.file_list{k},'size');
                        aInfo = audioinfo(BT.test_data.file_list{k});
                        dur(k) = aInfo.Duration;
                        [~, fname, ~] = fileparts(BT.test_data.file_list{k});
                        if (dur(k) > 600) && ~procLongFiles    % Piece shorter than 10 min
                            fprintf('File too large: %s\n', BT.test_data.file_list{k})
                        else
                            try
                                results = BT.do_inference(k);
                                timePerFile(k) = results{end};
                                BT.save_results(results, Params.results_path, fname);
                            catch
                                fprintf('Did not process file, some error: %s\n', BT.test_data.file_list{k})
                            end
                        end
                        % close all
                        fprintf('%s: %s tala with %d patterns: %s: %d/%d expt, %d/%d fold, %d/%d done...\n',...
                            [Params.dataset '_' Params.system '_' Params.store_name], ...
                            Params.meter_names{1}, Params.R, fname, ex, numExp, fld,...
                            folds, k, length(BT.test_data.file_list));
                    end
                    expVals.results_path = Params.results_path;
                    expVals.trainLab = Params.trainLab;
                    expVals.testLab = Params.testLab;
                    expVals.pieceDur = dur;
                    expVals.runTime = timePerFile;
                    expVals.Meff = BT.model.Meff;
                    expVals.N = BT.model.N;
                    expVals.M = BT.model.M;
                    expVals.minN = BT.model.minN;
                    expVals.maxN = BT.model.maxN;
                    if Params.inferenceMethod(1:2) == 'PF'
                        expVals.nParticles = Params.nParticles;
                    end
                    Params.storeExpVals{ex,fld} = expVals;
                    clear dur timePerFile;
                    close all
                    fclose('all');
                end
            end
            Params.simEndTime = datestr(clock, timeFormat);
            [path1, fname, ~] = fileparts(Params.results_path);
            save(fullfile(path1,'Parameters.mat'),'Params');
            clear Params expVals BT results
        end
    end
    save(fullfile(allParams.results_path, allParams.dataset, allParams.system, allParams.store_name,'allParameters.mat'),'allParams');
    clear allParams
end
if serverFlag
    exit;
end
% Server: 20/01/2016
% configList = {{'CMCMDa_small', 0, 1, 0, 'bar'}, ...
%     {'CMCMDa_small', 0, 2, 0, 'bar'}, ...
%     {'BallroomDataset', 0, 1, 0, 'bar'}, ... % Did not run
%     {'CMCMDa_small', 0, 1, 0, 'section'}};
% 28/01/2016
% configList = {{'HMDl', 0, 2, 0, 'bar'}, ...
%    {'HMDs', 0, 2, 0, 'bar'}, ...
%    {'HMDf', 0, 2, 0, 'bar'}};
% 29/01/2016
% configList = {{'CMCMDa_small', 0, 1, 3, 'bar'}, ...
%     {'HMDs', 0, 1, 3, 'bar'}, ...
%     {'BallroomDataset', 0, 1, 3, 'bar'}};

% Local: 20/01/2016
% configList = {{'CMCMDa_small', 1, 1, 0, 'bar'}, ...
%     {'HMDl', 1, 1, 0, 'bar'}, ...
%     {'HMDs', 1, 1, 0, 'bar'}, ...
%     {'HMDf', 1, 1, 0, 'bar'}};
% 22/01/2016
% configList = {{'CMCMDa_small', 2, 1, 0, 'bar'}, ...
%     {'HMDl', 2, 1, 0, 'bar'}, ...
%     {'HMDs', 2, 1, 0, 'bar'}, ...
%     {'HMDf', 2, 1, 0, 'bar'}};
% 23/01/2016
% configList = {{'CMCMDa_small', 2, 1, 0, 'section'}, ...
%     {'HMDl', 2, 1, 0, 'section'}, ...
%     {'HMDs', 2, 1, 0, 'section'}, ...
%     {'HMDf', 2, 1, 0, 'section'}};
% 28/01/2016
% configList = {{'CMCMDa_small', 1, 1, 0, 'section'}, ...
%    {'HMDl', 1, 1, 0, 'section'}, ...
%    {'HMDs', 1, 1, 0, 'section'}, ...
%    {'HMDf', 1, 1, 0, 'section'}};
% 28/01/2016
% configList = {{'CMCMDa_small', 0, 3, 0, 'bar'}, {'HMDs', 0, 3, 0, 'bar'}, ...
%   {'BallroomDataset', 0, 3, 0, 'bar'}};
