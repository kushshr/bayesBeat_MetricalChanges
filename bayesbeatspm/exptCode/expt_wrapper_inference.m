% clear
close all
clc
if size(whos('serverFlag'),1) && serverFlag
    % Running on server
    disp('Running on server...')
    addpath('../../bayesbeatSPM');
    base _path = '/homedtic/amurthy/UPFWork_Server/PhD';
else
    clear; 
    serverFlag = 0;        % Default is that we are running on local machine
    addpath('../../bayesbeatSPM');
    base_path = 'D:\trainSet';
end
numExp = 1;
folds = 1;
procLongFiles = 1;      % Set it to 1 if you wish to process longer files
timeFormat = 'HH:MM:SS.FFF dd-mmm-yyyy';
frame_length = 0.02;
numPatts = [1];   % Number of patterns per meter
% The order of variables is this: 
% Variable order in configList: dataset, system, patt_trans_opt, peakInfMode, pattern_size
% dataset: 'CMCMDa_small', 'CMCMDa_v2', 'HMDf', 'HMDl', 'HMDs', 'BallroomDataset', 'Cretan'
% system: Tracking (0), tempoInfTracking (1), tempoSamaInfTracking (2), Inference (3)
% patt_trans_opt: noTrans (0), prior (1), mix (2), acc (3)
% peakInfMode: allHop (0) peakHop (1) fixHop (2) obsHop (3)
% pattern_size: 'bar', 'section', 'beat'
% configList = {{'BallroomDataset', 3, 1, 0, 'bar'}, ...
%    {'HMDs', 3, 1, 0, 'bar'}, ...
%    {'HMDl', 3, 1, 0, 'bar'}, ...
%    {'HMDf', 3, 1, 0, 'bar'}};

% configList = {{'CMCMDa_small', 3, 1, 0, 'bar', 'HMM'},...
%     {'BallroomDataset', 3, 1, 0, 'bar', 'HMM'}, ...
%     {'HMDs', 3, 1, 0, 'bar', 'HMM'}, ...
%     {'HMDf', 3, 1, 0, 'bar', 'HMM'}, ...
%     {'HMDl', 3, 1, 0, 'bar', 'HMM'}};

configList = {{'Metrical', 3, 1, 0, 'bar'}, ...
    };
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
        % Copy params 
        Params = allParams;  %  Copy all params and then modify whatever is needed
        Params.simStartTime = datestr(clock, timeFormat);
        Params.meters = allParams.meters;
        Params.meter_names = allParams.meter_names;
        Params.R = numPatts(r);
        Params.sections = allParams.sections;
        Params.sectionLens = allParams.sectionLens;
        Params.section_names = allParams.section_names;
        Params.min_tempo = allParams.min_tempo;
        Params.max_tempo = allParams.max_tempo;
        % Some param values need to change for this iteration
        if strcmp(Params.pattern_size,'section')
            secLens = Params.sectionLens{:} ./ Params.meters(2);
            Params.M = Params.Minit * 16 * max([secLens{:}]);
            clear secLens
        else
            % Get the longest cycle and scale it to that
            Params.M = Params.Minit * 4 * max(Params.meters(:,1)./Params.meters(:,2));
        end
        % Set a name to store the results
        if Params.inferenceMethod(1:2) == 'HM'
            disp('An exact inference using HMM chosen');
        elseif Params.inferenceMethod(1:2) == 'PF'
            fprintf('Approximate inference using a Particle Filter: %s\n', Params.store_name);
            if strcmp(Params.pattern_size,'section')
                Params.nParticles = 1500*Params.R*length(Params.sections{:}); % 1500 particles/pattern/section
            else
                Params.nParticles = 1500*Params.R*length(Params.meter_names); % 1500 particles/pattern/section
            end
        end
        % Now for many trials with diff folds
        for ex = 1:numExp
            for fld = 1:folds
                % specify a simulation id
                sim_id = 1000*ex+fld;
                Params.results_path = fullfile(allParams.results_path, Params.dataset, ...
                    Params.system, Params.store_name, ['nPatts_' num2str(Params.R)],...
                    num2str(sim_id));
                if ~isdir(Params.results_path)
                    mkdir(Params.results_path);
                end
                Params.train_set = ['train_' num2str(fld)];
                % Path to lab files
                Params.trainLab = fullfile(Params.base_path, 'Data', Params.dataset, ...
                    ['train_' num2str(fld) '.lab']);
                Params.testLab = fullfile(Params.base_path, 'Data', Params.dataset, ...
                    ['test_' num2str(fld) '.lab']);
                Params.testLab = 'D:\testSet\test_1.lab';
                Params.trainLab = 'D:\trainSet\train_1.lab';
                % CLUSTERING THE DATASET
                data_save_path = Params.results_path;
                Clustering = RhythmCluster(Params.trainLab, Params.feat_type, frame_length,...
                    data_save_path, Params.pattern_size);
                % cluster the dataset according to the meter of each file
                % Params.clusterIdFln = Clustering.make_cluster_assignment_file('meter');
                Clustering.make_feats_per_patt(Params.whole_note_div);
                %Clustering.make_feats_per_patt(Params.whole_note_div);
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
                            % dry run
                            results = BT.do_inference(k);
                            timePerFile(k) = results{end};
                            BT.save_results(results, Params.results_path, fname);
                        catch
                            fprintf('Did not process file, some error: %s\n', BT.test_data.file_list{k})
                        end
                    end
                    % close all
                    fprintf('%s: Inference with %d patterns: %s: %d/%d expt, %d/%d fold, %d/%d done...\n',...
                        [Params.dataset '_' Params.system '_' Params.store_name], ...
                        Params.R, fname, ex, numExp, fld,...
                        folds, k, length(BT.test_data.file_list));
                end
                expVals.results_path = Params.results_path;
                expVals.trainLab = Params.trainLab;
                expVals.testLab = Params.testLab;
                expVals.pieceDur = dur;
                expVals.runTime = timePerFile;
                expVals.minN = BT.model.minN;
                expVals.maxN = BT.model.maxN;
                expVals.Meff = BT.model.Meff;
                expVals.N = BT.model.N;
                expVals.M = BT.model.M;
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
    save(fullfile(allParams.results_path, allParams.dataset, allParams.system, allParams.store_name,'allParameters.mat'),'allParams');
    clear allParams
end
if serverFlag
    exit;
end
% 
% Done code
% configList = {{'BallroomDataset', 3, 1, 0, 'bar'}, ...
%     {'HMDs', 3, 1, 0, 'bar'}, ...
%     {'HMDl', 3, 1, 0, 'bar'}, ...
%     {'HMDf', 3, 1, 0, 'bar'}};

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD CODE %%%%%%%%%%%%%%%%%%%%%%%%%%
% % clear
% close all
% clc
% if size(whos('serverFlag'),1) && serverFlag
%     % Running on server
%     disp('Running on server...')
%     addpath('../../bayesbeatSPM');
%     basepath = '/homedtic/amurthy/UPFWork_Server/PhD';
% else
%     clear; 
%     serverFlag = 0;        % Default is that we are running on local machine
%     addpath('../../bayesbeatSPM');
%     basepath = '/media/Code/UPFWork/PhD';
% end
% numExp = 3;
% folds = 2;
% procLongFiles = 1;      % Set it to 1 if you wish to process longer files
% timeFormat = 'HH:MM:SS.FFF dd-mmm-yyyy';
% frame_length = 0.02; 
% for ex = 1:numExp
%     for fld = 1:folds
%         Params.startTime = datestr(clock, timeFormat);
%         % specify a simulation id
%         sim_id = 1000*ex+fld;
%         % Run config file based on the inference scheme
%         % For HMM
%         % Params = HMM_config(basepath);
%         % For PF
%         Params = PF_config(basepath);
%         if strcmp(Params.pattern_size,'section')
%             secLens = cellfun(@(x,y) (x./y), Params.sectionLens, num2cell(Params.meters(:,2)'), 'uniformOutput', 0);
%             Params.M = Params.Minit * 16 * max([secLens{:}]);
%             clear secLens
%         else
%             % Get the longest cycle and scale it to that
%             Params.M = Params.Minit * 4 * max(Params.meters(:,1)./Params.meters(:,2));
%         end
%         % Set a name to store the results
%         Params.store_name = [Params.store_name '_1500pp'];
%         if Params.inferenceMethod(1:2) == 'HM'
%             disp('An exact inference using HMM chosen');
%         elseif Params.inferenceMethod(1:2) == 'PF'
%             fprintf('Approximate inference using a Particle Filter: %s\n', Params.store_name);
%         end
%         Params.results_path = fullfile(Params.results_path, Params.dataset,...
%               'Inference', Params.store_name, num2str(sim_id));
%         if ~isdir(Params.results_path)
%             mkdir(Params.results_path);
%         end
%         Params.train_set = ['train_' num2str(fld)];
%         % Path to lab files
%         Params.trainLab = fullfile(Params.base_path, 'Data', Params.dataset, ...
%                     ['train_' num2str(fld) '.lab']);
%         Params.testLab = fullfile(Params.base_path, 'Data', Params.dataset, ...
%                     ['test_' num2str(fld) '.lab']);
%         % CLUSTERING THE DATASET
%         data_save_path = Params.results_path;
%         Clustering = RhythmCluster(Params.trainLab, Params.feat_type, frame_length,...
%             data_save_path, Params.pattern_size);
%         % cluster the dataset according to the meter of each file
%         % Params.clusterIdFln = Clustering.make_cluster_assignment_file('meter');
%         Clustering.make_feats_per_patt(Params.whole_note_div);
%         Params.clusterIdFln = Clustering.do_clustering(Params.pattern_size, Params.R, ...
%             'meters', Params.meters, 'meter_names', Params.meter_names,...
%             'sections', Params.sections, 'section_names', Params.section_names,...
%             'save_pattern_fig', Params.fig_store_flag, ...
%             'plotting_path', Params.results_path);                
%         % TRAINING THE MODEL
%         % create beat tracker object
%         BT = BeatTracker(Params, sim_id);
%         % set up training data
%         BT.init_train_data();
%         % set up test_data
%         BT.init_test_data();
%         % initialize probabilistic model
%         BT.init_model();
%         % train model
%         BT.train_model();
%         % TEST THE MODEL
%         % do beat tracking
%         for k = 1:length(BT.test_data.file_list)
%             [samps Fs] = wavread(BT.test_data.file_list{k},'size');
%             dur(k) = samps(1)/Fs;
%             [~, fname, ~] = fileparts(BT.test_data.file_list{k});
%             if (dur(k) > 600) && ~procLongFiles    % Piece shorter than 10 min
%                 fprintf('File too large: %s\n', BT.test_data.file_list{k})
%                 timePerFile(k) = NaN;
%             else
%                 try
%                     results = BT.do_inference(k);
%                     timePerFile(k) = results{end};
%                     BT.save_results(results, Params.results_path, fname);
%                 catch
%                     fprintf('Did not process file, some error: %s\n', BT.test_data.file_list{k})
%                     timePerFile(k) = NaN;
%                 end
%             end
%             close all
%             fprintf('%s: %d/%d expt, %d/%d fold, %d/%d done...\n',...
%                 fname, ex, numExp, fld,...
%                 folds, k, length(BT.test_data.file_list));
%         end
%         Params.pieceDur = dur;
%         Params.runTime = timePerFile;
%         Params.endTime = datestr(clock);
%         Params.minN = BT.model.minN;
%         Params.maxN = BT.model.maxN;
%         Params.Meff = BT.model.Meff;
%         clear dur timePerFile;
%         storeParams{ex,fld} = Params;
%     end
% end
% [path1, fname, ~] = fileparts(Params.results_path);
% save(fullfile(path1,'Parameters.mat'),'storeParams');
% if serverFlag
%     exit;
% end
