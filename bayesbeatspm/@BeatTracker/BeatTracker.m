classdef BeatTracker < handle
    % Beat tracker Class
    properties
        input_fln               % input filename (.wav or feature file)
        model                   % probabilistic model
        inferenceMethod         % forward, viterbi, ...
        feature                 
        train_data
        test_data
        sim_dir                 % directory where results are saved
        temp_path
        init_model_fln          % fln of initial model to start with
        Params                  % parameters from config file
    end
    
    methods
        function obj = BeatTracker(Params, sim_id)
            % parse parameters and set defaults
            if ~isfield(Params, 'inferenceMethod')
                obj.inferenceMethod = 'HMM_viterbi';
            else
                obj.inferenceMethod = Params.inferenceMethod;
            end
            if exist('sim_id', 'var')
                % obj.sim_dir = fullfile(Params.results_path, num2str(sim_id));
                obj.sim_dir = fullfile(Params.results_path);
                if ~exist(obj.sim_dir, 'dir')
                    system(['mkdir ', obj.sim_dir]);
                end
            end
            if isfield(Params, 'temp_path')
                obj.temp_path = Params.temp_path;
            end
            obj.Params = Params;
            % Set default values if not specified otherwise
            if ~isfield(obj.Params, 'learn_tempo_ranges')
                obj.Params.learn_tempo_ranges = 1;
            end
            if ~isfield(obj.Params, 'same_tempo_per_meter')
                obj.Params.same_tempo_per_meter = 0;
            end
            if ~isfield(obj.Params, 'min_tempo')
                obj.Params.min_tempo = 60;
            end
            if ~isfield(obj.Params, 'max_tempo')
                obj.Params.max_tempo = 220;
            end
            if ~isfield(obj.Params, 'save_beats')
                obj.Params.save_beats = 1;
            end
            if ~isfield(obj.Params, 'save_downbeats')
                obj.Params.save_downbeats = 1;
            end
            if ~isfield(obj.Params, 'save_section_times')
                obj.Params.save_section_times = 1;
            end
            if ~isfield(obj.Params, 'save_tempo')
                obj.Params.save_tempo = 0;
            end
            if ~isfield(obj.Params, 'save_tempo_seq')
                obj.Params.save_tempo_seq = 0;
            end
            if ~isfield(obj.Params, 'save_rhythm')
                obj.Params.save_rhythm = 0;
            end
            if ~isfield(obj.Params, 'save_rhythm_seq')
                obj.Params.save_rhythm_seq = 0;
            end            
            if ~isfield(obj.Params, 'save_section_names')
                obj.Params.save_section_names = 0;
            end
            if ~isfield(obj.Params, 'save_section_seq')
                obj.Params.save_section_seq = 0;
            end            
            if ~isfield(obj.Params, 'save_meter')
                obj.Params.save_meter = 0;
            end
            if ~isfield(obj.Params, 'save_features_to_file')
                obj.Params.save_features_to_file = 0;
            end
            if ~isfield(obj.Params, 'load_features_from_file')
                obj.Params.load_features_from_file = 1;
            end
            if ~isfield(obj.Params, 'transition_model_type')
                obj.Params.transition_model_type = '2015';
            end
            if ~isfield(obj.Params, 'alpha')
                obj.Params.alpha = 100;
            end
            if ~isfield(obj.Params, 'use_silence_state')
                obj.Params.use_silence_state = 0;
            end
            if ~isfield(obj.Params, 'outlier_percentile')
                obj.Params.outlier_percentile = 5;
            end
        end
        
        
        function init_model(obj)
            if isfield(obj.Params, 'model_fln') && ~isempty(obj.Params.model_fln)
                if exist(obj.Params.model_fln, 'file')
                    c = load(obj.Params.model_fln);
                    fields = fieldnames(c);
                    obj.model = c.(fields{1});
                    obj.model = obj.convert_to_new_model_format(obj.model);
                    obj.init_model_fln = obj.Params.model_fln;
                    obj.feature = Feature(obj.model.obs_model.feat_type, ...
                        obj.model.frame_length);
                    if isfield(obj.Params, 'use_mex_viterbi')
                        obj.model.use_mex_viterbi = obj.Params.use_mex_viterbi;
                    end
                    fprintf('* Loading model from %s\n', obj.Params.model_fln);
                else
                    error('Model file %s not found', obj.Params.model_fln);
                end
            else
                obj.feature = Feature(obj.Params.feat_type, ...
                    obj.Params.frame_length);
                if isempty(obj.train_data) % no training data given -> set defaults
                    obj.train_data.rhythm_names = cellfun(@(x) num2str(x), ...
                        num2cell(1:obj.Params.R), 'UniformOutput', false);
                end
                
                switch obj.Params.inferenceMethod(1:2)
                    case 'HM'
                        obj.model = HMM(obj.Params, obj.train_data);
                    case 'PF'
                        obj.model = PF(obj.Params, obj.train_data);
                    otherwise
                        error('BeatTracker.init_model: inference method %s not known', ...
                            obj.Params.inferenceMethod);
                end
            end
            
        end
        
        
        function init_train_data(obj)
            fprintf('* Set up training data ...\n');
            obj.train_data = Data(obj.Params.trainLab, 1);
            if ~isfield(obj.Params, 'clusterIdFln'), return;  end
            obj.train_data = obj.train_data.read_pattern_bars(...
                obj.Params.clusterIdFln, obj.Params.pattern_size);
            if ~isempty(obj.train_data.pr)
                % if pr is given in cluster file save it here
                obj.Params.pr = obj.train_data.pr;
            end
            if ~isempty(obj.train_data.prprior)
                % if prprior is given in cluster file save it here
                obj.Params.prprior = obj.train_data.prprior;
            end
            if ~isempty(obj.train_data.B_pr)
                % if B_pr is given in cluster file save it here
                obj.Params.B_pr = obj.train_data.B_pr;
            end
            if ~isempty(obj.train_data.B_prprior)
                % if B_prprior is given in cluster file save it here
                obj.Params.B_prprior = obj.train_data.B_prprior;
            end
            % make filename of features
            [~, clusterFName, ~] = fileparts(obj.Params.clusterIdFln);
            featStr = '';
            for iDim = 1:obj.Params.featureDim
                featType = strrep(obj.Params.feat_type{iDim}, '.', '-');
                featStr = [featStr, featType];
            end
            % featuresFln = fullfile(obj.Params.data_path, ['features_', ...
            %    clusterFName, '_', featStr, '.mat']);
            featuresFln = fullfile(obj.Params.results_path, ...
                [clusterFName, '_', featStr, '.mat']);
            barGrid_eff = obj.Params.whole_note_div * ...
                (obj.train_data.rhythm2meter(:, 1) ./ ...
                obj.train_data.rhythm2meter(:, 2));
            obj.train_data = ...
                obj.train_data.extract_feats_per_file_pattern_barPos_dim(...
                obj.Params.whole_note_div, barGrid_eff, ...
                obj.Params.featureDim, featuresFln, obj.Params.feat_type, ...
                obj.Params.frame_length, ...
                obj.Params.reorganize_bars_into_cluster);
            % process silence data
            if obj.Params.use_silence_state
                fid = fopen(obj.Params.silence_lab, 'r');
                silence_files = textscan(fid, '%s\n'); silence_files = ...
                    silence_files{1};
                fclose(fid);
                obj.train_data.feats_silence = [];
                for iFile=1:length(silence_files)
                    obj.train_data.feats_silence = ...
                        [obj.train_data.feats_silence;  ...
                        obj.feature.load_feature(silence_files{iFile})];
                end
            end
        end
        
        function init_test_data(obj)
            % create test_data object
            obj.test_data = Data(obj.Params.testLab, 0);
            if isfield(obj.Params, 'test_annots_folder')
                obj.test_data = obj.test_data.set_annots_path(obj.Params.test_annots_folder);
            end
            % Read ground truth ann for tempo informed tracking
            obj.test_data.read_gtTempo;
            obj.test_data.read_gtAnn;
        end
        
        function train_model(obj)
            if isempty(obj.init_model_fln)
                if length(obj.Params.min_tempo) == size(obj.Params.meters, 1);
                    min_tempo_param_per_rhythm = obj.Params.min_tempo(obj.train_data.rhythm2meterID);
                elseif length(obj.Params.min_tempo) == 1
                    min_tempo_param_per_rhythm = repmat(obj.Params.min_tempo, 1, obj.Params.R);
                else
                    disp(strcat(['BeatTracker/train_model: min_tempo incorrectly specified.'],...
                        ['Specify min_tempo with one value per meter, or just one value for all meters.']));
                    disp('Taking only the first value of min_tempo...');
                    min_tempo_param_per_rhythm = repmat(obj.Params.min_tempo(1),1,obj.Params.R);
                end                
                if length(obj.Params.max_tempo) == size(obj.Params.meters, 1);
                    max_tempo_param_per_rhythm = obj.Params.max_tempo(obj.train_data.rhythm2meterID);
                elseif length(obj.Params.max_tempo) == 1
                    max_tempo_param_per_rhythm = repmat(obj.Params.max_tempo,1,obj.Params.R);
                else
                    disp(strcat(['BeatTracker/train_model: max_tempo incorrectly specified.'],...
                        ['Specify max_tempo with one value per meter, or just one value for all meters.']));
                    disp('Taking only the first value of max_tempo...');
                    max_tempo_param_per_rhythm = repmat(obj.Params.max_tempo(1),1,obj.Params.R);
                end
                if obj.Params.learn_tempo_ranges
                    % get tempo ranges from data for each file
                    [tempo_min_per_cluster, tempo_max_per_cluster] = ...
                        obj.train_data.get_tempo_per_cluster(...
                        obj.Params.outlier_percentile);
                    % find min/max for each pattern
                    tempo_min_per_cluster = min(tempo_min_per_cluster)';
                    tempo_max_per_cluster = max(tempo_max_per_cluster)';
                    % restrict ranges
                    lowInd = tempo_min_per_cluster < min_tempo_param_per_rhythm(:);
                    highInd = tempo_max_per_cluster > max_tempo_param_per_rhythm(:);
                    tempo_min_per_cluster(lowInd) = min_tempo_param_per_rhythm(lowInd);
                    tempo_max_per_cluster(highInd) = max_tempo_param_per_rhythm(highInd);
                else
                    tempo_min_per_cluster = min_tempo_param_per_rhythm;
                    tempo_max_per_cluster = max_tempo_param_per_rhythm;
                end
                
                if obj.Params.same_tempo_per_meter
                    nMeters = unique(obj.train_data.rhythm2meterID);
                    for iR = 1:length(nMeters)
                        tempo_min_per_cluster(obj.train_data.rhythm2meterID == iR) = ...
                            min(tempo_min_per_cluster(obj.train_data.rhythm2meterID == iR));
                        tempo_max_per_cluster(obj.train_data.rhythm2meterID == iR) = ...
                            max(tempo_max_per_cluster(obj.train_data.rhythm2meterID == iR));
                    end
                end
                tempo_min_per_cluster = tempo_min_per_cluster(:);
                tempo_max_per_cluster = tempo_max_per_cluster(:);
                % fprintf('%s\n',sprintf('%.2f ',tempo_min_per_cluster));
                % fprintf('%s\n',sprintf('%.2f ',tempo_max_per_cluster));  
                fprintf('* Set up transition model\n');
                obj = obj.train_transition_model(tempo_min_per_cluster, ...
                    tempo_max_per_cluster);
                fprintf('* Set up observation model\n');
                if obj.Params.use_silence_state
                    obj.model = obj.model.make_observation_model(obj.train_data);
                else
                    obj.model = obj.model.make_observation_model(obj.train_data);
                end
                obj.model = obj.model.make_initial_distribution();
                
                fln = fullfile(obj.sim_dir, 'model.mat');
                switch obj.inferenceMethod(1:2)
                    case 'HM'
                        hmm = obj.model;
                        save(fln, 'hmm');
                    case 'PF'
                        pf = obj.model;
                        save(fln, 'pf');
                end
                fprintf('* Saved model (Matlab) to %s\n', fln);
            end            
            %  obj.model.save_hmm_data_to_hdf5('~/diss/src/matlab/beat_tracking/bayes_beat/data/filip/');
            
            if isfield(obj.Params, 'viterbi_learning_iterations') && ...
                    obj.Params.viterbi_learning_iterations > 0
                %    obj.model.trans_model = TransitionModel(obj.model.M, obj.model.Meff, obj.model.N, obj.model.R, obj.model.pn, obj.model.pr, ...
                %       obj.model.rhythm2meter_state, ones(1, obj.model.R), ones(1, obj.model.R)*obj.model.N);
                obj.refine_model(obj.Params.viterbi_learning_iterations);
            end
        end
        
        
        function obj = train_transition_model(obj, tempo_min_per_cluster, ...
                tempo_max_per_cluster)
            pr = obj.train_data.pr;
            prprior = obj.train_data.prprior;
            switch obj.inferenceMethod(1:2)
                case 'HM'
                    % To be corrected!!
                    obj.model = obj.model.make_transition_model(...
                        tempo_min_per_cluster, tempo_max_per_cluster, ...
                        obj.Params.alpha, obj.Params.pn, pr, prprior);
                case 'PF'
                    obj.model = obj.model.make_transition_model(...
                        tempo_min_per_cluster, tempo_max_per_cluster, ...
                        obj.Params, obj.train_data);
            end
        end
        
        function test_file_ids = retrain_model(obj, test_files_to_exclude)
            fprintf('    Retraining observation model ');
            if length(test_files_to_exclude) == 1
                % leave one out: get pattern idx of test file to only
                % retrain the remaining patterns
                r_i = unique(obj.train_data.bar2cluster(...
                    obj.train_data.bar2file == exclude_test_file_id));
            else
                % retrain all patterns
                r_i = 1:obj.model.R;
            end
            % Get file indices of test files within the training set
            file_idx = zeros(length(obj.train_data.file_list), 1);
            for test_id = 1:length(test_files_to_exclude)
                for train_id = 1:length(obj.train_data.file_list)
                    if strfind(obj.train_data.file_list{train_id}, ...
                            test_files_to_exclude{test_id});
                        file_idx(train_id) = 1;
                        break;
                    end
                end
            end
            test_file_ids = find(file_idx)';
            % exclude test files from training:
            obj.model = ...
                obj.model.retrain_observation_model(...
                obj.train_data.feats_file_pattern_barPos_dim(~file_idx, :, ...
                :, :), r_i);
            fprintf('done\n');
        end
        
        function refine_model(obj, iterations)
            if ~isempty(obj.init_model_fln) && ~isempty(strfind(obj.init_model_fln, '-'))
                dash = strfind(obj.init_model_fln, '-');
                iter_start = str2double(obj.init_model_fln(dash(end)+1:end-4)) + 1;
            else
                iter_start = 1;
            end
            % read annotations and add them to train_data object:
            obj.train_data.read_beats;
            obj.train_data.read_meter;
            for i = iter_start:iter_start+iterations-1
                fprintf('* Viterbi training: iteration %i\n', i);
                [obj.model, bar2cluster] = obj.model.viterbi_training(obj.feature, obj.train_data);
                hmm = obj.model;
                save(fullfile(obj.sim_dir, ['hmm-', obj.train_data.dataset, '-', num2str(i), '.mat']), 'hmm');
                save(fullfile(obj.sim_dir, ['bar2cluster-', obj.train_data.dataset, '-', num2str(i), '.mat']), 'bar2cluster');
            end
        end
        
        function results = do_inference(obj, test_file_id, do_output)
            if ~exist('do_output', 'var')
                do_output = 1;
            end
            [~, fname, ~] = fileparts(obj.test_data.file_list{test_file_id});
            % load feature
            observations = obj.feature.load_feature(...
                obj.test_data.file_list{test_file_id}, ...
                obj.Params.save_features_to_file, ...
                obj.Params.load_features_from_file);
            % Read the ground truth tempo value if given, for tempo
            % informed tracking
            gtbpm.val = obj.test_data.gtTempo{test_file_id};
            gtbpm.range = obj.Params.tempoInfMargin;
            % Also read the ground truth sama examples given if any
            gtAnn = obj.test_data.gtAnn{test_file_id};
            
            % compute observation likelihoods
            tic;
            results = obj.model.do_inference(observations, fname, ...
                obj.inferenceMethod, do_output, gtbpm, gtAnn);
            results{end+1} = toc;
            if obj.model.save_inference_data
                % save state sequence of annotations to file
                annot_fln = strrep(strrep(obj.test_data.file_list{test_file_id}, ...
                    'wav', 'beats'), 'audio', 'annotations/beats');
                if exist(annot_fln, 'file')
                    annots = load(annot_fln);
                    meter = max(annots(:, 2));
                    r = find(obj.model.rhythm2meter(:, 1) == meter);
                    if isempty(r)
                        fprintf(['    Cannot compute true path, file not', ...
                            'in test_data included ...\n']);
                    else
                        [m, n] = HMM.getpath(obj.model.Meff(r), annots, ...
                            obj.model.frame_length, size(observations, 1));
                        anns = [m, n, ones(length(m), 1) * r];
                        save(['/tmp/', fname, '_anns.mat'], 'anns');
                    end
                end
            end
            if do_output
                fprintf('    Real time factor: %.2f\n', toc / (size(observations, 1) * obj.feature.frame_length));
            end
        end
        
        
        function load_model(obj, fln)
            temp = load(fln);
            names = fieldnames(temp);
            obj.model = temp.(names{1});
        end
        
        
        function [] = save_results(obj, results, save_dir, fname)
            if ~exist(save_dir, 'dir')
                system(['mkdir ', save_dir]);
            end
            if obj.Params.save_beats
                BeatTracker.save_beats(results{1}, fullfile(save_dir, ...
                    [fname, '.beats.txt']));
            end
            if obj.Params.save_downbeats
                BeatTracker.save_downbeats(results{1}, fullfile(save_dir, ...
                    [fname, '.downbeats.txt']));
            end
            if obj.Params.save_section_times
                BeatTracker.save_section_times(results{1}, fullfile(save_dir, ...
                    [fname, '.section.txt']));
            end
            if obj.Params.save_tempo
                BeatTracker.save_tempo(results{2}, fullfile(save_dir, ...
                    [fname, '.bpm.txt']));
            end
            if obj.Params.save_tempo_seq
                BeatTracker.save_tempo_seq(results{6}, results{2}, fullfile(save_dir, ...
                    [fname, '.bpm.seq']));
            end
            if obj.Params.save_meter
                BeatTracker.save_meter(results{3}, fullfile(save_dir, ...
                    [fname, '.meter.txt']));
            end
            if obj.Params.save_rhythm
                BeatTracker.save_rhythm(results{4}, fullfile(save_dir, ...
                    [fname, '.rhythm.txt']), obj.model.rhythm_names);
            end
            if obj.Params.save_section_names
                BeatTracker.save_rhythm(results{5}, fullfile(save_dir, ...
                    [fname, '.secname.txt']), obj.model.pattInfo.name);
            end
            if (sum(results{1}(:,2) == 1) >= 2)         % Atleast two full cycles need to be present to get sequence of rhythm and sections 
                if obj.Params.save_rhythm_seq
                    BeatTracker.save_rhythm_seq(results{6}, results{4}, ...
                        results{1}(results{1}(:,2) == 1,1), ...
                        fullfile(save_dir, [fname, '.rhythm.seq']), obj.model.rhythm_names);
                end
                if obj.Params.save_section_seq
                    BeatTracker.save_section_seq(results{6}, results{5}, results{1}, ...
                        fullfile(save_dir, [fname, '.section.seq']), obj.model.pattInfo.name);
                end
            end
        end
    end
    
    methods(Static)
        function [] = save_beats(beats, save_fln)
            fid = fopen(save_fln, 'w');
            fprintf(fid, '%.3f\t%i\n', beats(:,1:2)');
            fclose(fid);
        end
        
        function [] = save_downbeats(beats, save_fln)
            fid = fopen(save_fln, 'w');
            fprintf(fid, '%.3f\n', beats(beats(:, 2) == 1)');
            fclose(fid);
        end
        
        function [] = save_section_times(beats, save_fln)
            sectionTimes = beats(find(diff(beats(:,3)))+1,[1 3]);
            fid = fopen(save_fln, 'w');
            fprintf(fid, '%.3f\t%i\n', sectionTimes');
            fclose(fid);
        end
        
        function [] = save_tempo(tempo, save_fln)
            fid = fopen(save_fln, 'w');
            fprintf(fid, '%.2f\n', median(tempo));
            fclose(fid);
        end
        
        function [] = save_tempo_seq(ts, tempo, save_fln)
            dlmwrite(save_fln, [ts(:) tempo(:)], 'precision', '%.2f');
        end
        
        function [] = save_rhythm(rhythm, save_fln, rhythm_names)
            rhythm(rhythm < 0) = [];
            r = unique(rhythm);
            fid = fopen(save_fln, 'w');
            for i=1:length(r)
                fprintf(fid, '%s\n', rhythm_names{r(i)});
            end
            fclose(fid);
        end
        
        function [] = save_rhythm_seq(ts, rhythm, downbeats, save_fln, rhythm_names)
            ind = (downbeats(1:end-1) + downbeats(2:end))/2;
            fid = fopen(save_fln, 'w');
            for i=1:length(ind)
                [dummy ii] = min(abs(ts - ind(i)));
                r = rhythm(ii);
                fprintf(fid, '%s, %d\n', rhythm_names{r}, r);
            end
            fclose(fid);
        end
        
        function [] = save_section_names(sec, save_fln, sec_names)
            sec(sec < 0) = [];
            r = unique(sec);
            fid = fopen(save_fln, 'w');
            for i=1:length(r)
                fprintf(fid, '%s\n', sec_names{r(i)});
            end
            fclose(fid);
        end
        
        function [] = save_section_seq(ts, sec, beats, save_fln, sec_names)
            sectionTimes = beats(find(diff(beats(:,3)))+1,[1 3]);
            ind = (sectionTimes(1:end-1,1) + sectionTimes(2:end,1))/2;
            fid = fopen(save_fln, 'w');
            for i=1:length(ind)
                [dummy ii] = min(abs(ts - ind(i)));
                r = sec(ii);
                fprintf(fid, '%s, %d\n', sec_names{r}, r);
            end
            fclose(fid);
        end
        
        function [] = save_meter(meter, save_fln)
            meter(:,meter(1,:) == -1) = [];
            m = unique(meter', 'rows')';
            fid = fopen(save_fln, 'w');
            fprintf(fid, '%i/%i\n', m(1), m(2));
            fclose(fid);
        end
        
        
        function smoothedBeats = smooth_beats_sequence(inputBeats, win)
            % smooth_inputBeats(inputBeatFile, outputBeatFile, win)
            %   smooth beat sequence according to Dixon et al., Perceptual Smoothness
            %   of Tempo in Expressively Performed Music (2004)
            % ----------------------------------------------------------------------
            % INPUT Parameter:
            %   win                 :
            %
            % OUTPUT Parameter:
            %   Params            : structure array with beat tracking Paramseters
            %
            % 11.06.2012 by Florian Krebs
            % ----------------------------------------------------------------------
            if win < 1
                smoothedBeats = inputBeats;
                return
            end
            d = diff(inputBeats);
            % to correct for missing values at the ends, the sequence d is extended by
            % defining
            d = [d(1+win:-1:2); d; d(end-1:-1:end-win)];
            dSmooth = zeros(length(d), 1);
            for iBeat = 1+win:length(d)-win
                dSmooth(iBeat) = sum(d(iBeat-win:iBeat+win)) / (2*win+1);
            end
            dSmooth=dSmooth(win+1:end-win);
            smoothedBeats = inputBeats;
            smoothedBeats(2:end) = inputBeats(1) + cumsum(dSmooth);
            
        end
        
        function new_model = convert_to_new_model_format(old_model)
            new_model = old_model.convert_old_model_to_new();
        end
        
    end
    
end
