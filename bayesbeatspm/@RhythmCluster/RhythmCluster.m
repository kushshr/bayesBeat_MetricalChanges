classdef RhythmCluster < handle
    % WORKFLOW:
    %
    % if patterns are given via labels
    % 1) RhythmCluster
    % 2) make_cluster_assignment_file
    
    properties
        feature             % feature object
        feat_matrix_fln     % fln where feature values of dataset are stored
        clusters_fln        % fln where cluster ids are stored
        dataset             % training dataset on which clustering is performed
        train_lab_fln       % lab file with training data files
        train_file_list     % list of training files
        data_save_path      % path where cluster ids per bar are stored
        exclude_songs_fln   % vector of file ids that contain more than one bar and have supported meter
        n_clust_per_patt    % Number of clusters per pattern class (input parameter)
        n_cluster_tot       % Total number of pattern clusters
        pattern_size        % size of one rhythmical pattern {'beat', 'bar'}
        data_per_patt       % [nPatt x feat_dim*bar_grid]
        data_per_song       % [nSongs x feat_dim*bar_grid]
        patt2file           % [1 x nPatts] vector
        pattInfo            % Pattern info, including meter, section, length
        patt_2_cluster      % Pattern to cluster mapping, can take 1:n_cluster_tot
        patt_2_meter        % Pattern to meter mapping, can take 1:n_cluster_tot
        patt_2_section      % Pattern to section mapping
        patt_2_pattclass    % Pattern to pattern class mapping
        patt_2_len          % Pattern length in beats
        file_2_meter        % [nFiles x 2]
        rhythm2meter        % [n_clust_tot x 2]
        rhythm2meterID      % [n_clust_tot x 1] storing the mapping between cluster ID and meter ID
        rhythm2section      % [n_clust_tot x 1] storing the mapping between cluster ID and section ID
        rhythm2pattclass    % [n_clust_tot x 1] storing the mapping between cluster ID and combined pattern ID
        rhythm_names        % {n_clust_tot x 1} strings
        rhythm2len_num      % [n_clust_tot x 1] storing the mapping between cluster ID and section length numerator
        rhythm2len_den      % [n_clust_tot x 1] storing the mapping between cluster ID and section length denominator
        cluster_transition_matrix  % The transition matrix between clusters
        cluster_prior              % Cluster prior probability vector
        section_transition_matrix  % The transition matrix between sections
        section_prior              % Cluster prior probability vector
%         cluster_transitions_fln  % Filename to store cluster transition matrix
%         cluster_prior_fln        % Filename to store cluster prior vector
    end
    
    methods
        function obj = RhythmCluster(dataset, feat_type, frame_length, ...
                data_save_path, pattern_size)
            %  obj = RhythmCluster(dataset, feat_type, frame_length, data_save_path, pattern_size)
            %  Construct Rhythmcluster object
            % ----------------------------------------------------------------------
            % INPUT parameter:
            % dataset                 : path to lab file (list of training files)
            %
            % OUTPUT parameter:
            % obj
            %
            % 09.04.2013 by Florian Krebs
            % ----------------------------------------------------------------------
            if nargin == 1
                feat_type = {'lo230_superflux.mvavg', 'hi250_superflux.mvavg'};
                frame_length = 0.02;
                data_save_path = './data';
                pattern_size = 'bar';
            end
            obj.feature = Feature(feat_type, frame_length);
            if ~exist(data_save_path, 'dir')
                system(['mkdir ', data_save_path]);
            end
            obj.data_save_path = data_save_path;
            if exist('pattern_size', 'var')
                obj.pattern_size = pattern_size;
            else
                obj.pattern_size = 'bar';
            end
            % load list of training files
            if strfind(dataset, '.lab')
                obj.train_lab_fln = dataset;
                [~, obj.dataset, ~] = fileparts(dataset);
            else
                obj.train_lab_fln = fullfile('~/diss/data/beats/lab_files', [dataset, '.lab']);
                obj.dataset = dataset;
            end
            if exist(obj.train_lab_fln, 'file')
                train_data = Data(obj.train_lab_fln, 1);
                obj.train_file_list = train_data.file_list;
            else
                error('ERROR RhythmCluster.m: %s not found\n', obj.train_lab_fln);
            end
        end
        
        
        function make_feats_per_song(obj, whole_note_div)
            % make_feats_per_song: Computes mean of features for all bar
            % position within a song
            %
            % INPUT:    whole_note_div : number of bar positions for a whole
            %           note
            % OUTPUT:   obj.data_per_song [nSongs x feat_dim*bar_grid]
            % ------------------------------------------------------------
            
            % extract features per bar position and create obj.bar2file,
            % obj.data_per_bar, obj.file_2_meter, and obj.bar_2_meter
            obj.make_feats_per_bar(whole_note_div);
            % compute mean per song
            obj.data_per_song = NaN(length(obj.train_file_list), ...
                size(obj.data_per_bar, 2));
            for iCol=1:size(obj.data_per_bar, 2) % loop over (feat_dim*bar_grid)
                obj.data_per_song(:, iCol) = accumarray(obj.bar2file', ...
                    obj.data_per_bar(:, iCol), [], @mean);
            end
            % find songs that contain more than one bar and have
            % allowed meter
            exclude_song_ids = ~ismember(1:length(obj.train_file_list), ...
                unique(obj.bar2file));
            % save features (mean of all bars of one song)
            obj.feat_matrix_fln = fullfile(obj.data_save_path, ['onsetFeat-', ...
                num2str(obj.feature.feat_dim), 'd-', obj.dataset, '-songs.txt']);
            dlmwrite(obj.feat_matrix_fln, obj.data_per_song(~exclude_song_ids, :), ...
                'delimiter', '\t', 'precision', 4);
            fprintf('Saved data per song to %s\n', obj.feat_matrix_fln);
            % Save list of excluded files to file
            exclude_song_ids = find(exclude_song_ids);
            if ~isempty(exclude_song_ids)
                obj.exclude_songs_fln = fullfile(obj.data_save_path, [obj.dataset, ...
                    '-exclude.txt']);
                fid = fopen(obj.exclude_songs_fln, 'w');
                for i=1:length(exclude_song_ids)
                    fprintf(fid, '%s\n', obj.train_file_list{exclude_song_ids(i)});
                end
                fclose(fid);
                fprintf('Saved files to be excluded (%i) to %s\n', ...
                    length(exclude_song_ids), obj.exclude_songs_fln);
            end
        end
        
        
        function make_feats_per_patt(obj, whole_note_div)
            if exist(obj.train_lab_fln, 'file')
                fprintf('    Found %i files in %s\n', length(obj.train_file_list), ...
                    obj.train_lab_fln);
                dataPerPatt = [];
                for iDim =1:obj.feature.feat_dim
                    % organise features into bar position grid
                    Output = Data.extract_patts_from_feature(obj.train_file_list, ...
                        obj.feature.feat_type{iDim}, whole_note_div, ...
                        obj.feature.frame_length, obj.pattern_size, 1);
                    % compute mean feature of each bar position cell. 
                    dataPerPatt = [dataPerPatt, cellfun(@mean, Output.dataPerPatt)];
                end
                % save data to object
                obj.patt2file = Output.patt2file(:);
                obj.data_per_patt = dataPerPatt;
                obj.file_2_meter = Output.file2meter;
                obj.patt_2_meter = Output.file2meter(obj.patt2file, :);
                obj.patt_2_section = Output.patt2sec(:);  
                obj.patt_2_len = Output.patt_len(:);
                % save features organised by bars positions to file
                obj.feat_matrix_fln = fullfile(obj.data_save_path, ['onsetFeat_', ...
                    num2str(obj.feature.feat_dim), 'd_', obj.dataset, '.txt']);
                dlmwrite(obj.feat_matrix_fln, dataPerPatt, 'delimiter', '\t', 'precision', 4);
                fprintf('    Saved data per bar to %s\n', obj.feat_matrix_fln);
            else
                error('    %s not found', obj.train_lab_fln)
            end
        end
        
        function load_feats_per_bar(obj)
            obj.feat_matrix_fln = fullfile(obj.data_save_path, ['onsetFeat_', ...
                num2str(obj.feature.feat_dim), 'd_', obj.dataset, '.txt']);
            obj.data_per_bar = dlmread(obj.feat_matrix_fln, '\t');
        end
        
        function load_feats_per_song(obj)
            obj.feat_matrix_fln = fullfile(obj.data_save_path, ['onsetFeat-', ...
                num2str(obj.feature.feat_dim), 'd-', obj.dataset, '-songs.txt']);
            obj.data_per_song = dlmread(obj.feat_matrix_fln, '\t');
        end
        
        function [ca_fln, clust_trans_fln, clust_prior_fln] = do_clustering(obj, pattern_size, n_clust_per_patt, ...
                varargin)
            % Clusters the bars using a kmeans algorithm
            %
            % ------------------------------------------------------------------------
            % INPUT parameters:
            % n_clust_per_patt   number of clusters per pattern
            % ((pattern_scope))  removed the option; average for song does
            %                    is not menaingful when pattern lengths can
            %                    vary in a song
            % 'meters'           [num_meters x 2] meters that are expected
            %                    to be in the data
            % 'meter_names'      cell array [num_meters x 1]: assigns a name
            %                    for each meter. 
            % 'sections'         The [nMeters x 1] sections in the data, 
            %                    specified using the start beat of the 
            %                    section, for each meter
            % 'section_names'    The section names: [num_meter x 1] cell
            %                    array. Each cell is [nSections x 1]
            %                    cell array of names
            % 'save_pattern_fig' [default=1] save pattern plot
            % 'plotting_path'    [default='/tmp/']
            %
            % OUTPUT parameters:
            % ca_fln             file with cluster assignment data
            % clust_trans_fln
            %
            %
            % 06.09.2012 by Florian Krebs
            %
            % changelog:
            % 03.01.2013 modularize code into subfunctions
            % 29.01.2015 added config file name input parameter
            % ------------------------------------------------------------------------
            % parse arguments
            parser = inputParser;
            % set defaults
            % default_scope = 'song';
            % valid_scope = {'bar','song'};
            % check_scope = @(x) any(validatestring(x, valid_scope));
            default_plotting_path = '/tmp/';
            % add inputs
            addRequired(parser, 'obj', @isobject);
            addRequired(parser, 'n_clusters', @isnumeric);
            % addOptional(parser, 'pattern_scope', default_scope, check_scope);
            addParameter(parser, 'meter_names', '', @iscell);
            addParameter(parser, 'meters', -1, @isnumeric);
            addParameter(parser, 'section_names', '', @iscell);
            addParameter(parser, 'sections', '', @iscell);
            addParameter(parser, 'save_pattern_fig', 1, @isnumeric);
            addParameter(parser, 'plotting_path', default_plotting_path, @ischar);
            parse(parser, obj, n_clust_per_patt, varargin{:});
            % -------------------------------------------------------------
            obj.n_clust_per_patt = n_clust_per_patt;
            S = obj.data_per_patt;
            meter_per_item = obj.patt_2_meter;
            % if strcmpi(pattern_scope, 'bar') % Pattern scope now deprecated
            %     S = obj.data_per_bar;
            %    meter_per_item = obj.bar_2_meter;
            % else
            %     S = obj.data_per_song;
            %     meter_per_item = obj.file_2_meter;
            % end
            % check if meter found in the data corresponds to meter given
            % in the input arguments [num_meters x 2]
            meter_data = unique(meter_per_item, 'rows');
            if size(meter_data, 2) ~= 2
                meter_data = meter_data';
            end
            if parser.Results.meters == -1
                meters = meter_data;
                % create meter_names from time signature
                for i_meter=1:size(meters, 1)
                    meter_names{i_meter} = [num2str(meters(i_meter, 1)), ...
                        '-', num2str(meters(i_meter, 2))];
                end
            else
                same_num_meters = (size(meter_data, 1) == ...
                    size(parser.Results.meters, 1));
                same_content = ~sum(~ismember(meter_data, parser.Results.meters, ...
                    'rows'));
                if ~(same_num_meters && same_content)
                    error(['ERROR RhythmCluster.do_clustering: Number of ',...
                        'meters in data does not match number of meters', ...
                        'specified in the function input argument']);
                end
                meters = unique(meter_per_item, 'rows');
                meter_names = parser.Results.meter_names;
            end
            % Assuming sections have been sanely generated! No sanity check
            sec_per_item = obj.patt_2_section;
            if strcmp(pattern_size, 'section')
                sections = parser.Results.sections;
                section_names = parser.Results.section_names;
            elseif strcmp(pattern_size, 'bar')
                sections = cellfun(@(x) [x(1)], parser.Results.sections,'unif',0);
                section_names = cell(size(sections));
                section_names(:) = {{'full'}};
            end
            % Ajay/Andre: distribute patterns equally among meters and
            % sections
            % n_clusters_per_meter = ceil(n_clusters/size(meters,1))*ones(1,size(meters,1));
            
            % Some housekeeping to keep track of patterns, generate their
            % names and an ID for each
            [patt_class IA patt_id] = unique([meter_per_item sec_per_item], 'rows');
            obj.pattInfo.ID = 1:size(patt_class,1);
            for iPatt = 1:length(obj.pattInfo.ID)
                meterID = find(ismember(meters, patt_class(iPatt,1:2), 'rows'));
                sec = sections{meterID};
                secLen = diff([sec sec(1)+patt_class(iPatt,1)]);
                obj.pattInfo.n_clusters(iPatt) = n_clust_per_patt;
                obj.pattInfo.class(iPatt,1:3) = patt_class(iPatt,:);
                obj.pattInfo.name{iPatt} = [num2str(patt_class(iPatt,1)) '_' ...
                    num2str(patt_class(iPatt,2)) '_' num2str(patt_class(iPatt,3)) '_' ...
                    meter_names{meterID} '_' section_names{meterID}{patt_class(iPatt,3)}];
                obj.pattInfo.len(iPatt) = secLen(patt_class(iPatt,3));
                obj.pattInfo.class(iPatt,4) = sec(patt_class(iPatt,3));
                % Compute how many items we have per pattern
                n_items_per_patt(iPatt) = sum(patt_id == iPatt);
                % Update the patt_2_pattclass variable 
                obj.patt_2_pattclass(ismember(obj.patt_2_meter, patt_class(iPatt,1:2), 'rows') & ...
                    (patt_class(iPatt,3) == obj.patt_2_section)) = iPatt;
            end
            obj.n_cluster_tot = n_clust_per_patt * size(patt_class,1);
            % normalise feature to treat them all equally when clustering
            bar_pos = size(S, 2) / obj.feature.feat_dim;
            for i_dim = 1:obj.feature.feat_dim
                vals = S(:, (i_dim-1)*bar_pos+1:bar_pos*i_dim);
                vals = vals - nanmean(vals(:));
                vals = vals / nanstd(vals(:));
                S(:, (i_dim-1)*bar_pos+1:bar_pos*i_dim) = vals;
            end
            opts = statset('MaxIter', 200);
            % cluster different meter separately
            ctrs = zeros(obj.n_cluster_tot, size(S, 2));
            cidx = zeros(size(S, 1), 1);
            p = 1;
            % replace nans because the MATLAB kmeans ignores all datapoints
            % that contain nans
            S(isnan(S)) = -9999;
            obj.rhythm2meter = zeros(obj.n_cluster_tot, 2);
            obj.rhythm2meterID = zeros(obj.n_cluster_tot, 1);
            obj.rhythm2section = zeros(obj.n_cluster_tot, 1);
            obj.rhythm2pattclass = zeros(obj.n_cluster_tot, 1);
            obj.rhythm2len_num = zeros(obj.n_cluster_tot, 1);
            obj.rhythm2len_den = zeros(obj.n_cluster_tot, 1);
            
            for iPatt = 1:length(obj.pattInfo.ID)
                idx_i = (patt_id == iPatt);
                if obj.pattInfo.n_clusters(iPatt) > 1
                    [cidx(idx_i), ctrs(p:p+obj.pattInfo.n_clusters(iPatt)-1, :)] = ...
                        kmeans(S(idx_i, :), obj.pattInfo.n_clusters(iPatt), 'Distance', ...
                        'sqEuclidean', 'Replicates', 5, 'Options', opts);
                    cidx(idx_i) = cidx(idx_i) + p - 1;
                    
                    obj.rhythm2meter(p:p+obj.pattInfo.n_clusters(iPatt)-1, :) ...
                        = repmat(obj.pattInfo.class(iPatt, 1:2), ...
                           obj.pattInfo.n_clusters(iPatt), 1);
                    obj.rhythm2meterID(p:p+obj.pattInfo.n_clusters(iPatt)-1, 1) ...
                        = repmat(find(ismember(meters, obj.pattInfo.class(iPatt, 1:2), 'rows')), ...
                        obj.pattInfo.n_clusters(iPatt), 1);
                    obj.rhythm2section(p:p+obj.pattInfo.n_clusters(iPatt)-1, 1) ...
                        = repmat(obj.pattInfo.class(iPatt, 3), ...
                        obj.pattInfo.n_clusters(iPatt), 1);
                    obj.rhythm2pattclass(p:p+obj.pattInfo.n_clusters(iPatt)-1, 1) = iPatt;
                    obj.rhythm2len_den(p:p+obj.pattInfo.n_clusters(iPatt)-1, 1) ...
                        = repmat(obj.pattInfo.class(iPatt, 2), ...
                        obj.pattInfo.n_clusters(iPatt), 1);
                    obj.rhythm2len_num(p:p+obj.pattInfo.n_clusters(iPatt)-1, 1) = obj.pattInfo.len(iPatt);
                    p = p + obj.pattInfo.n_clusters(iPatt);
                else
                    ctrs(p, :) = mean(S(idx_i, :));
                    cidx(idx_i) = p;
                    obj.rhythm2meter(p,:) = obj.pattInfo.class(iPatt, 1:2);
                    obj.rhythm2meterID(p) = find(ismember(meters, obj.pattInfo.class(iPatt, 1:2), 'rows'));
                    obj.rhythm2section(p) = obj.pattInfo.class(iPatt, 3);
                    obj.rhythm2pattclass(p) = iPatt;
                    obj.rhythm2len_num(p) = obj.pattInfo.len(iPatt);
                    obj.rhythm2len_den(p) = obj.pattInfo.class(iPatt, 2);
                    p=p+1;
                end
            end
            
            
            % reintroduce nans
            ctrs(ctrs==-9999) = nan;
            if parser.Results.save_pattern_fig
                % plot patterns and save plot to png file
                obj.plot_patterns(cidx, ctrs, bar_pos, ...
                    parser.Results.plotting_path);
            end
            % save cluster assignments
            obj.clusters_fln = fullfile(obj.data_save_path, ['ca-', ...
                obj.dataset, '-', num2str(obj.feature.feat_dim), ...
               'd-', num2str(obj.n_cluster_tot), 'R-kmeans.mat']);
            % make up a name for each newly created rhythm pattern. It should
            % contain both the meter, section and pattern ID
            for i_R = 1:obj.n_cluster_tot
                obj.rhythm_names{i_R} = [obj.dataset, '-', ...
                    obj.pattInfo.name{obj.rhythm2pattclass(i_R)}, '_', ...
                    num2str(mod(i_R-1, obj.n_clust_per_patt)+1)];
            end
            obj.patt_2_cluster = cidx;
            % create rhythm transitions
            obj.compute_cluster_transitions();
            % save all variables related to a cluster assignment to file
            obj.save_cluster_alignment_file();
            ca_fln = obj.clusters_fln;
        end
        
        function [] = compute_cluster_transitions(obj)
            % To find section transitions
            B = zeros(length(obj.pattInfo.ID));
            B_prior = zeros(1,length(obj.pattInfo.ID));
            for iFile=1:length(obj.train_file_list)
                patts = find(obj.patt2file==iFile);
                for iPatt=1:length(patts(1:end-1))
                    B(obj.patt_2_pattclass(patts(iPatt)), obj.patt_2_pattclass(patts(iPatt+1))) = ...
                        B(obj.patt_2_pattclass(patts(iPatt)), obj.patt_2_pattclass(patts(iPatt+1))) + 1;
                    B_prior(obj.patt_2_pattclass(patts(iPatt))) = B_prior(obj.patt_2_pattclass(patts(iPatt))) + 1;
                end
            end
            B_prior = B_prior ./ sum(B_prior);
            B = bsxfun(@rdivide, B , sum(B , 2));
            % To find cluster transitions 
            A = zeros(obj.n_cluster_tot, obj.n_cluster_tot);
            A_prior = zeros(1,obj.n_cluster_tot);
            for iFile=1:length(obj.train_file_list)
                patts = find(obj.patt2file==iFile);
                for iPatt=1:length(patts(1:end-1))
                    A(obj.patt_2_cluster(patts(iPatt)), obj.patt_2_cluster(patts(iPatt+1))) = ...
                        A(obj.patt_2_cluster(patts(iPatt)), obj.patt_2_cluster(patts(iPatt+1))) + 1;
                    A_prior(obj.patt_2_cluster(patts(iPatt))) = A_prior(obj.patt_2_cluster(patts(iPatt))) + 1;
                end
            end
            A_prior = A_prior ./ sum(A_prior);
            % normalise transition matrix
            A = bsxfun(@rdivide, A , sum(A , 2));
            % Now need to regularize A and prior if needed
            if ~all(A_prior)
                A_prior = A_prior + 0.01;
                A_prior = A_prior ./ sum(A_prior);
                fprintf('    Pattern prior vector has zero counts, back off model used...\n');
            end
            if any(isnan(sum(A,2)))
                indxA = find(isnan(sum(A,2)));
                A(indxA,:) = 0;
                for pp = 1:length(indxA)
                    A(indxA(pp),indxA(pp)) = 1;
                end
                fprintf('    Pattern transition matrix rank deficient, and regularized...\n');
            end
            % obj.cluster_transition_matrix = eye(size(A));
            obj.cluster_transition_matrix = [0.98 0.02 0.02; 0.02 0.98 0.02; 0.02 0.02 0.98];
            obj.cluster_prior = A_prior;
            obj.section_transition_matrix = [0.98 0.02 0.02; 0.02 0.98 0.02; 0.02 0.02 0.98];
            obj.section_prior = B_prior;
        end
        
        function [ca_fln] = make_cluster_assignment_file(obj, clusterType, rhythm_names)
            % [bar2rhythm] = make_cluster_assignment_file(trainLab)
            %   Creates vector bar2rhythm that assigns each
            %   bar in trainLab to the pattern specified by the dancestyle annotation.
            % ----------------------------------------------------------------------
            %INPUT parameter:
            % trainLab          : filename of labfile (e.g., 'boeck.lab' or 'ballroom.lab').
            %                       a labfile is a textfile with paths to files that are analyzed
            % clusterType         : {'meter', 'dancestyle', 'rhythm', 'none'}
            %                   'meter': bars are clustered according to the meter (functions reads .meter file);
            %                   'dancestyle', according to the genre (functions reads .dancestyle file))
            %                   'none' : put all bars into one single
            %                   cluster
            % clustAssFln       : if clusterType='auto', specify a textfile, with a
            %                       pattern id for each file (in the same order as in trainLab)
            % rhythm_names      : cell array of strings
            %
            %OUTPUT parameter:
            % bar2rhythm        : [nBars x 1] assigns each bar to a pattern
            % file_2_cluster    : [nFiles x 1] assigns each file to a pattern
            % n_clusters        : number of patterns
            % bar2file          : [nBars x 1] assigns each bar a file id
            % file2nBars        : [nFiles x 1] number of bars per file
            % rhythm_names      : {R x 1} string
            % rhythm2meter      : [R x 2] meter for each pattern
            % rhythm2meterID      : [R x 1] meterID for each pattern
            % 09.07.2013 by Florian Krebs
            % ----------------------------------------------------------------------
            if nargin == 1
                clusterType = 'auto';
            end
            if strcmpi(clusterType, 'auto')
                songClusterIds = load(obj.clusters_fln, '-ascii');
            end
            bar2rhythm = [];
            obj.bar2file = [];
            %             dancestyles = {'ChaCha', 'Jive', 'Quickstep', 'Rumba', 'Samba', 'Tango', ...
            %               'VienneseWaltz', 'Waltz'};
            % check if there are songs that should be excluded
            if isempty(obj.exclude_songs_fln)
                ok_songs = 1:length(obj.train_file_list);
            else
                fid = fopen(obj.exclude_songs_fln, 'r');
                exclude_songs = textscan(fid, '%s');
                fclose(fid);
                exclude_songs = exclude_songs{1};
                ok_songs = find(~ismember(obj.train_file_list, exclude_songs));
            end
            meter = zeros(length(ok_songs), 2);
            meters = []; % unique meters present in the data [T x 2]
            fileCounter = 0;
            nBars = zeros(length(ok_songs), 1);
            obj.file_2_cluster = zeros(length(ok_songs), 1);
            for iFile = 1:length(ok_songs)
                [~, fname, ~] = fileparts(obj.train_file_list{ok_songs(iFile)});
                fprintf('%i) %s\n', iFile, fname);
                [beats, error] = Data.load_annotations_bt(...
                    obj.train_file_list{ok_songs(iFile)}, 'beats');
                if error, error('no beat file found\n'); end
                if size(beats, 2) < 2
                    error('Downbeat annotations missing for %s\n', ...
                        obj.train_file_list{ok_songs(iFile)});
                end
                % load meter from file
                meter(fileCounter + 1, :) = Data.load_annotations_bt(...
                    obj.train_file_list{ok_songs(iFile)}, 'meter');
                % get pattern id of file
                switch lower(clusterType)
                    case 'meter'
                        if isempty(meters)
                            meter_present = 0;
                        else
                            meter_present = ismember(meters, ...
                                meter(fileCounter + 1, :), 'rows');
                        end
                        if sum(meter_present) == 0
                            meters = [meters; meter(fileCounter + 1, :)];
                        end
                        patternId = find(ismember(meters, ...
                            meter(fileCounter + 1, :), 'rows'));
                    case 'dancestyle'
                        style = Data.load_annotations_bt(...
                            obj.train_file_list{ok_songs(iFile)}, 'dancestyle');
                        patternId = find(strcmp(rhythm_names, style));
                        if isempty(patternId)
                            fprintf('Please add %s to the rhythm_names\n', style);
                        end
                        obj.file_2_cluster(iFile) = patternId;
                    case 'auto'
                        patternId = songClusterIds(iFile);
                    case 'rhythm'
                        style = Data.load_annotations_bt(...
                            strrep(obj.train_file_list{ok_songs(iFile)}, ...
                            'audio', 'annotations/rhythm'), 'rhythm');
                        patternId = find(strcmp(rhythm_names, style));
                    case 'none'
                        patternId = 1;
                end
                if strcmp(obj.pattern_size, 'bar')
                    obj.rhythm2meter(patternId, :) = ...
                        meter(fileCounter + 1, :);
                    obj.rhythm2meterID(patternId, :) = patternId;
                elseif strcmp(obj.pattern_size, 'beat')
                    obj.rhythm2meter(patternId, :) = [1, 4];
                    obj.rhythm2meterID(patternId, :) = 1;
                else
                    error('Pattern size unknown!')
                end
                fileCounter = fileCounter + 1;
                % determine number of bars
                if strcmp(obj.pattern_size, 'beat')
                    nBars(iFile) = size(beats, 1) - 1;
                else
                    [nBars(iFile), ~] = Data.get_full_bars(beats);
                end
                bar2rhythm = [bar2rhythm; ones(nBars(iFile), 1) * patternId];
                obj.bar2file = [obj.bar2file; ones(nBars(iFile), 1) * ok_songs(iFile)];
            end
            obj.n_clusters = max(bar2rhythm);
            % BUG HERE?: The meters have to be ordered in the increasing
            % order, fails otherwise! Check!
            obj.bar_2_cluster = bar2rhythm;
            obj.clusters_fln = fullfile(obj.data_save_path, ['ca-', obj.dataset, ...
                '-', num2str(obj.feature.feat_dim), 'd-', num2str(obj.n_clusters), ...
                'R-', clusterType, '.mat']);
            ca_fln = obj.clusters_fln;
            if ~exist('rhythm_names', 'var')
                for i = unique(bar2rhythm(:))'
                    rhythm_names{i} = [clusterType, num2str(i)];
                end
            end
            obj.rhythm_names = rhythm_names;
            obj.file_2_nBars = nBars;
            obj.compute_cluster_transitions();
            obj.save_cluster_alignment_file();
        end
    end
    
    methods (Access=protected)
        function [] = save_cluster_alignment_file(obj)
            rhythm2meter = obj.rhythm2meter;
            rhythm2meterID = obj.rhythm2meterID;
            rhythm2section = obj.rhythm2section;
            rhythm2pattclass = obj.rhythm2pattclass;
            patt2file = obj.patt2file;
            patt2rhythm = obj.patt_2_cluster;
            patt2sec = obj.patt_2_section;
            patt2meter = obj.patt_2_meter;
            patt2pattclass = obj.patt_2_pattclass;
            patt2len = obj.patt_2_len;
            rhythm_names = obj.rhythm_names;
            pattInfo = obj.pattInfo;
            pr = obj.cluster_transition_matrix;
            prprior = obj.cluster_prior;
            B_pr = obj.section_transition_matrix;
            B_prprior = obj.section_prior;
            rhythm2len_num = obj.rhythm2len_num;
            rhythm2len_den = obj.rhythm2len_den;
            save(obj.clusters_fln, '-v7.3', ...
                'rhythm2meter', 'rhythm2meterID', 'rhythm2section', 'rhythm2pattclass',...
                'rhythm2len_num', 'rhythm2len_den', 'patt2file', 'patt2rhythm', ...
                'patt2sec', 'patt2meter', 'patt2pattclass', 'patt2len', 'rhythm_names', ...
                'pattInfo', 'pr', 'prprior', 'B_pr', 'B_prprior');
            fprintf('    Saved pattern variables to %s\n', obj.clusters_fln);
        end
        
        function [] = plot_patterns(obj, cidx, ctrs, bar_pos, ...
                plotting_path)
            plot_cols = ceil(sqrt(obj.n_cluster_tot));
            h = figure( 'Visible','off');
            set(h, 'Position', [100 100 obj.n_cluster_tot*100 obj.n_cluster_tot*100]);
            items_per_cluster = hist(cidx, obj.n_cluster_tot);
            col = hsv(obj.feature.feat_dim);
            for c = 1:obj.n_cluster_tot
                subplot(ceil(obj.n_cluster_tot/plot_cols), plot_cols, c)
                hold on
                for fdim = 1:obj.feature.feat_dim
                    data = ctrs(c, (fdim-1)*bar_pos+1:fdim*bar_pos);
                    data = data - min(data);
                    data = data / max(data);
                    data = data + fdim;
                    stairs([data, data(end)], 'Color', col(fdim, :));
                end
                if obj.feature.feat_dim == 1
                    y_label = obj.feature.feat_type{1};
                else
                    y_label = sprintf('Bottom: %s', ...
                        strrep(obj.feature.feat_type{1}, '_', '\_'));
                end
                ylabel(sprintf('%s', y_label));
                xlabel('bar position')
                title(sprintf('cluster %i (%i examples)', c, ...
                    items_per_cluster(c)));
                xlim([1 length(data)])
            end
            outfile = fullfile(plotting_path, ['patterns-', ...
                obj.dataset, '-kmeans-', ...
                num2str(obj.n_cluster_tot), '.png']);
            fprintf('    Writing patterns to %s\n', outfile);
            % save to png
            print(h, outfile, '-dpng');
        end
    end
end

