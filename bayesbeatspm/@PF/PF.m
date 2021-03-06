classdef PF < handle
    % Hidden Markov Model Class
    properties (SetAccess=private)
        M                           % number of positions in a 4/4 bar
        Meff                        % number of positions per rhythm [R x 1]
        N                           % number of tempo states
        R                           % number of rhythmic pattern states
        V                           % number of section pattern states
        T                           % number of meters
        pr                          % probability of a switch in rhythmic pattern
        prprior                     % prior probability of a pattern
        B_pr                        % probability of a switch in a section
        B_prprior                   % prior probability of a section
        frames_per_beat             % cell arrays with tempi relative to the framerate
        % frames_per_beat{1} is a row vector with tempo
        % values in [audio frames per beat] for pattern 1
        rhythm2meter                % assigns each rhythmic pattern to a meter [R x 1]
        rhythm2meterID              % assigns each rhythmic pattern to a meter state (better if present)
        rhythm2pattclass            % assigns each rhythmic pattern to a pattern class 
        rhythm2section              % assigns each rhythmic pattern to a section
        rhythm2len_num              % assigns each rhythmic pattern to the length (numerator)
        rhythm2len_den              % assigns each rhythmic pattern to the length (denominator)
        pattInfo                    % Pattern info structure
        rhythm2meter_state          % assigns each rhythmic pattern to a
        % meter state (1, 2, ...)  - var not needed anymore but keep due to
        % compatibility
        meter_state2meter           % specifies meter for each meter state
        % (9/8, 8/8, 4/4) [2 x nMeters] - var not needed anymore but keep due to
        % compatibility
        barGrid                     % number of different observation model params per bar (e.g., 64)
        minN                        % min tempo (n_min) for each rhythmic pattern
        maxN                        % max tempo (n_max) for each rhythmic pattern
        frame_length                % audio frame length in [sec]
        dist_type                   % type of parametric distribution
        sample_trans_fun            % transition model
        disc_trans_mat              % transition matrices for discrete states
        obs_model                   % observation model
        initial_m                   % initial value of m for each particle
        initial_n                   % initial value of n for each particle
        initial_r                   % initial value of r for each particle
        initial_v                   % initial value of v for each particle
        nParticles                  % number of particles
        particles                   % A structure to hold the particle filter outputs
        particles_grid              % particle grid for viterbi
        sigma_N                     % standard deviation of tempo transition
        bin2decVector               % vector to compute indices for disc_trans_mat quickly
        ratio_Neff                  % reample if NESS < ratio_Neff * nParticles
        resampling_interval         % fixed resampling interval [samples]
        resampling_scheme           % type of resampling employed
        warp_fun                    % function that warps the weights w to a compressed space
        save_inference_data         % save intermediate output of particle filter for visualisation
        state_distance_coefficients % distance factors for k-means clustering [1 x nDims]
        cluster_merging_thr         % if distance < thr: merge
        cluster_splitting_thr       % if spread > thr: split
        inferenceMethod             % 'PF'
        n_max_clusters              % If number of clusters > n_max_clusters, kill cluster with lowest weight
        n_initial_clusters          % Number of cluster to start with
        rhythm_names                % cell array of rhythmic pattern names
        train_dataset                % dataset, which PF was trained on
        pattern_size                % size of one rhythmical pattern {'beat', 'bar'}
        use_silence_state           % state that describes non-musical content
        n_depends_on_r              % Flag to say if tempo depends on the style/pattern/rhythm
        patt_trans_opt              % Flag to say how pattern transitions need to be done
        pattWtMatInit               % Pattern weight matrix, init
        currBarStartInit            % A pointer to the frame the current bar started, init
        peakInfMode                 % If inference is to be done only at peaks
        peakParam                   % Peak finder params
        tempoInfMode                % tempo informed inference 
        storevar                    % A place to backup initial values
    end
    
    methods
        function obj = PF(Params, data)
            obj.M = Params.M;
            obj.N = Params.N;
            obj.R = Params.R * length(data.pattInfo.ID);
            obj.V = length(data.pattInfo.ID);
            obj.T = size(Params.meters, 2);
            obj.nParticles = Params.nParticles;
            %             obj.sigma_N = Params.sigmaN; % moved this to
            %             make_transition_model
            obj.barGrid = max(Params.whole_note_div * data.rhythm2len_num ...
                ./ data.rhythm2len_den);
            obj.frame_length = Params.frame_length;
            obj.dist_type = Params.observationModelType;
            obj.rhythm2meter = data.rhythm2meter;
            obj.rhythm2meterID = data.rhythm2meterID;
            obj.rhythm2pattclass = data.rhythm2pattclass;
            obj.rhythm2section = data.rhythm2section;
            obj.rhythm2len_num = data.rhythm2len_num;
            obj.rhythm2len_den = data.rhythm2len_den;
            obj.pr = data.pr;
            obj.prprior = data.prprior;
            obj.B_pr = data.B_pr;
            obj.B_prprior = data.B_prprior;
            obj.pattInfo = data.pattInfo;
            patt_durations = data.rhythm2len_num ./ data.rhythm2len_den;
            obj.Meff = round((patt_durations) ...
                * (Params.M ./ (max(patt_durations)))); 
            obj.Meff = obj.Meff(:);
            obj.ratio_Neff = Params.ratio_Neff;
            obj.resampling_scheme = Params.resampling_scheme;
            obj.warp_fun = Params.warp_fun;
            obj.save_inference_data = Params.save_inference_data;
            obj.state_distance_coefficients = Params.state_distance_coefficients;
            obj.cluster_merging_thr = Params.cluster_merging_thr;
            obj.cluster_splitting_thr = Params.cluster_splitting_thr;
            obj.n_max_clusters = Params.n_max_clusters;
            obj.n_initial_clusters = Params.n_initial_clusters;
            obj.rhythm_names = data.rhythm_names;
            obj.pattern_size = Params.pattern_size;
            obj.resampling_interval = Params.res_int;
            obj.use_silence_state = Params.use_silence_state;
            obj.n_depends_on_r = Params.n_depends_on_r;
            obj.patt_trans_opt = Params.patt_trans_opt;
            obj.peakInfMode = Params.peakInfMode;
            obj.tempoInfMode = Params.tempoInfMode;
            obj.peakParam = Params.peak;
            tempVar = ver('matlab');
            if str2double(tempVar.Release(3:6)) < 2011
                % Matlab 2010
                disp('MATLAB 2010');
                RandStream.setDefaultStream(RandStream('mt19937ar','seed', ...
                    sum(100*clock)));
            else
                rng('shuffle');
            end
        end
        
        function obj = make_initial_distribution(obj)
            
            % TODO: Implement prior initial distribution for tempo
            obj.initial_m = zeros(obj.nParticles, 1);
            obj.initial_n = zeros(obj.nParticles, 1);
            obj.initial_r = zeros(obj.nParticles, 1);
            obj.initial_v = zeros(obj.nParticles, 1);
            
            % use pseudo random monte carlo
            % n_grid = min(obj.minN):max(obj.maxN);
            % n_m_cells = floor(obj.nParticles / length(n_grid));
            m_grid_size_init = min(obj.Meff ./ obj.rhythm2len_num)/8;  % 8th of a beat resolution for initial grid
            m_grid_points = ceil(sum(obj.Meff) / m_grid_size_init);
            n_grid_points = floor(obj.nParticles / m_grid_points);
            n_inv_grid = linspace(1/max(obj.maxN), 1/min(obj.minN), n_grid_points);
            
            n_m_cells = floor(obj.nParticles / length(n_inv_grid));
            m_grid_size = sum(obj.Meff) / n_m_cells;
            r_m = rand(obj.nParticles, 1) - 0.5; % between -0.5 and +0.5
            r_n = rand(obj.nParticles, 1) - 0.5;
            nParts = zeros(obj.R, 1);
            c=1;
            for iR = 1:obj.R
                % create positions between 1 and obj.Meff(obj.rhythm2meter_state(iR))
                m_grid = (1+m_grid_size/2):m_grid_size:...
                    (obj.Meff(iR)-m_grid_size/2);
                % n_grid_iR = obj.minN(iR):obj.maxN(iR);
                n_inv_grid_iR = linspace(1/obj.maxN(iR), 1/obj.minN(iR), n_grid_points);
                n_inv_grid_unit_iR = diff(n_inv_grid_iR(1:2));
                nParts(iR) = length(m_grid) * length(n_inv_grid_iR);
                temp = repmat(m_grid, length(n_inv_grid_iR), 1);
                obj.initial_m(c:c+nParts(iR)-1) = temp(:)+ ...
                    r_m(c:c+nParts(iR)-1) * m_grid_size;
                % Randomize tempo values based on the value of tempo
                temp = repmat(n_inv_grid_iR, 1, length(m_grid))' + ...
                    (n_inv_grid_unit_iR * r_n(c:c+nParts(iR)-1)); 
                temp = flipud(1./temp);
                % redistribute points outside the allowed ranges
                temp(temp > obj.maxN(iR)) = rand(sum(temp > obj.maxN(iR)),1) * ...
                    (obj.maxN(iR) - obj.minN(iR)) + obj.minN(iR);
                temp(temp < obj.minN(iR)) = rand(sum(temp < obj.minN(iR)),1) * ...
                    (obj.maxN(iR) - obj.minN(iR)) + obj.minN(iR);
                obj.initial_n(c:c+nParts(iR)-1) = temp;                
                obj.initial_r(c:c+nParts(iR)-1) = iR;
                c = c + nParts(iR);
            end
            if sum(nParts) < obj.nParticles
                % add remaining particles randomly
                % random pattern assignment
                obj.initial_r(c:end) = round(rand(obj.nParticles+1-c, 1)) ...
                    * (obj.R-1) + 1;
                n_between_0_and_1 = (r_n(c:end) + 0.5);
                m_between_0_and_1 = (r_m(c:end) + 0.5);
                % map n_between_0_and_1 to allowed tempo range
                obj.initial_n(c:end) = n_between_0_and_1 .* obj.maxN(...
                    obj.initial_r(c:end)) + (1 - n_between_0_and_1) .* ...
                    obj.minN(obj.initial_r(c:end));
                % map m_between_0_and_1 to allowed position range
                obj.initial_m(c:end) = m_between_0_and_1 .* ...
                    (obj.Meff(obj.initial_r(c:end))-1) + 1;
            end
            obj.initial_v = obj.rhythm2pattclass(obj.initial_r);
            obj.storevar.initial_m = obj.initial_m;
            obj.storevar.initial_r = obj.initial_r;
            obj.storevar.initial_n = obj.initial_n;
            obj.storevar.initial_v = obj.initial_v;
            if (obj.patt_trans_opt == 2 || obj.patt_trans_opt == 3)
                obj.pattWtMatInit = ones(obj.nParticles,1) * obj.prprior;  % Use from prior specified
                obj.currBarStartInit = ones(obj.nParticles,1);  % The curr bar started just now
            end
        end
        
        function obj = make_transition_model(obj, minTempo, maxTempo, ...
                params, data)
            position_states_per_beat = obj.Meff ./ obj.rhythm2len_num(:, 1);
            obj.sigma_N = params.sigmaN;
            % save pattern change probability and save as matrix [RxR]
            pr = data.pr;
            if (length(pr(:)) == 1) && (obj.R > 1)
                % expand pr to a matrix [R x R]
                % transitions to other patterns
                pr_mat = ones(obj.R, obj.R) * (pr / (obj.R-1));
                % pattern self transitions
                pr_mat(logical(eye(obj.R))) = (1-pr);
                obj.pr = pr_mat;
            elseif (size(pr, 1) == obj.R) && (size(pr, 2) == obj.R)
                % ok, do nothing
                obj.pr = pr;
            else
                error('p_r has wrong dimensions!\n');
            end
            
            if obj.patt_trans_opt == 0 && strcmp(obj.pattern_size, 'bar')     % No transitions allowed, so make pr a block identity matrix
                % TODO here
                obj.pr = eye(size(obj.pr));
            end
            
            obj.prprior = data.prprior;
            % convert from BPM into barpositions / audio frame
            % obj.minN = floor(position_states_per_beat .* obj.frame_length .* minTempo ./ 60);
            obj.minN = 0.8*position_states_per_beat .* obj.frame_length .* minTempo ./ 60; % 20% allowance on lower side
            obj.minN = obj.minN(:);
            % obj.minN(obj.minN > 1) = floor(obj.minN(obj.minN > 1)); % Floor for minN > 1
            % obj.minN(obj.minN < 1) = obj.minN(obj.minN < 1)*0.8; % A small 20% allowance for minN < 1
            obj.minN(obj.minN < 0.2) = 0.2;  % THIS IS THE MINIMUM TEMPO POSSIBLE, approx 6 bpm. 
            % obj.maxN = ceil(position_states_per_beat .* obj.frame_length .* maxTempo ./ 60);
            obj.maxN = 1.2*position_states_per_beat .* obj.frame_length .* maxTempo ./ 60; % 20% allowance on the higher side
            obj.maxN = obj.maxN(:);
            if max(obj.maxN) ~= obj.N
                fprintf('    N should be %i instead of %i -> corrected\n', ...
                    max(obj.maxN), obj.N);
                obj.N = max(obj.maxN);
            end
            % fprintf('    minN = %s\n', sprintf('%.2f ',minN))
            % fprintf('    maxN = %s\n', sprintf('%.2f ',maxN))

            % if ~n_depends_on_r % no dependency between n and r
            %     obj.minN = ones(1, obj.R) * min(obj.minN);
            %     obj.maxN = ones(1, obj.R) * max(obj.maxN);
            %     obj.N = max(obj.maxN);
            % end
            obj.frames_per_beat = cell(obj.R, 1);
            for ri = 1:obj.R
                obj.frames_per_beat{ri} = position_states_per_beat(ri) ./ ...
                    ([obj.minN(ri) obj.maxN(ri)]);
            end
            for r_i = 1:obj.R
                bpms = 60 ./ (obj.frames_per_beat{r_i} * obj.frame_length);
                fprintf('    R=%i: Tempo limited to %.1f - %.1f bpm, with minN = %.2f and maxN = %.2f\n', ...
                    r_i, bpms(1), bpms(end), obj.minN(r_i),obj.maxN(r_i));
            end
            % Transition function to propagate particles
            obj.sample_trans_fun = @(x) obj.propagate_particles_pf(obj, x);
            % Just a backup
            obj.storevar.minN = obj.minN;
            obj.storevar.maxN = obj.maxN;
            obj.storevar.N = obj.N;
        end
        
        function obj = make_observation_model(obj, train_data)
            
            % Create observation model
%             obj.obs_model = ObservationModel(obj.dist_type, obj.rhythm2meter, ...
%                 obj.M, obj.N, obj.R, obj.barGrid, obj.Meff, ...
%                 train_data.feat_type, obj.use_silence_state);
            obj.obs_model = ObservationModel(obj.dist_type, train_data, ...
                 obj.M, obj.N, obj.R, obj.barGrid, obj.Meff, obj.use_silence_state);
            
            % Train model
            obj.obs_model = obj.obs_model.train_model(train_data);
            obj.train_dataset = train_data.dataset;
        end
        
        function obj = retrain_observation_model(obj, data_file_pattern_barpos_dim, pattern_id)
            %{
            Retrains the observation model for all states corresponding to <pattern_id>.

            :param data_file_pattern_barpos_dim: training data
            :param pattern_id: pattern to be retrained. can be a vector, too.
            :returns: the retrained hmm object

            %}
            obj.obs_model = obj.obs_model.retrain_model(data_file_pattern_barpos_dim, pattern_id);
            %             obj.obs_model.learned_params{pattern_id, :} = ...
            %                 obj.obs_model.fit_distribution(data_file_pattern_barpos_dim(:, pattern_id, :, :));
            
        end
        
        function results = do_inference(obj, y, fname, inference_method, varargin)
            do_output = 1;
            gtTempo.bpm = [];
            gtTempo.range = [];
            gtSama = [];
            if length(varargin) == 3
                do_output = varargin{1};
                gtTempo.bpm = varargin{2}.val;
                gtTempo.range = varargin{2}.range;
                gtSama = varargin{3};
            elseif length(varargin) == 2
                do_output = varargin{1};
                gtTempo.bpm = varargin{2}.val;
                gtTempo.range = varargin{2}.range;
            elseif length(varargin) == 1
                do_output = varargin{1};
            end
            if isempty(strfind(inference_method, 'PF'))
                error('Inference method %s not compatible with PF model\n', inference_method);
            end
            % Initialize differently if doing a tempo informed or sama
            % informed tracking
            % Default case, do not
            startOffset = 1;      % No offset in frames
            startOffsetTime = 0;  % No offset in time
            % Copy from backup
            obj.initial_n = obj.storevar.initial_n;
            obj.initial_m = obj.storevar.initial_m;
            obj.initial_r = obj.storevar.initial_r;
            obj.initial_v = obj.storevar.initial_v;
            obj.minN = obj.storevar.minN;
            obj.maxN = obj.storevar.maxN;
            obj.N = obj.storevar.N;
            % Tempo informed tracking
            if sum(ismember([1,2],obj.tempoInfMode)) && ~isempty(gtTempo.bpm)
                % Now assign new values for tempo 
                position_states_per_beat = obj.Meff ./ obj.rhythm2len_num(:, 1);
                medTempo = position_states_per_beat .* obj.frame_length .* ...
                    gtTempo.bpm ./ 60; 
                obj.minN = (1-gtTempo.range) .* medTempo; % allowance on lower side
                obj.minN = obj.minN(:);
                obj.minN(obj.minN < 0.2) = 0.2;  % THIS IS THE MINIMUM TEMPO POSSIBLE, approx 6 bpm. 
                obj.maxN = (1+gtTempo.range) .* medTempo; % allowance on higher side
                obj.maxN = obj.maxN(:);
                obj.N = max(obj.maxN);
                for r_i = 1:obj.R
                    fprintf('    Tempo informed tracking, narrowing tempo range for ');
                    fprintf('    R=%i: Tempo limited to minN = %.2f and maxN = %.2f\n', ...
                    r_i, obj.minN(r_i),obj.maxN(r_i));
                end
                % Normalize the org range into the new range of tempo 
                obj.initial_n = (obj.initial_n - min(obj.initial_n)) ./ ...
                    (max(obj.initial_n) - min(obj.initial_n));
                vecminN = obj.minN(obj.initial_r);
                vecmaxN = obj.maxN(obj.initial_r);
                obj.initial_n = obj.initial_n .* (vecmaxN - vecminN) + vecminN; 
                clear vecminN vecmaxN
            end
            % Tempo sama informed tracking or sama informed tracking
            if sum(ismember([2,4],obj.tempoInfMode)) && ~isempty(gtTempo.bpm) && ~isempty(gtSama)
                startOffset = round(gtSama(1)/obj.frame_length);
                startOffsetTime = gtSama(1);
                yOrg = y; % Optional backup if needed
                y = y(startOffset:end,:);   % Store the signal from the first sama onwards only
                obj.initial_m = ones(size(obj.initial_m));  % Set initial_m to one, it starts from there
                chgInds = find(~ismember(obj.initial_r, find(obj.rhythm2section == 1)));
                obj.initial_r(chgInds) = randi(find(obj.rhythm2section == 1),length(chgInds),1);
                obj.initial_v = obj.rhythm2pattclass(obj.initial_r);
                clear chgInds 
            end
            % compute observation likelihoods
            obs_lik = obj.obs_model.compute_obs_lik(y);
            obs_lik(obs_lik < eps) = eps;       % Prevent underflow
            % Set up peak Inference: Assume no peak inference
            infInd = 1:size(y,1);   % Indices at which inference is done
            obsInd = zeros(1,size(y,1));   % Index mask at which obs model is to be updated
            % If peak inference mode or gated observation model is ON: Find peaks
            if (obj.peakInfMode == 1) || (obj.peakInfMode == 3)
                % We need onset peaks here
                feat1 = y(:,obj.peakParam.featID);
                feat1 = feat1 + abs(min(feat1));
                feat1 = feat1./max(feat1);
                [onsPeakInd dummy] = findPeaks(feat1, 'M', obj.peakParam.mode, ...
                    'wdTol',obj.peakParam.wdTol, 'ampTol',obj.peakParam.ampTol,...
                    'prominence',obj.peakParam.prominence);
                onsPeakInd = onsPeakInd(:)';
                if (onsPeakInd(1) ~= 1)    % Hack to include first frame in the peak list
                    onsPeakInd = [1 onsPeakInd];
                end
                if (onsPeakInd(end) ~= size(y,1)) % Hack to include last frame in the peak list
                    onsPeakInd = [onsPeakInd size(y,1)];
                end
            end
            % Select indices for peak Inf modes
            if obj.peakInfMode == 1    
                infInd = onsPeakInd;
                obsInd(infInd) = 1;
            elseif obj.peakInfMode == 2
                infInd = 1:obj.peakInfSkip:size(y,1);
                obsInd(infInd) = 1;
            elseif obj.peakInfMode == 3
                infInd = 1:size(y,1);
                obsInd(onsPeakInd) = 1;
            else
                obsInd(:) = 1;
            end
            
            % Now to call forward filtering
            obj = obj.forward_filtering(obs_lik, fname, infInd, obsInd);
            [m_path, n_path, r_path, v_path] = obj.path_via_best_weight();
            
            % strip of silence state
            if obj.use_silence_state
                idx = logical(r_path<=obj.R);
            else
                idx = true(length(r_path), 1);
            end
            % Compute the section change times, beat times, downbeats,
            % patterns
            meter = zeros(2, length(r_path));
            meter(:, idx) = obj.rhythm2meter(r_path(idx), :)';
            % compute beat times and bar positions of beats
            beats = obj.find_beat_section_times(m_path, r_path, v_path, infInd);
            if strcmp(obj.pattern_size, 'bar') || strcmp(obj.pattern_size, 'section')
                tempArr = (obj.Meff(r_path(idx)) * obj.frame_length);
                tempo = meter(1, idx)' .* 60 .* n_path(idx)' ./ ...
                    tempArr(:);
                clear tempArr
            else
                tempo = 60 .* n_path(idx) / (obj.M * obj.frame_length);
            end
            tempo = tempo(:)';
            infInd = [(1:(startOffset-1)) (infInd+startOffset-1)];
            beats(:,1) = beats(:,1) + startOffsetTime;
            if sum(ismember([2,4],obj.tempoInfMode)) && ~isempty(gtTempo.bpm) && ~isempty(gtSama)  % hack
                if beats(1,2) == 2  % If the first beat found was the second beat
                    % The instance of sama that was given must be included in the output! 
                    beats = [startOffsetTime 1 1; beats];
                    disp('Added the input sama to the output beat sequence');
                end
            end
            results{1} = beats;
            tempo = [-1*ones(1,startOffset-1) tempo];
            results{2} = tempo;
            meter = [-1*ones(2,startOffset-1) meter];
            results{3} = meter;
            r_path = [-1*ones(1,startOffset-1) r_path];
            results{4} = r_path;
            v_path = [-1*ones(1,startOffset-1,1) v_path];
            results{5} = v_path;
            ts = (infInd - 0.5) * obj.frame_length;
            results{6} = ts;           
            % Clear old variables
            clear startOffset startOffsetTime gtTempo gtSama
        end
        
        function [ joint_prob ] = compute_joint_of_sequence(obj, state_sequence, obs_lik)
            %obs_lik [R x barPos x nFrames]
            nFrames = size(obs_lik, 3);
            eval_lik = @(x, y) obj.compute_obs_lik(x, y, obs_lik, obj.M / obj.barGrid);
            joint_prob.obslik = zeros(nFrames, 1);
            n_diffs = diff(double(state_sequence(:, 2))); % nDiff(1) = s_s(2) - s_s(1)
            joint_prob.trans = log(normpdf(n_diffs/obj.M, 0, obj.sigma_N));
            for iFrame = 1:nFrames
                joint_prob.obslik(iFrame) = log(eval_lik([state_sequence(iFrame, 1), state_sequence(iFrame, 3)], iFrame));
            end
            joint_prob.sum = sum(joint_prob.trans) + sum(joint_prob.obslik);
        end
        
        function obj = convert_old_model_to_new(obj)
            % check dimensions of member variables. This function might be removed
            % in future, but is important for compatibility with older models
            % (in old models Meff and
            % rhythm2meter_state are row vectors [1 x K] but should be
            % column vectors)
            obj.Meff = obj.Meff(:);
            obj.rhythm2meter_state = obj.rhythm2meter_state(:);
            % In old models, pattern change probability was not saved as
            % matrix [RxR]
            if (length(obj.pr(:)) == 1) && (obj.R > 1)
                % expand pr to a matrix [R x R]
                % transitions to other patterns
                pr_mat = ones(obj.R, obj.R) * (obj.pr / (obj.R-1));
                % pattern self transitions
                pr_mat(logical(eye(obj.R))) = (1-obj.pr);
                obj.pr = pr_mat;
            end
            if isempty(obj.rhythm2meter)
                obj.rhythm2meter = obj.meter_state2meter(:, ...
                    obj.rhythm2meter_state)';
            end
        end
        
    end
    
    methods (Access=protected)
        
        
        function lik = compute_obs_lik(obj, states_m_r, iFrame, obslik, m_per_grid)
            % states_m_r:   is a [nParts x 2] matrix, where (:, 1) are the
            %               m-values and (:, 2) are the r-values
            % obslik:       likelihood values [R, barGrid, nFrames]
            subind = floor((states_m_r(:, 1)-1) / m_per_grid) + 1;
            obslik = obslik(:, :, iFrame);
            try
                ind = sub2ind([obj.R, obj.barGrid], states_m_r(:, 2), subind(:));
                lik = obslik(ind);
                lik = lik(:);
            catch exception
                fprintf('dimensions R=%i, barGrid=%i, states_m=%.2f - %.2f, subind = %i - %i, m_per_grid=%.2f\n', ...
                    obj.R, obj.barGrid, min(states_m_r(:, 1)), max(states_m_r(:, 1)), ...
                    min(subind), max(subind), m_per_grid);
            end
        end
        
        % Compute observation likelihood of a block
        function lik = compute_obs_lik_blk(obj, states_m, iFrame, obslik, m_per_grid)
            % states_m_r:   is a [nParts x 1] matrix, where (:, 1) are the
            %               m-values 
            % R is the total number of patterns
            % obslik:       likelihood values [R, barGrid, nFrames]
            % lik is [nParts x R] matrix of likelihood for all patterns 
            nParts = size(states_m,1);
            subind = floor((states_m-1) / m_per_grid) + 1;
            obslik = obslik(:, :, iFrame);
            %             r_ind = bsxfun(@times, (1:obj.R)', ones(1, obj.nParticles));
            try
                states_r = ones(nParts,1)*[1:obj.R];
                ind = sub2ind([obj.R, obj.barGrid], states_r(:), repmat(subind(:), obj.R,1));
                lik = reshape(obslik(ind),nParts,obj.R);
            catch exception
                fprintf('dimensions R=%i, barGrid=%i, states_m=%.2f - %.2f, subind = %i - %i, m_per_grid=%.2f\n', ...
                    obj.R, obj.barGrid, min(states_m(:, 1)), max(states_m(:, 1)), ...
                    min(subind), max(subind), m_per_grid);
            end
        end
        function [m_path, n_path, r_path, v_path] = path_via_best_weight(obj)
            % use particle with highest weight
            % ------------------------------------------------------------
            [~, bestParticle] = max(obj.particles.weight);
            m_path = obj.particles.m(bestParticle, :);
            r_path = obj.particles.r(bestParticle, :);
            n_path = obj.particles.n(bestParticle, :)';
            v_path = obj.particles.v(bestParticle, :)';
            
            m_path = m_path(:)';
            n_path = n_path(:)';
            r_path = r_path(:)';
            v_path = v_path(:)';
        end
        
        function obj = forward_filtering(obj, obs_lik, fname, infInd, obsInd)
            nTotFrames = size(obs_lik, 3);  % Number of total observation frames
            nFrames = length(infInd);      % This is the number of peaks
            if obj.save_inference_data
                logP_data_pf = log(zeros(obj.nParticles, 5, nFrames, 'single'));
            end
            % initialize particles
            iFrame = 1;
            m = zeros(obj.nParticles, nFrames, 'single');
            n = zeros(obj.nParticles, nFrames, 'single');
            r = zeros(obj.nParticles, nFrames, 'single');
            v = zeros(obj.nParticles, nFrames, 'single');
            m(:, iFrame) = obj.initial_m;
            n(:, iFrame) = obj.initial_n;
            r(:, iFrame) = obj.initial_r;
            v(:, iFrame) = obj.initial_v;
            % Initializing values
            if  obj.patt_trans_opt == 3     % Full model
                % pattWtMat = obj.pattWtMatInit;
                % But equal weight is not the best, the pattWtMat should be
                % set based on the initialization of particles also
                maskMat = kron(eye(obj.V),ones(obj.R/obj.V));
                pattParticleMask = maskMat(obj.initial_r,:);
                weightMat = log(double(pattParticleMask > 0) .* obj.pattWtMatInit);
                weightMat = normalizeLogProb(weightMat) - log(obj.nParticles);
                % Or start with flat weights
                % weightMat = ones(obj.nParticles,obj.R) .* log(1./(obj.nParticles*obj.R));
                weight = logsumexp(weightMat,2);
                currBarStart = obj.currBarStartInit;
            elseif obj.patt_trans_opt == 2  % ISMIR'15: mixture obs model
                % pattParticleMask = (obj.pr(obj.initial_r,:) > 0);
                maskMat = kron(eye(obj.V),ones(obj.R/obj.V));
                pattParticleMask = maskMat(obj.initial_r,:);
                weight = obj.prprior(obj.initial_r);
                weight = log(weight/sum(weight));
                weight = weight(:);
                % set based on the initialization of particles also
                % weight = log(ones(obj.nParticles,1)/obj.nParticles);
            else
                weight = log(obj.prprior(obj.initial_r));
                weight = weight(:);
                % Or start with flat weights
                % weight = log(ones(obj.nParticles,1)/obj.nParticles);
                weight = normalizeLogspace(weight);
            end
            % observation probability
            eval_lik = @(x, y) obj.compute_obs_lik(x, y, obs_lik, obj.M / ...
                obj.barGrid);
            eval_lik_blk = @(x, y) obj.compute_obs_lik_blk(x, y, obs_lik, obj.M /...
                obj.barGrid);
            
            if obj.patt_trans_opt == 3      % Full model
                obsMat = eval_lik_blk(m(:, iFrame), iFrame);
                obsMat = obsMat.*pattParticleMask;
                weightMat = weightMat + log(obsMat);    % Add to prior
                weightMat = normalizeLogspace2D(weightMat);   % Normalize
                weight = logsumexp(weightMat,2);
            elseif obj.patt_trans_opt == 2  % ISMIR'15: mixture obs model
                obsMat = eval_lik_blk(m(:, iFrame), iFrame);
                obsMat = obsMat.*pattParticleMask;
                weight = weight + log(sum(obsMat,2));    % Add to prior
                weight = normalizeLogspace(weight);
            else
                obs = eval_lik([obj.initial_m, obj.initial_r], iFrame);
                weight = normalizeLogspace(weight + log(obs(:)));
            end

            if obj.resampling_scheme > 1
                % divide particles into clusters by kmeans
                groups = obj.divide_into_fixed_cells([m(:, iFrame), n(:, iFrame), r(:, iFrame)], [obj.M; obj.N; obj.R], obj.n_initial_clusters);
                n_clusters = zeros(nFrames, 1);
                n_clusters(iFrame) = length(unique(groups));
            else
                groups = ones(obj.nParticles, 1);
            end
            if obj.save_inference_data
                % save particle data for visualizing
                logP_data_pf(:, 1, iFrame) = m(:, iFrame);
                logP_data_pf(:, 2, iFrame) = n(:, iFrame);
                logP_data_pf(:, 3, iFrame) = r(:, iFrame);
                logP_data_pf(:, 4, iFrame) = weight;
                logP_data_pf(:, 5, iFrame) = groups;
            end
            resampling_frames = zeros(nFrames, 1);
            perc = round(0.1*nFrames);
            fprintf('Particle filtering started');
            for iFrame=2:nFrames
                % transition from iFrame-1 to iFrame
                m(:, iFrame) = m(:, iFrame-1) + n(:, iFrame-1)*(infInd(iFrame)-infInd(iFrame-1));
                m(:, iFrame) = mod(m(:, iFrame) - 1, ...
                    obj.Meff(r(:, iFrame-1))) + 1;
                % Change the ones for which the bar changed
                newBars = find(m(:, iFrame) < m(:, iFrame-1));
                % Pattern transitions to be handled here
                r(:, iFrame) = r(:, iFrame-1);
                v(:, iFrame) = v(:, iFrame-1);
                
                if obj.patt_trans_opt == 3          
                    % Full model
                    [~, maxInds] = max(weightMat(newBars,:),[],2);
                    for mInd = 1:length(newBars)
                        r(newBars(mInd),currBarStart(newBars(mInd)):iFrame-1) = maxInds(mInd);
                        v(newBars(mInd),currBarStart(newBars(mInd)):iFrame-1) = ...
                            obj.rhythm2pattclass(r(newBars(mInd),currBarStart(newBars(mInd)):iFrame-1));
                        currBarStart(newBars(mInd)) = iFrame;
                        r(newBars(mInd),iFrame) = randsample(obj.R, 1, true, obj.pr(r(newBars(mInd),iFrame-1),:));  % This is only temporary, to compute the m transitions
                        v(newBars(mInd),iFrame) = obj.rhythm2pattclass(r(newBars(mInd),iFrame)); % This as well is temporary, to compute the m transitions
                        pattParticleMask(newBars(mInd),:) = (obj.pr(r(newBars(mInd),iFrame-1),:) > 0);  % Careful here!
                    end
                    % currBarStart(newBars) = iFrame;
                    weightMat(newBars,:) = bsxfun(@plus,log(obj.pr(r(newBars,iFrame-1),:)), weight(newBars));
                    % Observation likelihoods
                    if obsInd(infInd(iFrame))
                        obsMat = eval_lik_blk(m(:, iFrame), infInd(iFrame));   %% CHANGE AS NEEDED
                        obsMat = obsMat.*pattParticleMask;
                        weightMat = weightMat + log(obsMat);    % Add to weights
                        weightMat = normalizeLogspace2D(weightMat);   % Normalize
                        weight = logsumexp(weightMat,2);
                    end
                elseif obj.patt_trans_opt == 2      % Obs_model, not possible with section length patterns
                    % Update the pattern weight matrix
                    if obsInd(infInd(iFrame))
                        obsMat = eval_lik_blk(m(:, iFrame), infInd(iFrame));        
                        obsMat = obsMat.*pattParticleMask;
                        weight = weight + log(sum(obsMat,2));
                        weight = normalizeLogspace(weight);
                    end
                elseif obj.patt_trans_opt == 1 || obj.patt_trans_opt == 0
                    % Sampling from prior
                    for rInd = 1:length(newBars)
                        r(newBars(rInd), iFrame) = ...
                            randsample(obj.R, 1, true, obj.pr(r(newBars(rInd),iFrame-1),:));
                        v(newBars(rInd), iFrame) = obj.rhythm2pattclass(r(newBars(rInd), iFrame));
                    end
                    if obsInd(infInd(iFrame))
                        obs = eval_lik([m(:, iFrame), r(:, iFrame)], infInd(iFrame));
                        weight = weight + log(obs(:));
                        weight = normalizeLogspace(weight);
                    end
                % else
                    % obs = eval_lik([m(:, iFrame), r(:, iFrame)], infInd(iFrame));
                    % weight = weight + log(obs(:));
                    % weight = normalizeLogspace(weight);
                end
                % Resampling
                % ------------------------------------------------------------
                if obj.resampling_interval == 0
                    Neff = 1/sum(exp(weight).^2);
                    do_resampling = (Neff < obj.ratio_Neff * obj.nParticles);
                else
                    do_resampling = (rem(infInd(iFrame), obj.resampling_interval) == 0);
                end
                if do_resampling && (iFrame < nFrames)
                    if obj.patt_trans_opt == 3
                        weightOld = weight;
                    end
                    resampling_frames(iFrame) = iFrame;
                    if obj.resampling_scheme == 2 || obj.resampling_scheme == 3 % MPF or AMPF
                        groups = obj.divide_into_clusters([m(:, iFrame), n(:, iFrame-1), r(:, iFrame)], [obj.M; obj.N; obj.R], groups);
                        n_clusters(iFrame) = length(unique(groups));
                    end
                    [weight, groups, newIdx] = obj.resample(weight, groups);
                    weight = weight(:);
                    m(:, 1:iFrame) = m(newIdx, 1:iFrame);
                    r(:, 1:iFrame) = r(newIdx, 1:iFrame);
                    n(:, 1:iFrame) = n(newIdx, 1:iFrame);
                    v(:, 1:iFrame) = v(newIdx, 1:iFrame);
                    if obj.patt_trans_opt == 3
                        currBarStart = currBarStart(newIdx);
                        pattParticleMask = pattParticleMask(newIdx,:);
                        weightMat = weightMat(newIdx,:);
                        weightOld = weightOld(newIdx);
                        weightMat = bsxfun(@plus, weightMat, weight - weightOld); 
                    elseif obj.patt_trans_opt == 2
                        pattParticleMask = pattParticleMask(newIdx,:);
                    end
                end
                % Tempo transition from iFrame-1 to iFrame
                % n(:, iFrame) = n(:, iFrame-1) + randn(obj.nParticles, 1) * obj.sigma_N * obj.M;
                % n(:, iFrame) = n(:, iFrame-1) + obj.sigma_N * randn(obj.nParticles, 1) ...
                %    .* obj.Meff(r(:, iFrame-1))/obj.M .* ...
                %    (obj.maxN(r(:, iFrame-1)) - obj.minN(r(:, iFrame-1)));
                n(:, iFrame) = n(:, iFrame-1) .* (1 + randn(obj.nParticles, 1) .* obj.sigma_N .* obj.Meff(r(:, iFrame-1))/obj.M);
                out_of_range = n(:, iFrame) > obj.maxN(r(:, iFrame));
                n(out_of_range, iFrame) = obj.maxN(r(out_of_range, iFrame));
                out_of_range = n(:, iFrame) < obj.minN(r(:, iFrame));
                n(out_of_range, iFrame) = obj.minN(r(out_of_range, iFrame));
                
                if obj.save_inference_data
                    % save particle data for visualizing
                    logP_data_pf(:, 1, iFrame) = m(:, iFrame);
                    logP_data_pf(:, 2, iFrame) = n(:, iFrame);
                    logP_data_pf(:, 3, iFrame) = r(:, iFrame);
                    logP_data_pf(:, 4, iFrame) = v(:, iFrame);
                    logP_data_pf(:, 5, iFrame) = weight;
                    logP_data_pf(:, 6, iFrame) = groups;
                end
                if rem(iFrame, perc) == 0
                    fprintf('.');
                end
            end
            fprintf('\n');
            
            obj.particles.m = m;
            obj.particles.n = n;
            obj.particles.r = r;
            obj.particles.v = v;
            obj.particles.ts = infInd;   % Extra, not needed
            obj.particles.weight = weight;
            
            fprintf('    Average resampling interval: %.2f frames\n', ...
                mean(diff(resampling_frames(resampling_frames>0))));
            if obj.resampling_scheme > 1
                fprintf('    Average number of clusters: %.2f clusters\n', mean(n_clusters(n_clusters>0)));
            end
            if obj.save_inference_data
                save(['./', fname, '_pf.mat'], 'logP_data_pf');
            end
        end
        
        function beats = find_beat_times(obj, positionPath, rhythmPath, infInd)
            %   Find beat times from sequence of bar positions of the PF beat tracker
            % ----------------------------------------------------------------------
            %INPUT parameter:
            % positionPath             : sequence of position states
            % rhythmPath                : sequence of meter states
            %                           NOTE: so far only median of sequence is used !
            % nBarPos                  : bar length in bar positions (=M)
            % framelength              : duration of being in one state in [sec]
            %
            %OUTPUT parameter:
            %
            % beats                    : [nBeats x 2] beat times in [sec] and
            %                           beatnumber
            %
            % 29.7.2012 by Florian Krebs
            % ----------------------------------------------------------------------
            numframes = length(positionPath);
            beatpositions = cell(obj.R, 1);
            for i_r=1:obj.R
                is_compund_meter = ismember(obj.rhythm2meter(i_r, 1), ...
                    [6, 9, 12]);
                if is_compund_meter
                    beatpositions{i_r} = round(linspace(1, obj.Meff(i_r), ...
                        obj.rhythm2meter(i_r, 1) / 3 + 1));
                else % simple meter
                    beatpositions{i_r} = round(linspace(1, obj.Meff(i_r), ...
                        obj.rhythm2meter(i_r, 1) + 1));
                end
                beatpositions{i_r} = beatpositions{i_r}(1:end-1);
            end
            beatno = [];
            beats = [];
            for i = 1:numframes-1
                if rhythmPath(i) == 0
                    continue;
                end
                for beat_pos = 1:length(beatpositions{rhythmPath(i)})
                    if positionPath(i) == beatpositions{rhythmPath(i)}(beat_pos)
                        % current frame = beat frame
                        beats = [beats; [infInd(i), beat_pos]];
                        break;
                    elseif ((positionPath(i+1) > beatpositions{rhythmPath(i)}(beat_pos)) && (positionPath(i+1) < positionPath(i)))
                        % bar transition between frame i and frame i+1
                        bt = interp1([positionPath(i); obj.M+positionPath(i+1)],[infInd(i); infInd(i+1)],obj.M+beatpositions{rhythmPath(i)}(beat_pos));
                        beats = [beats; [round(bt), beat_pos]];
                        break;
                    elseif ((positionPath(i) < beatpositions{rhythmPath(i)}(beat_pos)) && (beatpositions{rhythmPath(i)}(beat_pos) < positionPath(i+1)))
                        % beat position lies between frame i and frame i+1
                        bt = interp1([positionPath(i); positionPath(i+1)],[infInd(i); infInd(i+1)],beatpositions{rhythmPath(i)}(beat_pos));
                        beats = [beats; [round(bt), beat_pos]];
                        break;
                    end
                end
            end
            beats(:, 1) = beats(:, 1) * obj.frame_length;
        end
        
        
        function beats = find_beat_section_times(obj, positionPath, rhythmPath, sectionPath, infInd)
            %   Find beat times from sequence of bar positions of the PF beat tracker
            % ----------------------------------------------------------------------
            %INPUT parameter:
            % positionPath             : sequence of position states
            % rhythmPath                : sequence of rhythm states
            % sectionPath                : sequence of section states
            %                           NOTE: so far only median of sequence is used !
            % nBarPos                  : bar length in bar positions (=M)
            % framelength              : duration of being in one state in [sec]
            %
            %OUTPUT parameter:
            %
            % beats                    : [nBeats x 3] beat times in [sec], 
            %                           beatnumber, section number
            %
            % 29.7.2012 by Florian Krebs
            % ----------------------------------------------------------------------
            numframes = length(positionPath);
            
            beatpositions = cell(obj.R, 1);
            beatnum = cell(obj.R, 1);
            secnum = cell(obj.R, 1);
            
            for i_r=1:obj.R
                % is_compund_meter = ismember(obj.rhythm2meter(i_r, 1), ...
                %    [6, 9, 12]);
                %if is_compund_meter
                %    beatpositions{i_r} = round(linspace(1, obj.Meff(i_r), ...
                %        obj.rhythm2meter(i_r, 1) / 3 + 1));
                %else % simple meter
                %    beatpositions{i_r} = round(linspace(1, obj.Meff(i_r), ...
                %        obj.rhythm2meter(i_r, 1) + 1));
                % end
                beatpositions{i_r} = round(linspace(1, obj.Meff(i_r), ...
                        obj.rhythm2len_num(i_r) + 1));
                beatpositions{i_r} = beatpositions{i_r}(1:end-1);
                beatnum{i_r} = (1:length(beatpositions{i_r})) + ...
                    obj.pattInfo.class(obj.rhythm2pattclass(i_r),4) - 1;    % This has the beat number
                secnum{i_r} = ones(1,length(beatpositions{i_r})) * obj.rhythm2section(i_r); 
            end
            beats = [];
            for i = 1:numframes-1
                if rhythmPath(i) == 0
                    continue;
                end
                for beat_pos = 1:length(beatpositions{rhythmPath(i)})
                    %if positionPath(i) == beatpositions{rhythmPath(i)}(beat_pos)
                    %    % current frame = beat frame
                    %    beats = [beats; [infInd(i), beatnum{rhythmPath(i)}(beat_pos) secnum{rhythmPath(i)}(beat_pos)]];
                    %    break;
                    if ((positionPath(i+1) >= beatpositions{rhythmPath(i)}(beat_pos)) && (positionPath(i+1) < positionPath(i)))
                        % bar transition between frame i and frame i+1
                        bt = interp1([positionPath(i); obj.Meff(rhythmPath(i))+positionPath(i+1)],...
                            [infInd(i); infInd(i+1)],obj.Meff(rhythmPath(i))+beatpositions{rhythmPath(i)}(beat_pos));
                        beats = [beats; [round(bt), beatnum{rhythmPath(i+1)}(beat_pos) secnum{rhythmPath(i+1)}(beat_pos)]];
                        break;
                    elseif ((positionPath(i) < beatpositions{rhythmPath(i)}(beat_pos)) && (beatpositions{rhythmPath(i)}(beat_pos) < positionPath(i+1)))
                        % beat position lies between frame i and frame i+1
                        bt = interp1([positionPath(i); positionPath(i+1)],[infInd(i); infInd(i+1)],beatpositions{rhythmPath(i)}(beat_pos));
                        beats = [beats; [round(bt), beatnum{rhythmPath(i)}(beat_pos) secnum{rhythmPath(i)}(beat_pos)]];
                        break;
                    end
                end
            end
            beats(:, 1) = beats(:, 1) * obj.frame_length;
        end
        
        function [groups] = divide_into_clusters(obj, states, state_dims, groups_old)
            % states: [nParticles x nStates]
            % state_dim: [nStates x 1]
            % groups_old: [nParticles x 1] group labels of the particles
            %               after last resampling step (used for initialisation)
            warning('off');
            [group_ids, ~, IC] = unique(groups_old);
            %             fprintf('    %i groups >', length(group_ids));
            k = length(group_ids); % number of clusters
            
            % adjust the range of each state variable to make equally
            % important for the clustering
            points = zeros(obj.nParticles, length(state_dims)+1);
            points(:, 1) = (sin(states(:, 1) * 2 * pi ./ ...
                obj.Meff(states(:, 3))) + 1) * ...
                obj.state_distance_coefficients(1);
            points(:, 2) = (cos(states(:, 1) * 2 * pi ./ ...
                obj.Meff(states(:, 3))) + 1) * ...
                obj.state_distance_coefficients(1);
            points(:, 3) = states(:, 2) * ...
                obj.state_distance_coefficients(2);
            if obj.patt_trans_opt <= 2    
                points(:, 4) =(states(:, 3)-1) * ...
                    obj.state_distance_coefficients(3) + 1;
            else
                % No effect of patterns in mix and acc models
                points(:, 4) = (states(:, 3)-1) * 0 + 1;     
            end
            
            % compute centroid of clusters
            % TODO: vectorise!
            centroids = zeros(k, length(state_dims)+1);
            for i_dim=1:size(points, 2)
                %                 centroids(iCluster, :) = mean(points(groups_old == group_ids(iCluster) , :));
                centroids(:, i_dim) = accumarray(IC, points(:, i_dim), [], @mean);
            end
            % do k-means clustering
            options = statset('MaxIter', 1);
            [groups, centroids, total_dist_per_cluster] = kmeans(points, k, 'replicates', 1, ...
                'start', centroids, 'emptyaction', 'drop', 'Distance', 'sqEuclidean', 'options', options);
            %             remove empty clusters
            total_dist_per_cluster = total_dist_per_cluster(~isnan(centroids(:, 1)));
            centroids = centroids(~isnan(centroids(:, 1)), :);
            [group_ids, ~, j] = unique(groups);
            group_ids = 1:length(group_ids);
            groups = group_ids(j);
            groups = groups(:); % Line added to ensure groups is a column vector
            % Check if centroids are too close
            if size(centroids,1) == 1   % If there is only a single centroid, dont try to merge
                merging = 0;
            else % Try merging
                merging = 1;
            end
            merged = 0;
            while merging
                D = squareform(pdist(centroids, 'euclidean'), 'tomatrix');
                ind = (tril(D, 0) > 0);
                D(ind) = nan;
                D(logical(eye(size(centroids, 1)))) = nan;
                
                % find minimum distance
                [min_D, arg_min] = min(D(:));
                if min_D > obj.cluster_merging_thr,
                    merging = 0;
                else
                    [c1, c2] = ind2sub(size(D), arg_min);
                    groups(groups==c2(1)) = c1(1);
                    % new centroid
                    centroids(c1(1), :) = mean(points(groups==c1(1), :));
                    % squared Euclidean distance
                    total_dist_per_cluster(c1(1)) = sum(sum(bsxfun(@minus, points(groups==c1(1), :), centroids(c1(1), :)).^2));
                    if length(c1) == 1,  merging = 0;  end
                    % remove old centroid
                    centroids = centroids([1:c2(1)-1, c2(1)+1:end], :);
                    total_dist_per_cluster = total_dist_per_cluster([1:c2(1)-1, c2(1)+1:end]);
                    merged = 1;
                end
            end
            % check if cluster spread is too high
            split = 0;
            [group_ids, ~, j] = unique(groups);
            group_ids = 1:length(group_ids);
            groups = group_ids(j);
            groups = groups(:);   % Line added to ensure groups is a column vector
            n_parts_per_cluster = hist(groups, 1:length(group_ids));
            separate_cl_idx = find((total_dist_per_cluster ./ n_parts_per_cluster') > obj.cluster_splitting_thr);
            %             separate_cl_idx = find(mean_dist_per_cluster > obj.cluster_splitting_thr);
            for iCluster = 1:length(separate_cl_idx)
                % find particles that belong to the cluster to split
                parts_idx = find((groups == separate_cl_idx(iCluster)));
                % put second half into a new group
                groups(parts_idx(round(length(parts_idx)/2)+1:end)) = max(groups) + 1;
                % update centroid
                centroids(separate_cl_idx(iCluster), :) = mean(points(parts_idx(1:round(length(parts_idx)/2)), :), 1);
                % add new centroid
                centroids = [centroids; mean(points(parts_idx(round(length(parts_idx)/2)+1:end), :), 1)];
                split = 1;
            end
            if split || merged
                try
                    [groups, ~, ~] = kmeans(points, [], 'replicates', 1, 'start', centroids, 'emptyaction', 'drop', ...
                        'Distance', 'sqEuclidean', 'options', options);
                catch exception
                    centroids
                    error('Problem\n');
                end
                [group_ids, ~, j] = unique(groups);
                group_ids = 1:length(group_ids);
                groups = group_ids(j); 
                groups = groups(:); % Line added to ensure groups is a column vector
            end
            warning('on');
        end
        
        
        function [weight, groups, newIdx] = resample(obj, weight, groups)
            if obj.resampling_scheme == 0
                newIdx = obj.resampleSystematic(exp(weight));
                weight = log(ones(obj.nParticles, 1) / obj.nParticles);
            elseif obj.resampling_scheme == 1 % APF
                % warping:
                w = exp(weight);
                f = str2func(obj.warp_fun);
                w_warped = f(w);
                newIdx = obj.resampleSystematic(w_warped);
                w_fac = w ./ w_warped;
                weight = log( w_fac(newIdx) / sum(w_fac(newIdx)) );
            elseif obj.resampling_scheme == 2 % K-MEANS
                % k-means clustering
                [newIdx, weight, groups] = obj.resample_in_groups(groups, ...
                    weight, obj.n_max_clusters);
                weight = weight';
            elseif obj.resampling_scheme == 3 % APF + K-MEANS
                % apf and k-means
                [newIdx, weight, groups] = obj.resample_in_groups(groups, ...
                    weight, obj.n_max_clusters, str2func(obj.warp_fun));
                weight = weight';
            else
                fprintf('WARNING: Unknown resampling scheme!\n');
            end
        end
    end
    
    methods (Static)
        
        outIndex = resampleSystematic( w, n_samples );
        
        outIndex = resampleSystematicInGroups( w, n_samples );
        
        function [groups] = divide_into_fixed_cells(states, state_dims, nCells)
            % divide space into fixed cells with equal number of grid
            % points for position and tempo states
            % states: [nParticles x nStates]
            % state_dims: [nStates x 1]
            groups = zeros(size(states, 1), 1);
            n_r_bins = state_dims(3);
            n_n_bins = floor(sqrt(nCells/n_r_bins));
            n_m_bins = floor(nCells / (n_n_bins * n_r_bins));
            
            m_edges = linspace(1, state_dims(1) + 1, n_m_bins + 1);
            n_edges = linspace(min(states(:,2))*0.9, state_dims(2)*1.1, n_n_bins + 1);
            for iR=1:state_dims(3)
                % get all particles that belong to pattern iR
                ind = find(states(:, 3) == iR);
                [~, BIN_m] = histc(states(ind, 1), m_edges);
                [~, BIN_n] = histc(states(ind, 2), n_edges);
                for m = 1:n_m_bins
                    for n=1:n_n_bins
                        ind2 = intersect(ind(BIN_m==m), ind(BIN_n==n));
                        groups(ind2) = sub2ind([n_m_bins, n_n_bins, state_dims(3)], m, n, iR);
                    end
                end
            end
            if sum(groups==0) > 0
                error('group assignment failed\n')
            end
        end
        
        [outIndex, outWeights, groups] = resample_in_groups(groups, weights, n_max_clusters, warp_fun);
        
    end
    
    
    
end
