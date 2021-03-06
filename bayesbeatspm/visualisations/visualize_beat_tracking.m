function [] = visualize_beat_tracking( data_fln, blocksize, audio_fln)
        % VISUALIZE_BEAT_TRACKING visualize (filtering) posterior probability of the hidden
        % states given the observations
        %
        % INPUT:
        %   data_fln:  filename of fixed grid (FG) and/or particle grid (PG) data
        %           FG has the ending '_hmm.mat' or '_func.mat'
        %               '_hmm.mat' has the fields
        %                   'logP_data', 'M', 'N', 'R', 'frame_length', 'obs_lik'
        %                   logP_data is a sparse matrix of size [n_states x n_frames]
        %               '_func.mat' contains the variables 'data' and 'evaluate_p'
        %           PG has the ending '_pf.mat' and contains 'logP_data_pf' which
        %               is of size [n_particles, 5, n_frames]:
        %               where dim 2 is : m x n x r x w x group
        %   blocksize: blocksize [in frames] of the video, where one frame 
        %   corresponds to the frame_length given in the data structure
        %
        % OUTPUT
        %   video
        % -----------------------------------------------------------------------
        close all
        % -----------------------------------------------------------------------
        % parse input arguments
        % -----------------------------------------------------------------------
        synth_data_ok = 0; hmm_data_ok = 0; pf_data_ok = 0; anns_data_ok = 0;
        [path, fname, ~] = fileparts(data_fln);
        if isempty(path), path = './'; end
        fln = fullfile(path, [fname, '_synth.mat']);
        if exist(fln, 'file')
            load(fln);
            synth_data_ok = 1;
            R = 1; % only 1D data
            frame_length = 0.02;
            n_m_points = 70;
            n_n_points = 20;
            M = data.M; N = data.N;
            m_points = round(linspace(1,M,n_m_points));
            n_points = round(linspace(1, N, n_n_points));
            m_grid = repmat(m_points, length(n_points), 1)';
            n_grid = repmat(n_points, length(m_points), 1);
            grid = [m_grid(:), n_grid(:)];
            nFrames = size(data.weights, 2);
        end
        % load hmm data
        fln = fullfile(path, [fname, '_hmm.mat']);
        if exist(fln, 'file')
            load(fln);
            hmm_data_ok = 1;
            nFrames = size(logP_data, 2);
            [~, nPos, ~] = size(obs_lik);
            max_lik = max(obs_lik(:));
        else
            fprintf('%s not found\n', fln);
            return
        end
        % load pf data
        fln = fullfile(path, [fname, '_pf.mat']);
        if exist(fln, 'file')
            load(fln);
            pf_data_ok = 1;
        end
        % load annotations
        fln = fullfile(path, [fname, '_anns.mat']);
        if exist(fln, 'file')
            load(fln);
            anns_data_ok = 1;
        end
        % check if everything there
        if (synth_data_ok + hmm_data_ok + pf_data_ok) == 0
            error('No data found\n')
        end
        % -----------------------------------------------------------------------
        % initialize graphical display
        % -----------------------------------------------------------------------
        hf = figure;
        % get rid of an annoying warning: "RGB color data not yet supported in Painter's mode."
        set(gcf, 'renderer','zbuffer');
        col_weights = autumn;
        col_part_sets = [0,  0,  1;
            0,  1,  0;
            0,  1,  1;
            1,  0,  0;
            1,  0,  1;
            1,  1,  0];
        col_part_sets = [col_part_sets; col_part_sets*0.5; col_part_sets*0.25; col_part_sets*0.75];
        col_part_sets = [col_part_sets; col_part_sets*0.3; col_part_sets*0.6];
        colormap(gray);
        set(gca,'YDir','normal')
        nPlots = R * 2;
        h_sp = zeros(nPlots, 1);
        if exist('movie_fln', 'var')
            movie_fln = fullfile('./', movie_fln);
        else
            movie_fln = fullfile('./', [fname, '.avi']);
        end
        aviobj = avifile(movie_fln, 'fps', 1 / ( blocksize * frame_length));
        % aviobj = VideoWriter(movie_fln, 'Motion JPEG AVI');
        % aviobj.FrameRate = 1/(2*frame_length);
        % open(aviobj);
        % % ax=get(hf,'Position');
        
        % size/resolution of the figure
        set(hf, 'Position', [51 1 500 400]);
        
        % position for each rhythmic pattern:
        y_pos = [0.61; 0.17];
        
        % -----------------------------------------------------------------------
        % visualize
        % -----------------------------------------------------------------------
        % security check: use max of 1000 frames
        if nFrames/blocksize > 1000
            error(['Please reduce length of the data to be visualised or', ...
                ' increase blocksize\n']);
        end
        for iFrame = 1:blocksize:nFrames
            for iR=1:R
                if hmm_data_ok
                    plot_id = (iR-1)*R+1;
                    h_sp(plot_id) = subplot(nPlots, 1, plot_id);
                    start_ind = sub2ind([round(M/x_fac), N, R], 1, 1, iR);
                    end_ind = sub2ind([round(M/x_fac), N, R], round(M/x_fac), N, iR);
                    frame = reshape(logP_data(start_ind:end_ind, iFrame), round(M/x_fac), N);
                    if nnz(frame) > 0
                        max_h = max(frame(:));
                        min_h = max([min(frame(~isinf(frame(:)))), max_h - 50]);
                    else
                        frame(:) = min_h;
                    end
                    imagesc(frame');
                    if isinf(min_h)
                        caxis([-99 -98])
                    else
                        caxis([min_h max_h])
                    end
                end
                if synth_data_ok
                    p = evaluate_p(grid, iFrame);
                    p = p/sum(p);
                    p_map = reshape(p, n_m_points, n_n_points)';
                    imagesc(m_points, n_points, p_map);
                    caxis([0 0.025])
                end
                colorbar;
                hold on;
                if pf_data_ok
                    ind = (logP_data_pf(:, 3, iFrame) == iR);
                    if sum(ind) > 0
                        particle_sets = unique(logP_data_pf(:, 5, iFrame))';
                        for iCluster=particle_sets
                            is_in_cluster = (logP_data_pf(:, 5, iFrame) == iCluster);
                            valid_idx = ind & is_in_cluster;
                            if sum(valid_idx) > 0
                                col_bins = linspace(min(logP_data_pf(valid_idx, 4, iFrame)), max(logP_data_pf(valid_idx, 4, iFrame)), 64);
                                [~, bins] = histc(logP_data_pf(valid_idx, 4, iFrame), col_bins);
                                bins(bins < 1) = 1; bins(bins > 64) = 64;
                                if length(particle_sets) == 1
                                    scatter(logP_data_pf(valid_idx, 1, iFrame)./x_fac, logP_data_pf(valid_idx, 2, iFrame), ...
                                        10, col_weights(bins, :), 'filled');
                                else
                                    scatter(logP_data_pf(valid_idx, 1, iFrame)./x_fac, logP_data_pf(valid_idx, 2, iFrame), ...
                                        10, col_part_sets(iCluster, :), 'filled');
                                end
                            end
                        end
                    end
                end
                if anns_data_ok
                    if anns(1,3) == iR
                        scatter(anns(iFrame, 1)./x_fac, anns(iFrame, 2), 30, 'c', 'filled', 'MarkerEdgeColor','k');
                    end
                end
                if iR == 1, title(sprintf('Frame: %i', iFrame)); end
                set(gca,'YDir','normal')
                if ~synth_data_ok
                    h_sp(plot_id+1) = subplot(nPlots, 1, plot_id+1);
                    obs_normed = obs_lik(iR, :, iFrame)/sum(sum(obs_lik(:, :, iFrame)));
                    stairs(obs_normed);
                    xlim([1 nPos+1])
                    % Position and length of the plots !
                    set(h_sp(plot_id), 'Position', [0.1 y_pos(iR)-0.4 0.75 0.7]); % [xmin ymin xlenght ylength]);
                    set(h_sp(plot_id+1), 'Position', [0.1 y_pos(iR)-0.4-0.09 0.75 0.08]); % [xmin ymin xlenght ylength]);
                end
            end
            F = getframe(hf);
            aviobj = addframe(aviobj,F);
        end
        close(hf);
        close(aviobj);
        
        % prepare audio file and add it to the video
        if exist(audio_fln, 'file')
            [y, fs] = wavread(audio_fln);
            y = y((1:blocksize:nFrames)*frame_length*fs);
            wavwrite(y, fs, 'test.wav');
            system(['lame test.wav ', fullfile(path, [fname, '.mp3'])]);
            fprintf('Saved cut mp3 to %s\n', fullfile(path, [fname, '.mp3']));
            [path, fname, ~] = fileparts(movie_fln);
            cmd = ['ffmpeg -i ', movie_fln, ' -b 9600 -qscale 5 -acodec copy -i ', fullfile(path, [fname, '.mp3']), ...
                    ' ', fullfile(path, [fname, '_av.avi'])];
            fprintf('Executing %s\n', cmd);
            system(cmd);
        end
    end
