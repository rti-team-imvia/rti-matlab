function processRTISelective(rti_dir, fit_indices, relight_indices, params)
    % Import required packages
    import pkg_DMD.* 
    import pkg_PTM.* 
    import pkg_HSH.* 
    import pkg_fcns.*

    % Load LP file and initialize
    [Images, is_rgb] = initialize_and_load(rti_dir);
    
    % Process images based on type
    if is_rgb
        results = process_rgb_images(Images, fit_indices, relight_indices, params);
    else
        results = process_grayscale_images(Images, fit_indices, relight_indices, params);
    end
    
    % Save results and compute errors
    save_results_and_compute_errors(results, Images, is_rgb, rti_dir, relight_indices, params);
end

function [Images, is_rgb] = initialize_and_load(rti_dir)
    % First find and load the LP file
    lp_files = dir(fullfile(rti_dir, '*.lp'));
    if isempty(lp_files)
        error('No .lp file found in directory: %s', rti_dir);
    end
    fprintf('Found LP file: %s\n', lp_files(1).name);
    
    % Read LP file
    fid = fopen(fullfile(rti_dir, lp_files(1).name));
    temp_L = textscan(fid, '%s %f %f %f');
    fclose(fid);
    
    % Initialize Images structure
    Images = struct();
    Images.image_paths = temp_L{1}(2:end);  % First entry is header
    Images.LP1 = [temp_L{2}(2:end) temp_L{3}(2:end) temp_L{4}(2:end)];
    Images.nfiles = length(Images.image_paths);
    
    % Read first image to get dimensions and check if RGB
    first_img_path = fullfile(rti_dir, Images.image_paths{1});
    if ~exist(first_img_path, 'file')
        error('Cannot find first image: %s', first_img_path);
    end
    
    first_img = imread(first_img_path);
    Images.height = size(first_img,1);
    Images.width = size(first_img,2);
    is_rgb = size(first_img,3) == 3;
    
    % Load images based on type
    if is_rgb
        [Images.DataV_R, Images.DataV_G, Images.DataV_B] = load_rgb_images(rti_dir, Images);
    else
        Images.DataV = load_grayscale_images(rti_dir, Images);
    end
    
    fprintf('Loaded %d %s images of size %dx%d\n', Images.nfiles, ...
        conditional(is_rgb, 'RGB', 'grayscale'), Images.height, Images.width);
end

function [DataV_R, DataV_G, DataV_B] = load_rgb_images(rti_dir, Images)
    DataV_R = zeros(Images.nfiles, Images.height * Images.width, 'single');
    DataV_G = zeros(Images.nfiles, Images.height * Images.width, 'single');
    DataV_B = zeros(Images.nfiles, Images.height * Images.width, 'single');
    
    fprintf('Loading RGB images...\n');
    for ii = 1:Images.nfiles
        img_path = fullfile(rti_dir, Images.image_paths{ii});
        if ~exist(img_path, 'file')
            error('Cannot find image: %s', img_path);
        end
        current_img = imread(img_path);
        
        % Reshape each channel to row vector and normalize
        DataV_R(ii,:) = single(reshape(current_img(:,:,1), 1, [])) / 255;
        DataV_G(ii,:) = single(reshape(current_img(:,:,2), 1, [])) / 255;
        DataV_B(ii,:) = single(reshape(current_img(:,:,3), 1, [])) / 255;
        
        if mod(ii,10) == 0
            fprintf('Processed %d/%d images\n', ii, Images.nfiles);
        end
    end
end

function DataV = load_grayscale_images(rti_dir, Images)
    DataV = zeros(Images.nfiles, Images.height * Images.width, 'single');
    
    fprintf('Loading grayscale images...\n');
    for ii = 1:Images.nfiles
        img_path = fullfile(rti_dir, Images.image_paths{ii});
        if ~exist(img_path, 'file')
            error('Cannot find image: %s', img_path);
        end
        current_img = imread(img_path);
        if size(current_img,3) > 1
            current_img = rgb2gray(current_img);
        end
        DataV(ii,:) = single(reshape(current_img, 1, [])) / 255;
        
        if mod(ii,10) == 0
            fprintf('Processed %d/%d images\n', ii, Images.nfiles);
        end
    end
end

function results = process_rgb_images(Images, fit_indices, relight_indices, params)
    fprintf('Processing RGB channels...\n');
    
    % Process each channel
    fprintf('Processing Red channel...\n');
    results_R = fit_and_relight_channel(Images.DataV_R, Images.LP1, fit_indices, relight_indices, params, Images);
    
    fprintf('Processing Green channel...\n');
    results_G = fit_and_relight_channel(Images.DataV_G, Images.LP1, fit_indices, relight_indices, params, Images);
    
    fprintf('Processing Blue channel...\n');
    results_B = fit_and_relight_channel(Images.DataV_B, Images.LP1, fit_indices, relight_indices, params, Images);
    
    % Combine results
    results = combine_rgb_results(results_R, results_G, results_B, Images, length(relight_indices));
end

function results = process_grayscale_images(Images, fit_indices, relight_indices, params)
    fprintf('Processing grayscale images...\n');
    results = fit_and_relight_channel(Images.DataV, Images.LP1, fit_indices, relight_indices, params, Images);
end

function results = combine_rgb_results(results_R, results_G, results_B, Images, n_relight)
    results = struct();
    results.DMD = zeros(Images.height, Images.width, 3, n_relight);
    results.PTM = zeros(Images.height, Images.width, 3, n_relight);
    results.HSH = zeros(Images.height, Images.width, 3, n_relight);
    
    for i = 1:n_relight
        results.DMD(:,:,:,i) = cat(3, results_R.DMD(:,:,i), results_G.DMD(:,:,i), results_B.DMD(:,:,i));
        results.PTM(:,:,:,i) = cat(3, results_R.PTM(:,:,i), results_G.PTM(:,:,i), results_B.PTM(:,:,i));
        results.HSH(:,:,:,i) = cat(3, results_R.HSH(:,:,i), results_G.HSH(:,:,i), results_B.HSH(:,:,i));
    end
end

function channel_results = fit_and_relight_channel(channel_data, LP1, fit_indices, relight_indices, params, Images)
    % Extract fitting subset
    fit_data = channel_data(fit_indices, :);
    fit_LP = LP1(fit_indices, :);
    
    % 1. DMD Fitting
    fprintf('Fitting DMD model...\n');
    data.nb_modes = params.dmd_modes;
    [data.modal_basis, data.normal_elt, data.Up] = DMD_build_RMB(data.nb_modes, 'fixed');
    interp_Qi_all = DMD_getEigenModes(data.nb_modes, fit_LP, data.modal_basis, data.normal_elt, data.Up);
    data.coef_DMD = LeastSquares(fit_data, interp_Qi_all);
    
    % 2. PTM Fitting
    fprintf('Fitting PTM model...\n');
    dirs_fit = LP_xyz2phitheta(fit_LP);
    interp_Pi = PTM_terms(dirs_fit);
    coef_PTM = LeastSquares(fit_data, interp_Pi);
    
    % 3. HSH Fitting
    fprintf('Fitting HSH model...\n');
    hsh_order = params.hsh_order;
    intrp_hi = HSH_getHarmonics(hsh_order, dirs_fit, 'real');
    coef_HSH = LeastSquares(fit_data, intrp_hi);
    
    % Relight using specified indices
    fprintf('Generating relit images...\n');
    relight_LP = LP1(relight_indices, :);
    n_relight = length(relight_indices);
    
    % Initialize results structure
    channel_results = struct();
    channel_results.DMD = zeros(Images.height, Images.width, n_relight);
    channel_results.PTM = zeros(Images.height, Images.width, n_relight);
    channel_results.HSH = zeros(Images.height, Images.width, n_relight);
    
    % Generate relit images
    for i = 1:n_relight
        % DMD relighting
        cur_im_dirE = reshape(relight_LP(i,:), 1, 1, 3);
        interp_Qi = DMD_getEigenModes(data.nb_modes, cur_im_dirE, ...
            data.modal_basis(:,1:data.nb_modes), data.normal_elt, data.Up);
        image_dmd = interp_Qi * data.coef_DMD;
        channel_results.DMD(:,:,i) = reshape(image_dmd, Images.height, Images.width);
        
        % PTM relighting
        new_dir = LP_xyz2phitheta(relight_LP(i,:));
        interp_pi_new = PTM_terms(new_dir);
        image_ptm = interp_pi_new * coef_PTM;
        channel_results.PTM(:,:,i) = reshape(image_ptm, Images.height, Images.width);
        
        % HSH relighting
        intrp_hi_new = HSH_getHarmonics(hsh_order, new_dir, 'real');
        image_hsh = intrp_hi_new * coef_HSH;
        channel_results.HSH(:,:,i) = reshape(image_hsh, Images.height, Images.width);
    end
end

function save_results_and_compute_errors(results, Images, is_rgb, rti_dir, relight_indices, params)
    % Save results if requested
    if params.save_results
        save_all_results(results, rti_dir, relight_indices, is_rgb, Images);
    end
    
    % Compute errors if requested
    if params.compute_errors
        compute_and_display_errors(results, Images, relight_indices, is_rgb, params.save_results, rti_dir);
    end
end

function save_all_results(results, rti_dir, relight_indices, is_rgb, Images)
    fprintf('Saving results...\n');
    save_dir = fullfile(rti_dir, 'results');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    % Create method-specific directories
    dmd_dir = fullfile(save_dir, 'dmd');
    ptm_dir = fullfile(save_dir, 'ptm');
    hsh_dir = fullfile(save_dir, 'hsh');
    
    if ~exist(dmd_dir, 'dir'), mkdir(dmd_dir); end
    if ~exist(ptm_dir, 'dir'), mkdir(ptm_dir); end
    if ~exist(hsh_dir, 'dir'), mkdir(hsh_dir); end
    
    for i = 1:length(relight_indices)
        if is_rgb
            imwrite(results.DMD(:,:,:,i), fullfile(dmd_dir, sprintf('dmd_relit_%03d.png', relight_indices(i))));
            imwrite(results.PTM(:,:,:,i), fullfile(ptm_dir, sprintf('ptm_relit_%03d.png', relight_indices(i))));
            imwrite(results.HSH(:,:,:,i), fullfile(hsh_dir, sprintf('hsh_relit_%03d.png', relight_indices(i))));
        else
            imwrite(results.DMD(:,:,i), fullfile(dmd_dir, sprintf('dmd_relit_%03d.png', relight_indices(i))));
            imwrite(results.PTM(:,:,i), fullfile(ptm_dir, sprintf('ptm_relit_%03d.png', relight_indices(i))));
            imwrite(results.HSH(:,:,i), fullfile(hsh_dir, sprintf('hsh_relit_%03d.png', relight_indices(i))));
        end
    end
end

function compute_and_display_errors(results, Images, relight_indices, is_rgb, save_results, save_dir)
    fprintf('Computing errors...\n');
    n_relight = length(relight_indices);
    errors = struct();
    errors.DMD = zeros(1, n_relight);
    errors.PTM = zeros(1, n_relight);
    errors.HSH = zeros(1, n_relight);
    
    for i = 1:n_relight
        if is_rgb
            % Get original RGB channels
            orig_r = reshape(Images.DataV_R(relight_indices(i),:), [Images.height, Images.width]);
            orig_g = reshape(Images.DataV_G(relight_indices(i),:), [Images.height, Images.width]);
            orig_b = reshape(Images.DataV_B(relight_indices(i),:), [Images.height, Images.width]);
            
            % DMD error
            recon_dmd_r = results.DMD(:,:,1,i);
            recon_dmd_g = results.DMD(:,:,2,i);
            recon_dmd_b = results.DMD(:,:,3,i);
            errors.DMD(i) = (compute_rmse(orig_r, recon_dmd_r) + ...
                            compute_rmse(orig_g, recon_dmd_g) + ...
                            compute_rmse(orig_b, recon_dmd_b)) / 3;
            
            % PTM error
            recon_ptm_r = results.PTM(:,:,1,i);
            recon_ptm_g = results.PTM(:,:,2,i);
            recon_ptm_b = results.PTM(:,:,3,i);
            errors.PTM(i) = (compute_rmse(orig_r, recon_ptm_r) + ...
                            compute_rmse(orig_g, recon_ptm_g) + ...
                            compute_rmse(orig_b, recon_ptm_b)) / 3;
            
            % HSH error
            recon_hsh_r = results.HSH(:,:,1,i);
            recon_hsh_g = results.HSH(:,:,2,i);
            recon_hsh_b = results.HSH(:,:,3,i);
            errors.HSH(i) = (compute_rmse(orig_r, recon_hsh_r) + ...
                            compute_rmse(orig_g, recon_hsh_g) + ...
                            compute_rmse(orig_b, recon_hsh_b)) / 3;
        else
            % For grayscale images
            orig = reshape(Images.DataV(relight_indices(i),:), [Images.height, Images.width]);
            errors.DMD(i) = compute_rmse(orig, results.DMD(:,:,i));
            errors.PTM(i) = compute_rmse(orig, results.PTM(:,:,i));
            errors.HSH(i) = compute_rmse(orig, results.HSH(:,:,i));
        end
    end
    
    % Save errors if requested
    if save_results
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        error_file = fullfile(save_dir, 'errors.mat');
        save(error_file, 'errors');
    end
    
    % Display average errors
    fprintf('\nAverage RMSE:\n');
    fprintf('DMD: %.4f\n', mean(errors.DMD));
    fprintf('PTM: %.4f\n', mean(errors.PTM));
    fprintf('HSH: %.4f\n', mean(errors.HSH));
end

function rmse = compute_rmse(orig, recon)
    % Helper function to compute RMSE
    rmse = sqrt(mean((orig(:) - recon(:)).^2));
end

function str = conditional(condition, true_str, false_str)
    if condition
        str = true_str;
    else
        str = false_str;
    end
end