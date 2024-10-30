function processRTISelective(rti_dir, fit_indices, relight_indices, params)
    % Import required packages
    import pkg_DMD.* 
    import pkg_PTM.* 
    import pkg_HSH.* 
    import pkg_fcns.*

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
    
    % Get image paths and light positions
    image_paths = temp_L{1}(2:end);  % First entry is often header
    
    % Initialize Images structure
    Images = struct();
    Images.LP1 = [temp_L{2}(2:end) temp_L{3}(2:end) temp_L{4}(2:end)];
    nfiles = length(image_paths);
    fprintf('Found %d images in LP file\n', nfiles);
    
    % Read first image to get dimensions
    first_img_path = fullfile(rti_dir, image_paths{1});
    if ~exist(first_img_path, 'file')
        error('Cannot find first image: %s', first_img_path);
    end
    first_img = imread(first_img_path);
    Images.height = size(first_img,1);
    Images.width = size(first_img,2);
    fprintf('Image dimensions: %dx%d\n', Images.height, Images.width);
    
    % Load and preprocess images
    Images.DataV = zeros(nfiles, Images.height * Images.width, 'single');
    fprintf('Loading and preprocessing images...\n');
    for ii = 1:nfiles
        img_path = fullfile(rti_dir, image_paths{ii});
        if ~exist(img_path, 'file')
            error('Cannot find image: %s', img_path);
        end
        currentimage = imread(img_path);
        if size(currentimage,3) > 1
            currentimage = rgb2gray(currentimage);
        end
        Images.DataV(ii,:) = single(currentimage(:))' / 255;
        if mod(ii,10) == 0
            fprintf('Processed %d/%d images\n', ii, nfiles);
        end
    end
    
    % Validate indices
    if max(fit_indices) > nfiles || max(relight_indices) > nfiles
        error('Indices exceed number of available images');
    end
    
    % Extract fitting subset
    fit_data = Images.DataV(fit_indices, :);
    fit_LP = Images.LP1(fit_indices, :);
    fprintf('Using %d images for fitting\n', length(fit_indices));
    
    % Fit models using subset
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
    relight_LP = Images.LP1(relight_indices, :);

%     fprintf('Using %d images for relighting\n', length(relight_indices));
    
    n_relight = length(relight_indices);
    
    % Preallocate output structures
    relit_images = struct();
    relit_images.DMD = zeros(Images.height, Images.width, n_relight);
    relit_images.PTM = zeros(Images.height, Images.width, n_relight);
    relit_images.HSH = zeros(Images.height, Images.width, n_relight);
    
    % Generate relit images
    for i = 1:n_relight
        fprintf('Relighting image %d/%d\n', i, n_relight);
        
        % DMD relighting
        cur_im_dirE = reshape(relight_LP(i,:), 1, 1, 3);
        interp_Qi = DMD_getEigenModes(data.nb_modes, cur_im_dirE, ...
            data.modal_basis(:,1:data.nb_modes), data.normal_elt, data.Up);
        image_dmd = interp_Qi * data.coef_DMD;
        relit_images.DMD(:,:,i) = reshape(image_dmd, Images.height, Images.width);
        
        % PTM relighting
        new_dir = LP_xyz2phitheta(relight_LP(i,:));
        interp_pi_new = PTM_terms(new_dir);
        image_ptm = interp_pi_new * coef_PTM;
        relit_images.PTM(:,:,i) = reshape(image_ptm, Images.height, Images.width);
        
        % HSH relighting
        intrp_hi_new = HSH_getHarmonics(hsh_order, new_dir, 'real');
        image_hsh = intrp_hi_new * coef_HSH;
        relit_images.HSH(:,:,i) = reshape(image_hsh, Images.height, Images.width);
    end
    
    % Save results if requested
    if params.save_results
        fprintf('Saving results...\n');
        save_dir = fullfile(rti_dir, 'results');
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        
        % Save coefficients
        coef_file = fullfile(save_dir, 'coefficients.mat');
        save(coef_file, 'data', 'coef_PTM', 'coef_HSH');
        
        for i = 1:n_relight
            % Create method-specific directories
            dmd_dir = fullfile(save_dir, 'dmd');
            ptm_dir = fullfile(save_dir, 'ptm');
            hsh_dir = fullfile(save_dir, 'hsh');
            
            if ~exist(dmd_dir, 'dir'), mkdir(dmd_dir); end
            if ~exist(ptm_dir, 'dir'), mkdir(ptm_dir); end
            if ~exist(hsh_dir, 'dir'), mkdir(hsh_dir); end
            
            % Save DMD images
            imwrite(relit_images.DMD(:,:,i), ...
                fullfile(dmd_dir, sprintf('dmd_relit_%03d.png', relight_indices(i))));
            
            % Save PTM images
            imwrite(relit_images.PTM(:,:,i), ...
                fullfile(ptm_dir, sprintf('ptm_relit_%03d.png', relight_indices(i))));
            
            % Save HSH images
            imwrite(relit_images.HSH(:,:,i), ...
                fullfile(hsh_dir, sprintf('hsh_relit_%03d.png', relight_indices(i))));
        end
    end
    % In the error calculation section, replace with:
    if params.compute_errors
        fprintf('Computing errors...\n');
        original_images = Images.DataV(relight_indices, :);
        errors = struct();
        errors.DMD = zeros(1, n_relight);
        errors.PTM = zeros(1, n_relight);
        errors.HSH = zeros(1, n_relight);
        
        for i = 1:n_relight
            orig = reshape(original_images(i,:), [Images.height, Images.width]);
            relit_dmd = relit_images.DMD(:,:,i);
            relit_ptm = relit_images.PTM(:,:,i);
            relit_hsh = relit_images.HSH(:,:,i);
            
            % Ensure all images are in column vectors for comparison
            orig_vec = orig(:);
            dmd_vec = relit_dmd(:);
            ptm_vec = relit_ptm(:);
            hsh_vec = relit_hsh(:);
            
            % Calculate RMSE
            errors.DMD(i) = sqrt(mean((orig_vec - dmd_vec).^2));
            errors.PTM(i) = sqrt(mean((orig_vec - ptm_vec).^2));
            errors.HSH(i) = sqrt(mean((orig_vec - hsh_vec).^2));
        end
        
        % Save errors if results saving is enabled
        if params.save_results
            error_file = fullfile(save_dir, 'errors.mat');
            save(error_file, 'errors');
        end
        
        % Display average errors
        fprintf('\nAverage RMSE:\n');
        fprintf('DMD: %.4f\n', mean(errors.DMD));
        fprintf('PTM: %.4f\n', mean(errors.PTM));
        fprintf('HSH: %.4f\n', mean(errors.HSH));
    end
    
end