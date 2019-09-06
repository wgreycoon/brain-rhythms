function [ zm, zmc, real_t, perm_t, mci  ] = tfrstats( tfr, condition_mapping, n_permutes, voxel_pval, mcc_cluster_pval, nantoggle )
%% time-frequency response (TFR) stats w/ cluster-based correction for multiple comparisons 
% 
% -------------------------------------------------------------------------
% Cluster-based statistics for time-frequency responses (single-channel
% analysis). Given a set of TASK and BASELINE responses, calculate the
% statistical significance of changes relative to baseline by Monte
% Carlo-style permutations/label-shuffling. Task/Bsln labels are shuffled
% to generate distributions of TFRs. 
%
% this function is based on code from Mike X. Cohen's delightful book on
% "Analyzing Neural Time Series Data".  I highly recommend it.
%
% progress bar inspired by code from:
% www.mathworks.com/matlabcentral/fileexchange/4863-ascii-progress-bar
% ------------------------------------------------------------------------- 
%
% Dependencies
% ---------------------------
% functions:
%
% ---------------------------
%
% Input 
% ---------------------------
%     tfr                :    array of data (FREQ. x TIME x TRIAL)
%     condition_mapping  :    LOGICAL vector with indexes for TASK (1's) 
%                             and BSLN (0's)
%     n_permutes         :    Number of label shuffles to do
%     voxel_pval         :    threshold for first-pass (uncorrected) stats
%     mcc_cluster_pval   :    threshold for multiple comparisons 
%                             (cluster-corrected) stats
%     nantoggle          :    choose NaN-robust or not. Function is almost
%                             2x faster if using mean/std instead of 
%                             nanmean/nanstd
%
%     Example = tfrstats( tfr, condition_mapping, 1000, 0.05, 0.05 );
%
% ---------------------------
%
% Output
% ---------------------------
%     zm                 :    zmap with uncorrected stats
%     zmc                :    zmap with cluster-corrected stats
%     real_t             :    actual tmap from task v. baseline data
%
% ---------------------------
%
% Author(s):   William Coon   wcoon@mgh.harvard.edu
%
% Last edited: August 22, 2019 1339hrs
%
% -------------------------------------------------------------------------
%%

% compute actual t-test of difference (using unequal N and std)
if nantoggle
    tnum   = squeeze(nanmean(tfr(:,:,condition_mapping==1),3) - nanmean(tfr(:,:,condition_mapping==0),3));
    tdenom = sqrt( (nanstd(tfr(:,:,condition_mapping==1),0,3).^2)./sum(condition_mapping==1) + (nanstd(tfr(:,:,condition_mapping==0),0,3).^2)./sum(condition_mapping==0) );
else
    tnum   = squeeze(mean(tfr(:,:,condition_mapping==1),3) - mean(tfr(:,:,condition_mapping==0),3));
    tdenom = sqrt( (std(tfr(:,:,condition_mapping==1),0,3).^2)./sum(condition_mapping==1) + (std(tfr(:,:,condition_mapping==0),0,3).^2)./sum(condition_mapping==0) );
end

real_t = tnum./tdenom;

num_frex                    = size(tfr,1);
nTimepoints                 = size(tfr,2);

% initialize null hypothesis matrices
permuted_tvals              = zeros(n_permutes,num_frex,nTimepoints);
max_pixel_pvals             = zeros(n_permutes,2);
max_clust_info              = zeros(n_permutes,1);

%for status bar...
prog_bar_pos                = 0;
tfti                        = 0.01;

% generate pixel-specific null hypothesis parameter distributions
for permi = 1:n_permutes
    
    thistic                     = tic; %for time estimation in pbar fn.
    
    fake_condition_mapping      = sign(randn(size(tfr,3),1));
    fake_condition_mapping(fake_condition_mapping==-1) = 0;
    
    % compute t-map of null hypothesis
    if nantoggle
        tnum                    = squeeze(nanmean(tfr(:,:,fake_condition_mapping==1),3)-nanmean(tfr(:,:,fake_condition_mapping==0),3));
        tdenom                  = sqrt( (nanstd(tfr(:,:,fake_condition_mapping==1),0,3).^2)./sum(fake_condition_mapping==1) + (nanstd(tfr(:,:,fake_condition_mapping==0),0,3).^2)./sum(fake_condition_mapping==0) );
    else
        tnum                    = squeeze(mean(tfr(:,:,fake_condition_mapping==1),3)-mean(tfr(:,:,fake_condition_mapping==0),3));
        tdenom                  = sqrt( (std(tfr(:,:,fake_condition_mapping==1),0,3).^2)./sum(fake_condition_mapping==1) + (std(tfr(:,:,fake_condition_mapping==0),0,3).^2)./sum(fake_condition_mapping==0) );
    end
    tmap                        = tnum./tdenom;
    
    % save all permuted values
    permuted_tvals(permi,:,:)   = tmap;
    
    % save maximum pixel values
    max_pixel_pvals(permi,:)    = [ min(tmap(:)) max(tmap(:)) ];
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % note that here, clusters were obtained by parametrically thresholding
    % the t-maps
    tmap(abs(tmap)<tinv(1-voxel_pval,size(tfr,3)-1))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo                   = bwconncomp(tmap);
    
    max_clust_info(permi)       = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero accounts for empty maps
    
    %display status
    [tfti,prog_bar_pos]         = pbarlocal( n_permutes, permi, tfti, prog_bar_pos, thistic, 'CLUSTER STATS');    
         
end

% now compute Z-map
zmap = (real_t-squeeze(mean(permuted_tvals,1)))./squeeze(std(permuted_tvals));

% uncorrected stats contours...
zmapthreshuncorrected = zmap;
zmapthreshuncorrected(abs(zmapthreshuncorrected)<norminv(1-voxel_pval))=false;
zmapthreshuncorrected=logical(zmapthreshuncorrected);

% apply cluster-level corrected threshold
zmapthresh = zmap;

% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;

% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-mcc_cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end

zm      = zmapthreshuncorrected;
zmc     = zmapthresh;
mci     = max_clust_info;
perm_t  = permuted_tvals;
% end



    function [ time_for_iteration_out, progress_bar_position ] = pbarlocal( npermutes, permidx, time_for_iteration_in, progress_bar_position, thistic, headerstr )
    %% ASCII Progress Bar w/ Estimated time-to-completion 
    % 
    % -------------------------------------------------------------------------
    % Function for returning a progress bar, based on a defined loop size
    % ------------------------------------------------------------------------- 
    %
    % Dependencies
    % ---------------------------
    % functions:
    %
    % ---------------------------
    %
    % Input 
    % ---------------------------
    %     npermutes                  :    outer loop max
    %     permidx                    :    current loop index
    %     time_for_iteration_in      :    time elapsed in last processing loop
    %     progress_bar_position      :    last progress bar size
    %     thistic                    :    timestamp from start of last
    %                                     processing loop
    %     headerstr                  :    what to call it (what are you
    %                                    monitoring the progress of?)
    %
    % ---------------------------
    %
    % Output
    % ---------------------------
    %     time_for_iteration_out     :    time elapsed in this processing loop
    %     thistic                    :    updated timestamp
    %
    % ---------------------------
    %
    % Author(s):   William Coon   wcoon@mgh.harvard.edu
    %
    % Last edited: August 23, 2019 1904hrs
    %
    % -------------------------------------------------------------------------
    %%

    progress_bar_position = progress_bar_position + 1 / npermutes;
    clc;
    fprintf('|=================== %s ====================|\n',headerstr);
    progress_string='|';       
    for counter = 1:floor(progress_bar_position * 100 / 2),
       progress_string = [progress_string, '#'];
    end
    disp(progress_string);
    % disp(['|================== ',num2str(floor(progress_bar_position * 100)),'% completed ====================|']);
    fprintf('|=================== %02d%% completed ====================|\n',floor(progress_bar_position * 100));
    steps_remaining = npermutes - permidx;
    minutes = floor(time_for_iteration_in * steps_remaining / 60);
    seconds = rem(floor(time_for_iteration_in *  steps_remaining), 60);
    disp(' ');
    if (seconds > 9),
     disp(['            Estimated remaining time: ', num2str(minutes), ':', num2str(seconds)]);
    else
     disp(['            Estimated remaining time: ', num2str(minutes), ':0', num2str(seconds)]);
    end
    time_for_iteration_out = toc(thistic);


    end

end

