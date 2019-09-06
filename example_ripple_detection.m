%========================================================================
%% Example RIPPLE detection
%========================================================================
%---------------------------
% This script demonstrates ripple detection using the "detect_ripples" 
% function, which is based on the validated detector in Staba et al. 2002
% and adds spike rejection from Helfrich et al. 2019.  The script also
% generates average time-frequency and time series representations of each
% example channels' avg. ripple, and individual detections if desired.
%
% ======
% Author(s):   William Coon   wcoon@mgh.harvard.edu; wgreycoon@gmail.com
%
% Last edited: September 4, 2019 1356hrs
% --------------------------------------------

load('example_ripple_data')
addpath('./cbrewer')
addpath('./functions')

%data specifics (included in example *.mat file...)
% channel_labels          = {'BLA2-BLA1','HR_8-HR_7'};
% fs                      = 512; % sample data is sampled at 512 Hz

%function parameters
ncyclesthresh           = 6;   % default from Staba et al. is '6'. see function for documentation
rip_ampl_factor         = 5.0; % default from Staba et al. is '5' see function for documentation
ripthreshs              = [];  % ignore, legacy...
verify                  = false; % legacy for visualizing detections...

%set ripple band parameters for FIR filtering here
params.lowstop          = 250; 
params.lowpassfreq      = 200;     
params.highpassfreq     = 80;
params.highstop         = 70;


% these are for additional spike rejection. 'artrejmode' can be any of the
% following strings:
%       'none'   :    no ratio-based detection performed
%       'high'   :    a high band (160-190Hz), defined in the
%                     detect_ripples function, is used for broadband
%                     artifact rejection
%       'low'    :    a low band (50-60 Hz), defined in the detect_ripples
%                     function, is used for broadband artifact rejection
%   'highlow'    :    both bands used
% I recommend using 'none', since it wasn't done in the Staba paper and in
% my experience didn't work very well for artifact rejection -- it threw
% out a lot of real detections along with spikes.
artrejmode              = 'none';

% 'artratio' is a scalar. If positive, this is the ratio of ripple power to
% artifact band power, defined in the 'artrejmode' variable above, that a
% detection must exceed to be counted as a ripple.  So for example if this
% was '3', then ripple power needs to be at least 3x greater than artifact
% band power to be counted as a detection. Setting this to '-1', as I did
% here, skips ratio-based spike rejection and instead applies the "3-peak
% rule" implemented in the Helfrich et al. 2019 paper (i.e., at the peak of
% a ripple candidate, the RAW trace must exhibit a ripple band oscillatory
% peak on either side of the detected peak, within a 40ms window -- so it
% is looking to find 3 peaks within 40ms centered on the peak of the ripple
artratio                = -1;

% detect!
[allRipples] = detect_ripples( double(ecog_signals)',...
                                    channel_labels,...
                                    fs,...
                                    ncyclesthresh,...
                                    rip_ampl_factor,...
                                    params,...
                                    verify,...
                                    ripthreshs,...
                                    artrejmode,...
                                    artratio);

%========================================================        



%%
%%

%% gather ripple time series, generate TF plots 

ripsavdir = './TF.ripples';
if ~exist(ripsavdir,'dir')
    mkdir(ripsavdir)
end

%Wavelet filtering RIPPLES
wavfreqs = exp(linspace(log(0.5),log(300),90));
wavfreqs = wavfreqs(1:84);

n        = 11;

%Getting some baseline estimates here...
pdall = wavbank( ecog_signals(1:fs*3000,:), wavfreqs, fs, n );

for hippch = 1:length(channel_labels)
    
    %========================================================================
    %% Gathering events (RIPPLES)
    twin = 1; %time window, in seconds, to gather spindle activity
    t = -(twin/2):1/fs:(twin/2);

    %bp signal
    ripsigbp = firbandpass2(fs, double(ecog_signals(:,hippch)), params.highpassfreq, params.lowpassfreq, params.highstop, params.lowstop);

    %collect coupled ripples
    riptracesPfilt = nan( twin*fs+1, allRipples(hippch).number );
    riptracesP     = nan( twin*fs+1, allRipples(hippch).number ); fprintf('[')
    for ripidx = 1:allRipples(hippch).number
        if allRipples(hippch).filtPeakIdx(ripidx)-(fs*twin/2) > 0 && allRipples(hippch).filtPeakIdx(ripidx)+(fs*twin/2) < size(ecog_signals,1)

        % NOTE: this locks to the PEAK of the ripple-filtered trace
        riptracesP(:,ripidx) = ecog_signals( allRipples(hippch).filtPeakIdx(ripidx)-(fs*twin/2):allRipples(hippch).filtPeakIdx(ripidx)+(fs*twin/2), hippch ); 
        riptracesPfilt(:,ripidx) = ripsigbp( allRipples(hippch).filtPeakIdx(ripidx)-(fs*twin/2):allRipples(hippch).filtPeakIdx(ripidx)+(fs*twin/2)); 

        end
    end
    %------------------------------------------------------
    
    %========================================================================
    % Wavelet filtering (RIPPLES)
    powerdata = wavbank( riptracesP, wavfreqs, fs, n );   

    %normalize 
    pd=powerdata;%powerdata=pd;
    for fidx = 1:size(powerdata,3)

        %get a baseline to norm against
        pow4norming = nanmean(reshape(pdall(:,hippch,fidx),[],1));

        %normalize 
        pd(:,:,fidx) = 10*log10(powerdata(:,:,fidx)./pow4norming);

    end
    %------------------------------------------------------------------------
    
    %========================================================================
    % Plotting

if size(riptracesPfilt,2)~=0    
    twin2plot = 0.5;
    t = -(twin2plot/2):1/fs:(twin2plot/2);
    fprintf('Creating Time-Frequency Plots ...\n')
    fprintf('[')
    figure('visible','on'); set(gcf,'position',[857          31         669        1062]); setPaperSize(gcf);    
    if ~exist(sprintf('%s/%s/indiv',ripsavdir,channel_labels{hippch}),'dir')
        mkdir(sprintf('%s/%s/indiv',ripsavdir,channel_labels{hippch}))
    end
%%% UNCOMMENT TO VISUALIZE INDIVIDUAL RIPPLES
%     for ripidx = 1:allRipples(hippch).number
%         
%         clf
%         
%         tidx1 = (twin/2)*fs - ((twin2plot/2)*fs);
%         tidx2 = (twin/2)*fs + ((twin2plot/2)*fs);
%         toplot  = squeeze(pd(tidx1:tidx2,ripidx,:));
%         
%         subplot(211)
%             contourf(t*1000,wavfreqs,(toplot)',80,'linecolor','none')
%             hold on
%             colormap(flipud(cbrewer('div','Spectral',256))); set(gca,'clim',[-40 6])
%             H=colorbar;
%             title(sprintf('Ripple %d/%d',ripidx,allRipples(hippch).number),'fontsize',16,'fontweight','bold')
%             ylabel('Frequency (Hz) ','fontsize',14)
%             xlabel('Time (ms)','fontsize',14)
%             ylabel(H,'10*log10(x/mean(bsl))','fontsize',12);
%             
%         subplot(212)
%             plot(t,riptracesPfilt(tidx1:tidx2,ripidx)); c=colorbar; set(c,'visible','off')
%             hold on; plot(t,abs(hilbert(riptracesPfilt(tidx1:tidx2,ripidx))),'r')
%              plot(t,riptracesP(tidx1:tidx2,ripidx),'g')
%              grid on;
%             ylabel('Amplitude (uV)','fontsize',14)
%             xlabel('Time (ms)','fontsize',14)
%                         
%             print(sprintf('%s/TF.ripples/%s/indiv/ripple_%03d.png',figsavedir,channel_labels{hippch},ripidx),'-dpng','-r72')
%         
%     end

    %GRAND AVERAGE OF RIPPLES
    clf

    tidx1 = (twin/2)*fs - ((twin2plot/2)*fs);
    tidx2 = (twin/2)*fs + ((twin2plot/2)*fs);
    toplot  = squeeze(mean( (pd(tidx1:tidx2,:,:)), 2 ));

    subplot(211)
        contourf(t*1000,wavfreqs,(toplot)',80,'linecolor','none')
        colormap(flipud(cbrewer('div','Spectral',256)))  
        colorbar
        H = colorbar;
        title(sprintf('%s Ripple Avg. (n=%d)',channel_labels{hippch},allRipples(hippch).number),'fontsize',16,'fontweight','bold')
        ylabel('Frequency (Hz) ','fontsize',14)
        xlabel('Time (ms)','fontsize',14)
        ylabel(H,'10*log10(x/mean(bsl))','fontsize',12);
        if ismac, set(gca,'fontsize',22), else set(gca,'fontsize',16), end

    subplot(212)
        x=plot(t,mean(riptracesPfilt(tidx1:tidx2,:),2)); c=colorbar; set(c,'visible','off'); set(x,'linewidth',3)
        hold on; x=plot(t,nanmean(abs(hilbert(riptracesPfilt(tidx1:tidx2,:))),2),'r'); set(x,'linewidth',3)
        grid on;
        x=plot(t,mean(riptracesP(tidx1:tidx2,:),2),'g'); set(x,'linewidth',3)
        ylabel('Amplitude (uV)','fontsize',14)
        xlabel('Time (s)','fontsize',14)
        set(gca,'xlim',[-twin2plot/2 twin2plot/2])
        if ismac, set(gca,'fontsize',22), else set(gca,'fontsize',16), end

        
    % print(sprintf('%sTF.ripples/_ripple_avg_%d_ripples.png',fdir,allRipples(hippch).number),'-dpng','-r72')
    print(sprintf('%s/%d_%s_broadband_ripple_avg_%d_ripples.png',ripsavdir,hippch,channel_labels{hippch},allRipples(hippch).number),'-dpng','-r150')
    print(sprintf('%s/%d_%s_broadband_ripple_avg_%d_ripples.pdf',ripsavdir,hippch,channel_labels{hippch},allRipples(hippch).number),'-dpdf','-r150','-painters')
    close gcf
    
end %end check that there is data to process here...


end
%% 
