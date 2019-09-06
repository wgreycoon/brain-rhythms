function [ time_for_iteration_out, progress_bar_position ] = pbar( npermutes, permidx, time_for_iteration_in, progress_bar_position, thistic, headerstr )
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
%     npermutes                 :    outer loop max
%     permidx                   :    current loop index
%     time_for_iteration_in     :    time elapsed in last processing loop
%     progress_bar_position     :    last progress bar size
%     thistic                   :    timestamp from start of last
%                                    processing loop
%     headerstr                 :    what to call it (what are you
%                                    monitoring the progress of?)
%
% ---------------------------
%
% Output
% ---------------------------
%     time_for_iteration_out    :    time elapsed in this processing loop
%     thistic                   :    updated timestamp
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
