function setPaperSize(fig_handle)     

    % Set figure and paper to use the same unit
    set(fig_handle, 'Units', 'centimeters')
    set(fig_handle, 'PaperUnits','centimeters');

    % Position of figure is of form [left bottom width height]
    % We only care about width and height
    pos = get(fig_handle,'Position');

    % Set paper size to be same as figure size
    set(fig_handle, 'PaperSize', [pos(3) pos(4)]);

    % Set figure to start at bottom left of paper
    % This ensures that figure and paper will match up in size
    set(fig_handle, 'PaperPositionMode', 'manual');
    set(fig_handle, 'PaperPosition', [0 0 pos(3) pos(4)]);  
    
end