function [] = plot_MFIs(colonies, varargin)
%PLOT_MFIS Plots a histogram of MFIs
%   Arguments:
%       COLONIES : cell array of colonies; each colony should be a table, e.g.
%           produced by PARSE_MFI_LINEAGE
%   Keyword arguments:
%       'tframes' : array of frame numbers to use. If empty (default), all frames
%           from the first colony will be plotted
%       'Channels': cell array of channel names
%       'Colors': optional cell array mapping channel names to colors. By
%           default chooses reasonable colors for GFP, RFP, mCherry, GFP, and
%           YFP, or automatic colors for others
%       'ColonyNames': cell array of colony names, with same dimensions as
%           COLONIES. If absent, colonies will be numbered `1:length(COLONIES)`


p = inputParser;
st = dbstack; p.FunctionName = st.name;
p.StructExpand = false;

addOptional(p,'tframes',[])
addParameter(p,'Channels',{ 'GFP', 'RFP' })
addParameter(p,'Colors',...
        {'GFP','green';
        'RFP','red';
        'mCherry','red';
        'CFP','cyan';
        'YFP','yellow'});
addParameter(p,'ColonyNames', {});

parse(p,varargin{:})
args = p.Results;
tframes = args.tframes;
channels = args.Channels;
colors = args.Colors;
colony_names = args.ColonyNames;


% if unspecified, default colony names to numbers
n_colonies = length(colonies);
if isempty(colony_names)
   colony_names = cellstr(num2str(1:n_colonies));
end
figure;
    
% for each colony
for ii = 1:n_colonies
    colony = colonies{ii};
    
    % if unspecified, default to # of frames in first colony
    if isempty(tframes)
        tframes = unique(colony.tframe);
    end
    n_tframes = length(tframes);
    
    colony_subplots = [];
    
    % for each frame
    for jj = 1:n_tframes
        colony_subplots(jj) = subplot(length(tframes),length(colonies), (ii-1)*n_tframes+jj );
        title(sprintf('%s, t = %d', colony_names{:}, tframes(jj)));
        
        colony_frame = colony(colony.tframe == tframes(jj), :);
        
        % for each channel, overlay a histogram
        for c = 1:length(channels)
            if ~ismember(channels{c}, {colors{:,1}})
                color = 'auto';
            else
                color = colors{ismember({colors{:,1}},channels{c}),2};
            end
            histogram(colony_frame.(['MFI_' channels{c}]), 'FaceAlpha',0.5, 'FaceColor',color);
            hold on;
        end
        hold off;
        % hide X ticks except for the last row
        if jj < n_tframes; set(colony_subplots(jj),'xticklabel',[]); end
    end
    % show X axis label for last row only; link all X axes to same
    % dimensions
    xlabel('MFI')
    linkaxes(colony_subplots,'x')
end

legend(channels);

