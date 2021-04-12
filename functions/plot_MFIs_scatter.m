function [] = plot_MFIs_scatter(colonies, xchannel, ychannel, varargin)
%PLOT_MFIS_SCATTER Makes a scatterplot of MFIs
%   Arguments:
%       COLONIES : cell array of colonies; each colony should be a table, e.g.
%           produced by PARSE_MFI_LINEAGE
%       XCHANNEL : channel to plot on the x axis
%       YCHANNEL : channel to plot on the y axis
%   Keyword arguments:
%       'tframes' : array of frame numbers to use. If empty (default), all frames
%           from the first colony will be plotted
%       'ColonyNames': cell array of colony names, with same dimensions as
%           COLONIES. If absent, colonies will be numbered `1:length(COLONIES)`


p = inputParser;
st = dbstack; p.FunctionName = st.name;
p.StructExpand = false;

addOptional(p,'tframes',[])
addParameter(p,'Colors',...
        {'GFP','green';
        'RFP','red';
        'mCherry','red';
        'CFP','cyan';
        'YFP','yellow'});
addParameter(p,'ColonyNames', {});
addParameter(p,'FrameDuration', []);

parse(p,varargin{:})
args = p.Results;
tframes = args.tframes;
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
    
    if ~isempty(args.FrameDuration)
        durations = tframes .* args.FrameDuration;
    else
        durations = [];
    end
    
    % for each frame
    for jj = 1:n_tframes
        colony_subplots(jj) = subplot(length(tframes),length(colonies), (ii-1)*n_tframes+jj );
        if ~isempty(colony_names)
            if isempty(durations)
                title(sprintf('%s, t = frame %d', colony_names{:}, tframes(jj)));
            else
                title(sprintf('%s, t = %d', colony_names{:}, durations(jj)));
            end
        else
            if isempty(durations)
                title(sprintf('t = frame %d', tframes(jj)));
            else
                title(sprintf('t = %d', durations(jj)));
            end
        end
        
        colony_frame = colony(colony.tframe == tframes(jj), :);
        plot(colony_frame.(['MFI_' xchannel]),colony_frame.(['MFI_' ychannel]),'.')

        % hide X ticks except for the last row
        if jj < n_tframes; set(colony_subplots(jj),'xticklabel',[]); end
        ylabel(ychannel)
    end
    % show X axis label for last row only; link all X axes to same
    % dimensions
    xlabel(xchannel)
    linkaxes(colony_subplots,'x')
end

% legend(channels);

