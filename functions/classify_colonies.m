function [onoff, fon] = classify_colonies(directory, classinfo, tframes, varargin)
% CLASSIFY_COLONIES  classifies cells in colonies as ON or OFF according to
% a threshold set using control samples.
%
% Arguments:
%    DIRECTORY : Folder containing .CSV files, one per microcolony, named
%        according to
%        `CSVFilePattern`
%    CLASSINFO : matrix or table from CALCULATE_THRESHOLD
%    TFRAMES : list of frames to calculate (optional; if ommitted, all
%        frames found in both the colony and `classinfo` table will be used)
%
% Keyword arguments:
%     'CSVFilePattern' : How do you want your CSV files to be named?
%                        Default: '^(?<image>.*)\.csv' (matches all .csv
%                        files)
%     'Channel' : Which channel to use from the CSV file? The column
%                 ['MFI_' Channel] from the CSV file will be used to get MFIs
%     'ShowPlots' : true to show some summary plots (default true)
%
% Returns:
%     ONOFF : NxMx2 matrix, where ONOFF(i,j,1) is the number of OFF cells 
%             in colony i, timepoint j, and ONOFF(i,j,2) is the number of 
%             ON cells
%       FON : NxM matrix, where FON(i,j) is the fraction of ON cells in
%             colony i, timepoint j

p = inputParser;
st = dbstack;
p.FunctionName = st.name; %'parse_MFI_lineage';
p.StructExpand = false;

addOptional(p,'tframes', []);
addParameter(p,'CSVFilePattern', '^(?<image>.*)\.csv');
addParameter(p,'Channel', 'GFP');
addParameter(p,'ShowPlots', true);

parse(p,tframes,varargin{:})
args = p.Results;
CSV_FILE_PATTERN = args.CSVFilePattern;
CHANNEL = args.Channel;
show_plots = args.ShowPlots;
tframes = args.tframes;

%directory containing original mat files
%directory ='F:\Dropbox\Christina_data\PA14 WT Pt-sfGFP in NTA\combined_4.4.19_4.14.19\'
csvfileobj=dir(fullfile(directory, '*.csv'));
[csvfilenames{1:length(csvfileobj)}]=csvfileobj(:).name;
numfiles=length(csvfileobj);

if istable(classinfo)
    % use instead of table2array to make sure columns are in right order
    classinfo = [classinfo.frame classinfo.threshold classinfo.LCI classinfo.UCI];
end

if isempty(tframes)
    class_tframes = classinfo(:,1);
end

if show_plots
    summary_figure = figure;
end

% loop through CSV files for each microcolony
for fnum=1:numfiles
    file_name_parts = regexp(csvfilenames{fnum}, CSV_FILE_PATTERN, 'names');
    csvfilename = fullfile(directory, csvfilenames{fnum});
    
    if (isempty(file_name_parts))
        fprintf('Skipping .csv file "%s" which does not fit the expected naming pattern\n', fnames{t})
        continue
    end
    colony_name = file_name_parts.image;
    
    data = readtable(csvfilename, detectImportOptions(csvfilename));
    data_mfis = data.(['MFI_' CHANNEL]);
    if isempty(tframes)
        colony_tframes = unique(data.tframe);
        colony_tframes = intersect(colony_tframes, class_tframes);
    else
        colony_tframes = tframes;
    end
    
    for t=1:length(colony_tframes)
        timepoint = colony_tframes(t);
        
        % Isolate GFP MFIs for this frame and classify each cell as
        % ON or OFF
        mfis{fnum,t} = data_mfis(data.tframe == timepoint);
        meanint(fnum,t) = nanmean(mfis{fnum,t});
        cellnum(fnum,t) = sum(isnan(mfis{fnum,t})==0);
        threshtmp=classinfo(find(classinfo(:,1)==colony_tframes(t)),2);
        
        % if the threshold classification dataset doesn't go far enough,
        % use the last available timepoint
        if isempty(threshtmp)==1
            threshtmp=classinfo(find(classinfo(:,1)==max(classinfo(:,1))),2);
        end
        thresh(fnum,t) = threshtmp;
        numon  = sum(sum(mfis{fnum,t} > thresh(t)));
        numoff = sum(sum(mfis{fnum,t} < thresh(t)));
        fon(fnum,t)=numon./(numon+numoff);
        
        onoff(fnum, t, 1) = numon;
        onoff(fnum, t, 2) = numoff;
    end
    
    
    if show_plots
        figure(summary_figure);
        subplot(3,1,1);
        hold on
        [F, x] = ecdf(mfis{fnum,:});

        %plot CDF of each colony
        plot(x,F,['-.'],'LineWidth',2)
        if fnum==1
            plot(thresh(fnum,end)*ones(size(F)),F,'k','LineWidth',2)
            hold on
        end
        xlabel([CHANNEL ' MFI'])
        ylabel('Empirical CDF')
        hold off


        %plot number on vs time
        subplot(3,1,3);
        hold on
        plot(colony_tframes, fon(fnum,:),'.-')
        xlabel('frame')
        ylabel('f_{on}')
        ylim([-.1 1.1])
        hold off

    end
    
end

if show_plots
    %plot distribution of fraction on
    subplot(3,1,2);
    hold on;
    histogram(fon(:),'BinMethod','sqrt','FaceColor','k','FaceAlpha',0.5)
    xlabel('fraction on')
    ylabel('Number colonies')
    hold off;
end




end