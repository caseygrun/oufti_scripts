function [outputs] = parse_MFI_lineage(directory, output_directory, varargin)
%PARSE_MFI_LINEAGE  use Oufti cell mesh and lineage data to export MFI data
%    T = PARSE_MFI_LINEAGE(DIRECTORY,OUTPUT_DIRECTORY,...) reads cell meshes
%        in .MAT files and extracts MFIs from .TIF files in DIRECTORY, then
%        exports MFI data to .CSV and .DAT files in OUTPUT_DIRECTORY. 
%        - The CSV file, named according to `OUT_FILE_PATTERN`, contains cell lineage
%          and MFI data. Each row is a timepoint for a particular cell, and columns
%          contain the following information: 
%            tframe : frame number
%            Cell : cell number from Oufti
%            Last_Ancestor : cell number of most recent ancestor of this cell
%            LeafStatus : 1 if this cell is a leaf (does not have any descendents)
%            MFI_{channel1} : mean fluorescence intensity of the cell in this
%              channel, after subtracting the average of the average intensity of
%              the background (parts of the image with no cells). This value is NOT
%              normalized and thus values depend on the range of intensities in the
%              input image. 
%            MFI_{channel2} : additional columns as necessary for other channels.
%        - The DAT files, named according to `DAT_FILE_PATTERN`, contain only MFI data
%          --- one .DAT file is created per channel per frame. 
% 
%    Various optional arguments can control the behavior of this function: 
%           'Channels'   which channels will you use? Your image files must be named using these
%                        channel names (case sensitive). You can have more than two channels. 
%                        Default: {'GFP', 'RFP'}
%     'MatFilePattern'   cell mesh and lineage data should be in an Oufti-produced .mat file. 
%                        how are your `.mat` files named? This should be a regular expression
%                        containing at least `channel` and `image` tokens/capture groups. See
%                        documentation for the `regexp` function. 
%                        Default: '^(?<image>.*)\.(?<ext>\w+)'
%                        (This particular pattern expects files like this:)
%                            some_long_experiment_name.mat
%     'TifFilePattern'   Phase and fluorescent channel images should be in separate (per-channel) 
%                        multi-page TIF files. 
%                        How are your `.tif` image files named? should be a format specifier
%                        suitable for `sprintf`. `%{channel}s` is the channel, `%{image}s` is the 
%                        base image name. 
%                        Default: '%{channel}s-%{image}s.tif'
%                        (This particular pattern expects files like this:
%                           GFP-some_long_experiment_name.tif
%                        where channel = 'GFP' and image =
%                        'some_long_experiment_name'.
%    'CSVFilePattern'    One CSV file will be produced per colony (e.g. per .MAT file). 
%                        How do you want your CSV file(s) named? should be a format string
%                        suitable for `sprintf` with placeholder `%{image}s` for the base 
%                        image name. 
%                        Default: '%{image}s.csv'
%    'DatFilePattern'    One DAT file will be produced per colony per channel per frame.
%                        How do you want your .DAT file(s) named? should be a format string
%                        suitable for `sprintf` with placeholders: `%{image}s` for the base 
%                        image name, %{channel} is the channel, and `%{frame}` will be the frame
%                        number. 
%                        Default: '%{image}s-%{channel}s-%{frame}d.csv'
%     
%    'FramesToExport'    which frames do you want to export? if this set to the empty array `[]`, all frames
%                        with 1+ cells will be exported. Note: the `Last_Ancestor` column in the CSV
%                        file may refer to cells birthed on frames outside this range. 
% 
%   Returns an Nx2 cell array `outputs`, where N = number of .MAT files.
%       outputs{i, 1} is the image base name; for example, for mesh file i named 'Ph-some_long_experiment_name.mat', 
%                     `outputs{i, 1} = 'some_long_experiment_name'`
%       outputs{i, 2} is the CSV file output for file i, as a MATLAB table. Each row is a single cell at a single time. 
%                     For example,
%                     `outputs{i, 1}.tframe`  is a vector indicating which frame number the row refers to
%                     `outputs{i, 1}.Cell`    is a vector indicating which cell number the row refers to
%                     `outputs{i, 1}.MFI_GFP` is a vector containing the GFP intensities of each cell at each frame
%                     You could plot a distribution of GFP MFIs from a single frame `t` from colony `i` like this:
%                         plot_gaussian(outputs{i, 2}.MFI_GFP(outputs{i, 2}.tframe == t))
% 
%   Daniel Lee, Christina Lin, and Casey Grun
%   MIT License
% 

p = inputParser;
st = dbstack;
p.FunctionName = st.name; %'parse_MFI_lineage';
p.StructExpand = false;

addParameter(p,'Channels',{ 'GFP', 'RFP' })
addParameter(p,'MatFilePattern', '^(?<image>.*)\.(?<ext>\w+)');
addParameter(p,'TifFilePattern', '%{channel}s-%{image}s.tif');
addParameter(p,'CSVFilePattern', '%{image}s.csv')
addParameter(p,'DatFilePattern', '%{image}s-%{channel}s-%{frame}d.dat')
addParameter(p,'FramesToExport', []);

parse(p,varargin{:})
args = p.Results;
channels = args.Channels;
OUT_FILE_PATTERN = args.CSVFilePattern;
DAT_FILE_PATTERN = args.DatFilePattern;
MAT_FILE_PATTERN = args.MatFilePattern;
TIF_FILE_PATTERN = args.TifFilePattern;
frames_to_export = args.FramesToExport;


% replace human-readable {placeholder} names with positional placeholders
% since MATLAB sprintf does not support named placeholders. 
TIF_FILE_PATTERN = strrep(TIF_FILE_PATTERN, '{channel}','1$'); 
TIF_FILE_PATTERN = strrep(TIF_FILE_PATTERN, '{image}','2$');
OUT_FILE_PATTERN = strrep(OUT_FILE_PATTERN, '{image}','1$');
DAT_FILE_PATTERN = strrep(DAT_FILE_PATTERN, '{image}','1$');
DAT_FILE_PATTERN = strrep(DAT_FILE_PATTERN, '{channel}','3$');
DAT_FILE_PATTERN = strrep(DAT_FILE_PATTERN, '{frame}','2$');

% create output directory
mkdir(output_directory)

% Load Oufti output .mat files from directory
matfileobj=dir(fullfile(directory, '*.mat'));
[matfilenames{1:length(matfileobj)}]=matfileobj(:).name;
numfiles=length(matfileobj);

outputs = cell(numfiles, 2);

% loop through colony files
for fnum=1:numfiles
    
    % check that this file matches the expected file name pattern, and 
    % extract the base image name. If this is not a Phase image, skip
    file_name_parts = regexp(matfilenames{fnum}, MAT_FILE_PATTERN, 'names');
    if (isempty(file_name_parts))
        fprintf('Skipping .mat file "%s" which does not fit the expected naming pattern\n', matfilenames{fnum})
        continue
    else 
        fprintf('Processing .mat file "%s"...\n', matfilenames{fnum})
    end
    outputs{fnum,1} = file_name_parts.image;
    
    % get complete path to file, extract base image name
    matfilename = fullfile(directory, matfilenames{fnum});
    img_base_name = file_name_parts.image; 
    
    % get path to TIF file for each channel
    imgfilenames = cellfun(...
        @(channel) fullfile(directory, ...
            sprintf(TIF_FILE_PATTERN,channel,img_base_name)),...
        channels, 'UniformOutput', false);
    fprintf('- Looking for images %s\n',strjoin(imgfilenames,','))
    
    
    % imgfilestr=['GFP-' matfilenames{fnum}(1:end-7) '.tif'];
    % rfpimgfilestr=['RFP-' matfilenames{fnum}(1:end-7) '.tif'];
    % imgfilename=[directory imgfilestr];
    % rfpfilename=[directory rfpimgfilestr];
    
    % load segmented meshes and other data from oufti
    load(matfilename);
    
    % which frames should we export? 
    % number of timepoints in movie
    tspan=length(cellListN);
    if (isempty(frames_to_export)) 
        % export all frames that have >0 cells
        tframes = find(cellListN>0); 
    else 
        tframes = frames_to_export; 
    end
    
    % for each frame
    for t = 1:length(tframes)
        imgs = cellfun(@(channel_image_file) ...
            double(imread(channel_image_file,tframes(t))), ...
            imgfilenames, 'UniformOutput', false); 
        
%         img=double(imread(imgfilename,tframes(t)));
%         rfpimg=double(imread(rfpfilename,tframes(t)));
        allIDlist=cellList.cellId;
        
        % Load CellIDs and segmentation mesh data.
        numcells=max(cell2mat(allIDlist(~cellfun('isempty',allIDlist))));
        cellIDs=cellList.cellId{t};
        meshinfo=cellList.meshData{t};
        
        % this mask will keep track of which cell occupies each pixel in
        % the image; wholemask(i, j) = cell ID of the cell at pixel (i, j),
        % or 0 if no cell at this position
        wholemask=zeros(size(imgs{1}));
        
        % cell array to hold MFIs for each channel
        % frame_MFIs{ii} is a 1 x (# of cells in frame) matrix
        frame_MFIs = cellfun(@(~) zeros(1,cellListN(tframes(t))), channels, ...
            'UniformOutput', false);
%         gfpmfi=zeros(1,cellListN(tframes(t)));
%         rfpmfi=zeros(1,cellListN(tframes(t)));

        % setup matrices to capture linage data
        ancestorlist = zeros(1,cellListN(tframes(t)));
        leafstatus = zeros(1,cellListN(tframes(t)));
        
        % for each cell
        for cellnum = 1:cellListN(tframes(t))
            meshcoord = meshinfo{cellnum}.mesh;
            if meshcoord == 0
                
            else
                
                % generate a mask identifying this cell, from its mesh
                borders=[flip(meshcoord(:,4)) flip(meshcoord(:,3)); meshcoord(:,2) meshcoord(:,1)];
                
                % remove points containing +/-infinity
                infrows = find(ismember(abs(borders),[Inf Inf],'rows')==1);
                
                if isempty(infrows)==0
                    borders(infrows,:)=[];
                end
                
                % create mask identifing the position of this cell
                bactmask=poly2mask(...
                    (double(borders(:,2))),...
                    (double(borders(:,1))),...
                    size(imgs{1},1),...
                    size(imgs{1},2));
                
                % for each channel, mask the channel image to only include
                % this cell, then take the mean value (MFI)
                for ii = 1:length(channels)
                    frame_MFIs{ii}(cellnum) = mean(imgs{ii}(bactmask==1));
                end

                % keep track of which cell occupies this pixel
                wholemask = wholemask + cellnum * bactmask;
                
                birthframe(cellIDs(cellnum))=meshinfo{cellnum}.birthframe;
                polarity(cellIDs(cellnum))=meshinfo{cellnum}.polarity;
                ancestors{cellIDs(cellnum)}=meshinfo{cellnum}.ancestors;
                descendants{cellIDs(cellnum)}=meshinfo{cellnum}.descendants;
                
                if isempty(ancestors{cellIDs(cellnum)})==0
                    ancestorlist(cellnum)=max(ancestors{cellIDs(cellnum)});
                else
                    ancestorlist(cellnum)=NaN;
                end
            end
        end
        
        % create a mask to identify which pixels contain any cells
        % (i.e. `wholemaskBWbin(i, j)` == 1 if pixel `i, j` contains a
        % cell, otherwise 0
        wholemaskBWbin=imdilate(wholemask>0,strel('disk',10));
        
        % subtract background fluorescence values for each channel
        for ii = 1:length(channels)
            % get an image of this channel masked to remove all the cells
            % (all pixels containing a cell are set to 0)
            backmask=imcomplement(wholemaskBWbin).*imgs{ii};
            
            % get the mean fluorescence of pixels that do not contain a
            % cell, and subtract that from the MFI
            backval = mean(imgs{ii}(backmask>0));
            frame_MFIs{ii} = frame_MFIs{ii} - backval;
            MFIs{ii}{t} = frame_MFIs{ii};
        end

%         backval=mean(img(backmask>0));
%         rfpbackval=mean(rfpimg(backmask>0));
%         
%         %Subtract background fluorescence values
%         gfpmfi=gfpmfi-backval;
%         rfpmfi=rfpmfi-rfpbackval;
        
%         gfpmficell{t}=gfpmfi;
%         rfpmficell{t}=rfpmfi;

        % record lineage data
        celllists{t}=double([cellList.cellId{tframes(t)}]);
        ancestorscell{t}=ancestorlist;
        timelabels{t}=tframes(t)*ones(size(frame_MFIs{ii}));
    end
    
    % Define nodes for lineage tree
    nodes = zeros(1,numcells);
    A = find(cellfun('isempty',ancestors) == 0);
    for i=1:length(A)
        nodes(A(i)) = max(ancestors{A(i)});
    end
    
    
    % collect variables for output file as column vectors
    allts=cell2mat(timelabels)';
    allcellids=cell2mat(celllists)';
    alllastancestors=cell2mat(ancestorscell)';
    
    [count, cellnums]=hist(nodes,unique(nodes));
    
    allIDlist(cellfun('isempty', allIDlist))=[];
    cellIDlist=unique(cell2mat(allIDlist));
    
    leafIDs=setdiff(cellIDlist,cellnums);
    leafnum=length(leafIDs);
    leafstatus=ismember(allcellids,leafIDs);
    
    % generate `table` of all parameters, including MFIs
    output=table(allts,allcellids,alllastancestors,leafstatus,'VariableNames',...
        {'tframe','Cell','Last_Ancestor','LeafStatus'});
    for ii = 1:length(channels)
        channel = channels{ii};
        output.(['MFI_' channel]) = cell2mat(MFIs{ii})';
    end

    outputs{fnum,2} = output;
    
    % write a CSV file to `outfilename` with full lineage information
    outfilename = fullfile(output_directory, sprintf(OUT_FILE_PATTERN, img_base_name));
    fprintf('- Writing CSV data to "%s"...\n',outfilename)
    writetable(output,outfilename);

    % optionally, output a `.dat` file for each frame, containing only the 
    % MFI data
    if DAT_FILE_PATTERN
        for ii = 1:length(channels)
            channel = channels{ii};
            channel_MFIs = MFIs{ii};
            for t = 1:length(tframes)
                if iscell(channel_MFIs{t}); mfi = cell2mat(channel_MFIs{t})'; 
                else; mfi = channel_MFIs{t}';
                end
                
                % output a .dat file containing intensity values for cells in this channel in this frame
                datfilename = fullfile(output_directory,sprintf(DAT_FILE_PATTERN, img_base_name, t, channel)); 
                fprintf('-- Writing DAT file to "%s"...\n',datfilename )
                %outfilename=fullfile(output_directory,['gfpmfi_' imgfilestr '.dat']);
                save(datfilename ,'mfi','-ascii')
            end
        end
    end
end

end