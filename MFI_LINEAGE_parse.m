
% Daniel Lee, Christina Lin, and Casey Grun
% MIT License

% This script uses cell mesh and (optionally) lineage data produced by
% Oufti, along with cell images, to produce a CSV file with cell lineage
% and MFI data. The results are saved as a CSV file, where each row is 
% a timepoint for a particular cell, and columns the following information: 
%   tframe : frame number
%   Cell : cell number from Oufti
%   Last_Ancestor : cell number of most recent ancestor of this cell
%   LeafStatus : 1 if this cell is a leaf (does not have any descendents)
%   MFI_{channel1} : mean fluorescence intensity of the cell in this
%     channel, after subtracting the average of the average intensity of
%     the background (parts of the image with no cells). This value is NOT
%     normalized and thus values depend on the range of intensities in the
%     input image. 
%   MFI_{channel2} : additional columns as necessary for other channels.

clear
close all

% input directory
directory='./Example Datasets/PA14 WT GFP RFP Colony 20/';

% output directory
outdirname='./Example Datasets/PA14 WT GFP RFP Colony 20/lineage';

% which channels will you use? Your image files must be named using these
% channel names (case sensitive). You can have more than two channels. 
channels = { 'GFP', 'RFP' };

% cell mesh and lineage data should be in an Oufti-produced .mat file. 
% how are your `.mat` files named? This should be a regular expression
% containing at least `channel` and `image` tokens/capture groups. See
% documentation for the `regexp` function. 
% this particular pattern expects files like this:
% Ph-some_long_experiment_name.mat
MAT_FILE_PATTERN = '^(?<channel>\w+)-(?<image>.*)\.(?<ext>\w+)';

% phase and channel images should be in separate (per-channel) multi-page
% TIF files. 
% how are your `.tif` image files named? should be a format specifier
% suitable for `sprintf`. `%{channel}s` is the channel, `%{image}s` is the 
% base image name. This particular pattern expects files like this:
% Ph-some_long_experiment_name.tif
TIF_FILE_PATTERN = '%{channel}s-%{image}s.tif';

% how do you want your output files named? should be a format string
% suitable for `sprintf` with on1 placeholders: `%{image}s` is the base 
% image name
OUT_FILE_PATTERN = '%{image}s.csv';

% which frames to export? if this set to the empty array `[]`, all frames
% with 1+ cells will be exported. Note: the Last_Ancestor column in the CSV
% file may refer to cells birthed on frames outside this range. 
frames_to_export = [];

%%
TIF_FILE_PATTERN = strrep(TIF_FILE_PATTERN, '{channel}','1$'); 
TIF_FILE_PATTERN = strrep(TIF_FILE_PATTERN, '{image}','2$');
OUT_FILE_PATTERN = strrep(OUT_FILE_PATTERN, '{image}','1$');

% create output directory
mkdir(outdirname)

% Load Oufti output .mat files from directory
matfileobj=dir(fullfile(directory, '*.mat'));
[matfilenames{1:length(matfileobj)}]=matfileobj(:).name;
numfiles=length(matfileobj);
%%

%loop through colony files
for fnum=1:numfiles
    
    % check that this file matches the expected file name pattern, and 
    % extract the base image name. If this is not a Phase image, skip
    file_name_parts = regexp(matfilenames{fnum}, MAT_FILE_PATTERN, 'names');
    if (isempty(file_name_parts))
        fprintf('Skipping .mat file "%s" which does not fit the expected pattern\n', matfilenames{fnum})
        continue
    elseif (~strcmp(file_name_parts.channel,'Ph')) 
        fprintf('Skipping .mat file "%s" which does not appear to be a phase contrast image\n', matfilenames{fnum})
        continue
    else 
        fprintf('Processing .mat file "%s"...\n', matfilenames{fnum})
    end
    
    % get complete path to file, extract base image name
    matfilename = fullfile(directory, matfilenames{fnum});
    img_base_name = file_name_parts.image; 
    
    % get path to TIF file for each channel
    imgfilenames = cellfun(...
        @(channel) fullfile(directory, ...
            sprintf(TIF_FILE_PATTERN,channel,img_base_name)),...
        channels, 'UniformOutput', false);
    
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
    
    %loop through timepoints
    for t=1:length(tframes)
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
        
        %Define GFP and RFP MFIs
        MFIs = cellfun(@(~) zeros(1,cellListN(tframes(t))), channels, ...
            'UniformOutput', false);
%         gfpmfi=zeros(1,cellListN(tframes(t)));
%         rfpmfi=zeros(1,cellListN(tframes(t)));

        % setup matrices to capture linage data
        ancestorlist=zeros(1,cellListN(tframes(t)));
        leafstatus=zeros(1,cellListN(tframes(t)));
        
        %loop through cells
        for cellnum=1:cellListN(tframes(t))
            meshcoord=meshinfo{cellnum}.mesh;
            if meshcoord==0
%                 gfpmfi(cellnum)=NaN;
%                 rfpmfi(cellnum)=NaN;
            else
                
                % generate a mask identifying this cell, from its mesh
                borders=[flip(meshcoord(:,4)) flip(meshcoord(:,3)); meshcoord(:,2) meshcoord(:,1)];
                infrows=find(ismember(borders,[Inf Inf],'rows')==1);
                if isempty(infrows)==0
                    borders(infrows,:)=[];
                end
                bactmask=poly2mask(...
                    (double(borders(:,2))),...
                    (double(borders(:,1))),...
                    size(imgs{1},1),...
                    size(imgs{1},2));
                
                % for each channel, mask the channel image to only include
                % this cell, then take the mean value (MFI)
                for ii = 1:length(channels)
                    MFIs{ii}(cellnum) = mean(imgs{ii}(bactmask==1));
                end
%                 gfpmfi(cellnum)=mean(img(bactmask==1));
%                 rfpmfi(cellnum)=mean(rfpimg(bactmask==1));
                
                % keep track of which cell occupies this pixel
                wholemask=wholemask+cellnum*bactmask;
                
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
            backval=mean(imgs{ii}(backmask>0));
            MFIs{ii}=MFIs{ii}-backval;
            MFIcell{ii}{t} = MFIs{ii};
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
        timelabels{t}=tframes(t)*ones(size(MFIs{ii}));
    end
    
    %Define nodes for lineage tree
    nodes=zeros(1,numcells);
    A=find(cellfun('isempty',ancestors)==0);
    for i=1:length(A)
        nodes(A(i))=max(ancestors{A(i)});
    end
    
    
    % collect variables for output file as column vectors
    
    allts=cell2mat(timelabels)';
    allcellids=cell2mat(celllists)';
    alllastancestors=cell2mat(ancestorscell)';
    
    [count cellnums]=hist(nodes,unique(nodes));
    
    allIDlist(cellfun('isempty', allIDlist))=[];
    cellIDlist=unique(cell2mat(allIDlist));
    
    leafIDs=setdiff(cellIDlist,cellnums);
    leafnum=length(leafIDs);
    leafstatus=ismember(allcellids,leafIDs);
    
    % generate `table` of all parameters, including MFIs
    output=table(allts,allcellids,alllastancestors,leafstatus,'VariableNames',{'tframe','Cell','Last_Ancestor','LeafStatus'});
    for ii = 1:length(channels)
        channel = channels{ii};
        output.(['MFI_' channel]) = cell2mat(MFIcell{ii})';
    end
    
    % write to `outfilename`
    outfilename=fullfile(outdirname, sprintf(OUT_FILE_PATTERN, img_base_name));
    writetable(output,outfilename);
end