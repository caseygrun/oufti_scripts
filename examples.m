% clean up the workspace and close figures
clear variables;
close all;

%%

% parse MFIs for a dataset with several single-frame images, GFP-only
outputs = parse_MFI_lineage('example_data/PA14_WT_12H','example_data/example_output/PA14_WT_12H_parsed', 'Channels', {'GFP'});

%%

% look at the results
fprintf('Found %d colonies: \n', length(outputs(:,1)))
disp(strjoin(outputs(:,1),'\n'))

% plot a histogram of GFP MFIs for all colonies
plot_MFIs({outputs{:,2}}, 'Channels', {'GFP'});

%%

% plot a single frame from one colony and fit a Gaussian mixture model to calculate the
% percent ON cells
[percent_on, model, ~, ~] = plot_gaussian(outputs{1,2}.MFI_GFP, 1, 25);
fprintf('In colony 1, %f percent T3SS on\n', percent_on*100) 

%%

% (this will take a few minutes to run)

% use some data from \Delta _exsA_ (constitutively OFF) and \Delta _exsD_
% (constitutively ON) cells to calculate a threshold MFI, above which 50% 
% of cells with be ON. This threshold will be calculated for each timepoint
% and saved in `classinfo`
% classinfo = calculate_threshold('example_data/control_data',...
%     {'PA14 exsD Pt-sfGFP s1 col1','PA14 exsD Pt-sfGFP s3 col2','PA14 exsD Pt-sfGFP s4 col3','PA14 exsD Pt-sfGFP s4 col4','PA14 exsD Pt-sfGFP s4 col5'},...
%     {'PA14 exsA Pt-sfGFP MinS_NTA s1 col1','PA14 exsA Pt-sfGFP MinS_NTA s1 col2','PA14 exsA Pt-sfGFP MinS_NTA s3 col3','PA14 exsA Pt-sfGFP MinS_NTA s3 col4','PA14 exsA Pt-sfGFP MinS_NTA s4 col5'},...
%     'example_data/classinfo.csv',...
%     'DatFilePattern','gfpmfi_t(?<frame>\d+)\.dat');

classinfo = calculate_threshold(...
    'example_data/control_data/PA14 exsA Pt-sfGFP MinS_NTA*.csv',...
    'example_data/control_data/PA14 exsD Pt-sfGFP*.csv',...
    'example_data/classinfo.csv');


%%

% classify our colonies from earlier as GFP-ON or OFF, using the exsA and
% exsD control data
[onoff, f_on] = classify_colonies('example_data/example_output/PA14_WT_12H_parsed',classinfo,'Channel','GFP');


%% 

% now parse MFIs for a dataset with a single microcolony with GFP and RFP
% channels
outputs_GFP_RFP = parse_MFI_lineage('example_data/PA14_WT_GFP_RFP_Colony_20',...
    'example_data/example_output/PA14_WT_GFP_RFP_Colony_20_parsed', 'Channels', {'GFP', 'RFP'});

%%
plot_MFIs({outputs_GFP_RFP{:,2}}, 'Channels', {'GFP','RFP'});

%% 
plot_MFIs_scatter({outputs_GFP_RFP{:,2}}, 'GFP', 'RFP', 'tframes', [1, 5, 10, 15, 20, 25]);