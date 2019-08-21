%% INPUT %%
% get the parent directory
warning off
uiwait(msgbox('Load folder containing image'));
d = uigetdir('');
matlab_folder = cd;
cd(d)

% ask the user for an ouput stamp
prompt = {'Linescan thickness [px] - assign 0 if no width required'};
title_prompt = 'Parameters';
dims = [1 35]; % set input box size
definput = {'10'};
user_answer = inputdlg(prompt,title_prompt,dims,definput); % get user answer
thickness = str2double(user_answer{1,1});

% create list of files with this name
listing_cell = dir('**/cb1_m.tif');
listing_nocb = dir('**/no_cb1_m.tif');

cd(matlab_folder)

%% set up %%
% for each file
% for file_list = 1:length(listing_cell)

file_list = 1;

% open file
file_cb = listing_cell(file_list).name;
directory_cb = listing_cell(file_list).folder;
file_nocb = listing_nocb(file_list).name;
directory_nocb = listing_nocb(file_list).folder;

nt = length(imfinfo(strcat([directory_cb,'/', file_cb])));

frame = 1;

% read image at current frame
im = double(imread(fullfile([directory_cb '/' file_cb]), frame));          % original movie
im_no_cb = double(imread(fullfile([directory_nocb '/' file_nocb]), frame));  % movie without cell body

% find difference between the two images: cell body and centroid
cb = logical(im) - logical(im_no_cb);
mask_cb = imbinarize(cb);
mask_cb = imfill(mask_cb, 'holes');

stats_cb = regionprops(mask_cb, 'Area', 'Centroid');
stats_area = [stats_cb.Area];
idx_max_area = find(stats_area == max(stats_area));

stats_cb = stats_cb(idx_max_area);
cb_centroid = cat(1, stats_cb.Centroid);

% find cell edge
edge_line = edge(logical(im), 'Canny'); % get cell edge line
[y, x] = find(edge_line);

% draw linescans from centroid to points on the cell edge
% measure intensity
frequency = 50; % every x pts on the edge
x_line = x(1:frequency:end);
y_line = y(1:frequency:end);

%% linescans %%
intensity_linescan = struct([]);
intensity_linescan_nocb = struct([]);
for ii = 1:size(x_line,1)
    
    % line extremes
    x1 = cb_centroid(1);
    y1 = cb_centroid(2);
    x2 = x_line(ii);
    y2 = y_line(ii);
    
    % distance between points [px]
    n_points = ceil(sqrt((x2 - x1).^2 + (y2 - y1).^2)) + 1;
    
    % determine x and y locations along the line
    x_values = round(linspace(x1, x2, n_points));
    y_values = round(linspace(y1, y2, n_points));
    
    % average linescans vertically across defined thickness
    linescan_thick = zeros(size(x_values,2),1);
    linescan_thick_nocb = zeros(size(x_values,2),1);
    for jj = 1:size(x_values,2) % for each point on the line
        
        % define perpendicular line across thickness
        x1v = x_values(jj);
        x2v = x_values(jj);
        y1v = y_values(jj) - thickness/2;
        y2v = y_values(jj) + thickness/2;
        
        n_points_v = ceil(sqrt((x2v - x1v).^2 + (y2v - y1v).^2)) + 1;
        
        x_temp = round(linspace(x1v, x2v, n_points_v));
        y_temp = round(linspace(y1v, y2v, n_points_v));
        
        % scan perpendicular line to get intensity values
        linescan_v = zeros(length(x_temp),1);
        linescan_v_nocb = zeros(length(x_temp),1);
        for jj_temp = 1:length(x_temp)
            linescan_v(jj_temp,1) = im(y_temp(jj_temp),x_temp(jj_temp));            % with cell body
            linescan_v_nocb(jj_temp,1) = im_no_cb(y_temp(jj_temp),x_temp(jj_temp)); % without cell body
        end
        
        % average intensity across thickness
        linescan_v(linescan_v == 0) = NaN;
        linescan_thick(jj,1) = nanmean(linescan_v);
        
        linescan_v_nocb(linescan_v_nocb == 0) = NaN;
        linescan_thick_nocb(jj,1) = nanmean(linescan_v_nocb); 
    end
    
    intensity_linescan(ii).int = linescan_thick;
    
    % find index of first not NaN (out of cell body)
    idx = find(~isnan(linescan_thick_nocb), 1);
    linescan_thick_nocb = linescan_thick_nocb(idx:end);
    intensity_linescan_nocb(ii).int = linescan_thick_nocb;
end

%% unpack structure and average linescans %%

% average linescan for this frame
% find length longest track
len = zeros(length(x_line),1);
len_nocb = zeros(length(x_line),1);
for ii = 1:length(x_line)
    len(ii,1) = length(intensity_linescan(ii).int);
    len_nocb(ii,1) = length(intensity_linescan_nocb(ii).int);
end
max_len = max(len);
max_len_nocb = max(len_nocb);

% initialise matrix
intensity_linescan_matrix = zeros(max_len, length(x_line));
intensity_linescan_nocb_matrix = zeros(max_len_nocb, length(x_line));

% pull data from structure to matrix
for kkk = 1:length(x_line)
    
    temp = intensity_linescan(kkk).int;
    temp(end+1:max_len, 1) = NaN;
    intensity_linescan_matrix(:,kkk) = temp;
    
    temp_nocb = intensity_linescan_nocb(kkk).int;
    temp_nocb(end+1:max_len_nocb, 1) = NaN;
    intensity_linescan_nocb_matrix(:,kkk) = temp_nocb;
end

%% remove last noisy bit (keep until 50% of lines are available) %%

howmany_notnan = sum(~isnan(intensity_linescan_matrix),2);
halfdata_available = round(0.5 * size(intensity_linescan_matrix,2));

howmany_notnan_nocb = sum(~isnan(intensity_linescan_nocb_matrix),2);
halfdata_available_nocb = round(0.5 * size(intensity_linescan_nocb_matrix,2));

idx_halfdata_available = find(howmany_notnan < halfdata_available);
idx_halfdata_available = idx_halfdata_available(1);

idx_halfdata_available_nocb = find(howmany_notnan_nocb < halfdata_available_nocb);
idx_halfdata_available_nocb = idx_halfdata_available_nocb(1);

intensity_linescan_matrix_half = intensity_linescan_matrix(1:idx_halfdata_available,:);
intensity_linescan_nocb_matrix_half = intensity_linescan_nocb_matrix(1:idx_halfdata_available_nocb,:);

av_int = nanmean(intensity_linescan_matrix_half, 2);
av_int_nocb = nanmean(intensity_linescan_nocb_matrix_half, 2);

%% save %%

save(fullfile(d, 'av_int_withcb.mat'), 'av_int');
save(fullfile(d, 'av_int_withoutcb.mat'), 'av_int_nocb');
clear
