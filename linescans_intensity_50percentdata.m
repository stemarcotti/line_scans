%% INPUT %%
% get the parent directory
warning off
uiwait(msgbox('Load parent folder'));
d = uigetdir('');
matlab_folder = cd;
cd(d)

% create list of files with this name
listing_cell = dir('**/cb1_m.tif');
listing_nocb = dir('**/no_cb1_m.tif');

cd(matlab_folder)

%% calculate intensity linescans %%
% for each file
% for file_list = 1:length(listing_cell)

file_list = 1;

% open file
file_cb = listing_cell(file_list).name;
directory_cb = listing_cell(file_list).folder;
file_nocb = listing_nocb(file_list).name;
directory_nocb = listing_nocb(file_list).folder;

nt = length(imfinfo(strcat([directory_cb,'/', file_cb])));

% for each frame
%     for frame = 1:nt-1

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

figure
imshow(im, [])
hold on;
for ii = 1:size(x_line,1)
    
    % line extrems
    x1 = cb_centroid(1);
    y1 = cb_centroid(2);
    x2 = x_line(ii);
    y2 = y_line(ii);
    
    % distance between points [px]
    n_points = ceil(sqrt((x2 - x1).^2 + (y2 - y1).^2)) + 1;
    
    % determine x and y locations along the line
    x_values = round(linspace(x1, x2, n_points));
    y_values = round(linspace(y1, y2, n_points));
    
    % apply line mask
    linescan = zeros(length(x_values),1);
    for iii = 1:length(x_values)
        % linescan on original image
        % linescan(iii,1) = im(y_values(iii),x_values(iii));
        % linescan on im_no_cb
        linescan(iii,1) = im_no_cb(y_values(iii),x_values(iii));
    end
    linescan(linescan == 0) = NaN;
    % find index of first not NaN (out of cell body)
    idx = find(~isnan(linescan), 1);
    linescan = linescan(idx:end);
    intensity_linescan(ii).int = linescan;
    plot(x_values, y_values, 'w.', 'markersize', 0.01)
end

% average linescan for this frame
% find length longest track
len = zeros(length(x_line),1);
for ii = 1:length(x_line)
    len(ii,1) = length(intensity_linescan(ii).int);
end
max_len = max(len);

% initialise matrix
intensity_linescan_matrix = zeros(max_len, length(x_line));

% pull data from structure to matrix
for kkk = 1:length(x_line)
    temp = intensity_linescan(kkk).int;
    temp(end+1:max_len, 1) = NaN;
    intensity_linescan_matrix(:,kkk) = temp;
end

%%
howmany_notnan = sum(~isnan(intensity_linescan_matrix),2);
halfdata_available = round(0.5 * size(intensity_linescan_matrix,2));

idx_halfdata_available = find(howmany_notnan < halfdata_available);
idx_halfdata_available = idx_halfdata_available(1);

intensity_linescan_matrix_half = intensity_linescan_matrix(1:idx_halfdata_available,:);

av_int = nanmean(intensity_linescan_matrix_half, 2);

y_min = min(av_int);
y_max = max(av_int);
y_norm = (av_int-y_min) ./ (y_max - y_min);

save(fullfile(d, 'average_normalised.mat'), 'y_norm');
clear