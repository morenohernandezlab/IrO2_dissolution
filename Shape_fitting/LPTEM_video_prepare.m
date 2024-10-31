%% Liquid Phase Electron Microscopy Video Analyzer
% Reads in-situ video files with the .dm4 format
% Stores relevant information in the video_info structure
% Written by Ivan A. Moreno-Hernandez, starting on 02/24/2022


clear all
clc
% select folder that contains several video folders
mainfolder = uigetdir('','*Select Folder with video sub-folders*'); %Read directory

%%
clearvars -except mainfolder
%% Create video_info structure with information about videos
try
    parpool(16); %Creates parallel pool
end

clc
tic
video_info = prepare_directory(mainfolder);

%%
video_info = last_frame_info(video_info);
video_info = legacy_fix(video_info);
video_info = find_true_fps(video_info);

%% Save in particular directory with a specific name
[filename,out_path] = uiputfile('.mat');
fullpath = convertCharsToStrings(append(out_path,filename));

save(fullpath, 'video_info')

%% Save temp video info
clear all

[filename,in_path] = uiputfile('.mat');
fullpath = convertCharsToStrings(append(in_path,filename));

load(fullpath)

%load Analysis_all_2023_04_20.mat video_info


%% Finding starting frame automatically
video_info = find_start(video_info,1:size(video_info,2));
video_info = lc_time(video_info,1:size(video_info,2));
time_run = toc

 %% Use this section to manually fix scale bars or start frames as needed

%video_info(1).scale_nm_pix = 0.1054060; %For 2K videos taken at 200kX
% Use this part to manually define start times
%video_info(1).start_frame = 1;
%% Run this section if you changed the start frame manually
for i = 1:size(video_info,2)
clear temp_time
clear index
index = video_info(i).start_frame
temp_time = video_info(i).time(index)

video_info(i).time = [video_info(i).time - temp_time]; %-[video_info(video_info(i).start_frame).time]

end

%% Use this section to view videos
clf
v = 1
sp = 120; %Speed factor, greater than 1
video_info = playback_video(video_info,v,sp)

%% Save in particular directory with a specific name
[filename,out_path] = uiputfile('.mat');
fullpath = convertCharsToStrings(append(out_path,filename))

save(fullpath, 'video_info')

%% Functions for code

%This function will find the true frames per second using the folder
%structure
function out = playback_video(x,v,sp)
out = x; %Make sure to keep data
f = figure
f.Position = [1400 100 800 800];
frames_index = round(linspace(1,size(x(v).file_dir,1),round(size(x(v).file_dir,1)/sp)));

for i = frames_index
I = read_frame(x(v).file_dir(i));
out(v).counts(i,:) = mean(I,"all");
I = I./max(I,[],'all');
I = imadjust(I);
imshow(I);
%f.Position = [100 100 540 400];
%pause(1)
title(['Frame: ', num2str(i)])
%i
drawnow
end
end


function x = find_true_fps(x)
for i = 1:size(x,2)
file_list =[];
file_list = vertcat(x(i).frame_list(:).folder);
count = 0;
for j = 1:size(file_list,2)
    if file_list(j,:) == file_list(1,:)
    count = count +1
    end
end
x(i).fps = count;
end
end


% GUI to select line of interest from the video
% Inputs: (video_info, video, number of frames to analyze)\
%[video_info, video #]
function x = prop_finder(x, v,net)
for i = v
    frames = x(i).start_frame:1:x(i).num_frames;
    clear props;
    parfor j = frames %Change to parfor if possible
        Ibw = lc_binarize_nn(x,i,j, net); % Change this line depending on segmentation algorithm
        props(j) = lc_props(Ibw);
        j
    end
    x(i).props = props;
end
end


function out = lc_binarize_nn(x,v,f,net) % Uses neural net
I = read_frame_crop(x(v).file_dir(f), x(v).crop_coord(f,:));
I_size = size(I); %Save original size of frame
Inew = imresize(imadjust(mat2gray(I)),[512,512]); %Prepare new frame for analysis

%Prepare 3-channel image for neural net
    Ich1 = uint8(Inew(:,:,1)*(2^8-1));
    Ibw = imbinarize(imcomplement(Ich1));
    Ich2 = uint8(Ibw*255);
    Idist = bwdist(imcomplement(Ibw));
    Idist = uint8(Idist/max(Idist,[],'all')*255);
    Ich3 = Idist;
    Irgb = uint8(zeros(512,512,3));
    Irgb(:,:,1) = Ich1;
    Irgb(:,:,2) = Ich2;
    Irgb(:,:,3) = Ich3;
% Analyze image with neural net
    C = semanticseg(Irgb, net);
    Ibw = C == 'nanocrystal';
    %Post neural net processing
    Ibw = bwareaopen(Ibw,500);
    %Slight erosion
    SE = strel("square",2);
    Ibw = imerode(Ibw,SE);
    % Fill border and then clear
    Ibw(1:11,:) = 1;
    Ibw(end-10:end,:) = 1;
    Ibw(:,1:11) = 1;
    Ibw(:,end-10:end) = 1;
    %Ibw = imclearborder(Ibw); %Removed for nanowires
    % Smooth image
    windowSize = 11;
    kernel = ones(windowSize) / windowSize ^ 2;
    blurryImage = conv2(single(Ibw), kernel, 'same');
    Ibw_final = blurryImage > 0.5; % Rethreshold
    Ibw_final = imfill(Ibw_final,'holes');
    out = imresize(Ibw_final,I_size); %Resize to original size
end


function out = lc_binarize1(x,v,f) % Not currently working
I = read_frame_crop(x(v).file_dir(f), x(v).crop_coord(f,:));
In = imadjust(mat2gray(I));
%I = imadjust(I);
Ibw = imbinarize(imcomplement(In));
Ibw = imclearborder(Ibw);
Ibw = bwareaopen(Ibw,500);
%Slight erosion

%Blurr sides
windowSize = 41;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(Ibw), kernel, 'same');
Ibw = blurryImage > 0.4;
Ibw = bwareaopen(Ibw,500);
Ibw = imfill(Ibw,'holes');
out = Ibw;
end


function out = crop_rotate(x,v,f,ang_offset, crop_dim);
I = imadjust(mat2gray(read_frame(x(v).file_dir(f))));
offset = ang_offset;
crop_dim_pix = crop_dim./x(v).scale_nm_pix;
crop_dim_pix(1) = - crop_dim_pix(1);
crop_dim_pix(2) = - crop_dim_pix(2);
crop_dim_pix = crop_dim_pix+[double(x(v).lines_center(f,1))-crop_dim_pix(3)*.5,double(x(v).lines_center(f,2))-crop_dim_pix(4)*.5,0,0];
crop_dim_pix = uint32(crop_dim_pix);
I_r = rotateAround(I,double(x(v).lines_center(f,2)),double(x(v).lines_center(f,1)),x(v).lines_angle(f,1)+offset,'bicubic');
out = imcrop(I_r,crop_dim_pix);
end


function x = loi_gui(x,i,n)
for j = i
framestart = x(j).start_frame;
frameend = x(j).num_frames;
frames = round(linspace(framestart,frameend,n));
for k = 1:n
I = read_frame(x(j).file_dir(frames(k)));
set(gcf,'name',x(j).name(frames(k)),'numbertitle','off')
imshow(imadjust(mat2gray(I)))
loi = drawline('InteractionsAllowed','all')
line = customWait(loi)
line_info(k,1:4) = [loi.Position(1,:),loi.Position(2,:)]; %[x1 y1 x2 y2]?
line_center(k,1:2) = [line_info(k,1)+line_info(k,3),line_info(k,2)+line_info(k,4)]/2;
line_angle(k,1) = atan2(line_info(k,4)-line_info(k,2),line_info(k,3)-line_info(k,1))*(180/pi);
%%[x(k),y(k)] = ginput(1)
end

line_info_interp = interp1(frames,line_info,framestart:frameend);
line_center_interp = interp1(frames,line_center,framestart:frameend);
line_angle_interp = interp1(frames,line_angle,framestart:frameend)';

% Fill in for all frames prior to the start frame 
line_info_interp_all = uint32([repelem(line_info_interp(1,:),framestart-1,1); line_info_interp]);
line_info_interp_all = single(line_info_interp_all);

line_center_interp_all = uint32([repelem(line_center_interp(1,:),framestart-1,1); line_center_interp]);

line_angle_interp_all = ([repelem(line_angle_interp(1,:),framestart-1,1); line_angle_interp]);
x(j).lines = line_info_interp_all;
x(j).lines_angle = line_angle_interp_all;
x(j).lines_center = line_center_interp_all;
end
close
end


function x = prop_trajectory(x,v)
for k = v
props = x(k).props;
time = [x(k).time];
for i = 1:size(props,2)
props(i).time = time(i);
end
trajectory = [];
count = 1;
for j = 1:size(props,2)
    
    if props(j).valid == 1
        trajectory(count,1) = props(j).time;
        trajectory(count,2) = props(j).area.*(x(k).scale_nm_pix).^2; % Convert area from pixels ^2 to nm^2
        count = count +1;
    end
end
x(k).trajectory = trajectory;
end
end


function out = prep_rgb(I)
I = mat2gray(I);
I = imadjust(I);
I_rgb(:,:,1) = I;
I_rgb(:,:,2) = I;
I_rgb(:,:,3) = I;
out = I_rgb;
end


function prop_view(x,v,net)
fig = uifigure('Position',[1800-800 100 800 900]);
%create first image
im = uiimage(fig);
I = read_frame_crop(x(v).file_dir(x(v).start_frame), x(v).crop_coord(x(v).start_frame,:));
im.ImageSource = prep_rgb(I);
im.Position = [50 350+150 350 350];
%create second image
im2 = uiimage(fig);
I_bw = lc_binarize_nn(x,v,x(v).start_frame,net);
I_bw_rgb = prep_rgb(I_bw);
I_bw_boundary = [x(v).props(x(v).start_frame).boundary];
I_bw_boundary(:,3) = 1.5; %Dimension of circle
I_bw_rgb = insertShape(I_bw_rgb,'filled-circle',I_bw_boundary);
im2.ImageSource = I_bw_rgb;
im2.Position = [400 350+150 350 350];
% Create text area
im_scale = x(v).crop_coord(x(v).start_frame,3)*x(v).scale_nm_pix;
txa = uitextarea(fig,'Value',{['Image dimensions on each side: ',num2str(im_scale),' nm'],['Displayed data: Frame ', num2str(x(v).start_frame),', Time ', num2str(x(v).time(x(v).start_frame)), ' s']});
txa.Position = [50 860 700 40];
%create plot area
ax = uiaxes(fig,'Position',[0 70 790 400]);
x_plot = [x(v).trajectory(:,1)];
y_plot = [x(v).trajectory(:,2)];
%ax.Box = 'on';
plot(ax,x_plot,y_plot,'Linewidth',3);
ax.YLim = [0,1.3*max(y_plot,[],'all')];
ax.XLim = [0,max(x_plot,[],'all')];
time_from_frame = interp1([x(v).start_frame:x(v).num_frames], ([x(v).start_frame:x(v).num_frames]-x(v).start_frame)./x(v).fps, x(v).start_frame);
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = 'Projected Area (nm^{2})';
ax.YLabel.FontSize = 16;
ax.XLabel.FontSize = 16;

hold(ax,'on')
plot(ax,[time_from_frame,time_from_frame],[0, 1.3*max(y_plot,[],'all')]);
hold(ax,'off')

% Initialize sliders
sld1 = uislider(fig,'Position',[68 40 720 3], 'ValueChangedFcn',@(sld1,event) update_sld(sld1,x,v,im,im2,txa,ax,net));
sld1.Limits = [x(v).start_frame, x(v).num_frames];

end


% Create ValueChangedFcn callback
function update_sld(sld1,x,v,im,im2,txa,ax,net)
value1 = round(sld1.Value);
sld1.Value = value1;
%first Image
I = read_frame_crop(x(v).file_dir(value1), x(v).crop_coord(value1,:));
im.ImageSource = prep_rgb(I);
%second image
I_bw = lc_binarize_nn(x,v,value1,net);
I_bw_rgb = prep_rgb(I_bw);
I_bw_boundary = [x(v).props(value1).boundary];
I_bw_boundary(:,3) = 1.5; %Dimension of circle
I_bw_rgb = insertShape(I_bw_rgb,'filled-circle',I_bw_boundary);
im2.ImageSource = I_bw_rgb;
im_scale = x(v).crop_coord(value1,3)*x(v).scale_nm_pix;
time_from_frame = interp1([x(v).start_frame:x(v).num_frames], ([x(v).start_frame:x(v).num_frames]-x(v).start_frame)./x(v).fps, value1);

%text
txa.Value = {['Image dimensions on each side: ',num2str(im_scale),' nm'],['Displayed data: Frame ', num2str(value1),', Time ', num2str(time_from_frame), ' s']};

%update plot area
x_plot = [x(v).trajectory(:,1)];
y_plot = [x(v).trajectory(:,2)];
plot(ax,x_plot,y_plot,'Linewidth',3);
ax.YLim = [0,1.3*max(y_plot,[],'all')];
hold(ax,'on')
plot(ax,[time_from_frame,time_from_frame],[0, 1.3*max(y_plot,[],'all')]);
hold(ax,'off')

end


%[video_info, video #, frame #]
function x = lc_time(x,v)
for i = v
    
    x(i).time = (single([1:x(i).num_frames]') - single(x(i).start_frame))./single(x(i).fps);
end

end


function out = lc_props(I)
try
    % Determine properties of area
    props = regionprops(I,{'Area','Perimeter','BoundingBox','Image','PixelIdxList','Orientation','Solidity','Perimeter','Centroid'});
    % Find index for particle closest to center of image
    coords_centroid = vertcat(props.Centroid)-size(I)/2;
    dist_centroid = (coords_centroid(:,1).^2 + coords_centroid(:,2).^2).^(1/2);
    [~,indx] = min(dist_centroid);
    
    %Determine outline of particle
    B = bwboundaries(I,'noholes'); %Find particles
    boundary = cell2mat(B(indx));  %output boundary
    coord_boundary = [boundary(:,2),boundary(:,1)].*[1,-1];
    coord_boundary = coord_boundary-mean(coord_boundary); %Coordinates in [X,Y], units in pixels
    %Save results to output
    out.area = props(indx).Area; %output area
    out.perimeter= props(indx).Perimeter; %output perimeter
    out.boundary = [boundary(:,2),boundary(:,1)]; %Flip coordinates to match image
    out.coord_boundary = coord_boundary;
    out.solidity = props(indx).Solidity; %output Solidity
    out.orientation = props(indx).Orientation; %output Orientation
    out.centroid = props(indx).Centroid; %output Centroid
    out.valid = 1; %report if property finder worked
catch
    out.area = 0; %output area
    out.perimeter= 0; %output perimeter
    out.boundary = [0,0];  %output boundary
    out.coord_boundary = [0,0];  %output boundary
    out.solidity = 0; %output Solidity
    out.orientation = 0; %output Solidity
    out.centroid = [0,0]; %output Centroid
    out.valid = 0; %report if property finder worked

end
end


% Find all folders within a folder
function out = folders_in(x)
d = dir(x);
% remove all files (isdir property is 0)
dfolders = d([d(:).isdir]==1);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
out = strings(length(dfolders),1);
for floop = 1:length(dfolders)
out(floop,1) = strcat(dfolders(floop).folder,'\',dfolders(floop).name);
end
end

% find all .dm3 frames in all folders and subfolders
function out = find_frames(x)
filelist = dir(fullfile(x, '**\*.dm4'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list
out = filelist;
end

% This part of the code is used to prepare the locations of all of the
% frames from the videos
function out = prepare_directory(x)
folders = folders_in(x);
for i = 1:size(folders,1)
    out(i).video_folder = folders(i);
    out(i).frame_list = find_frames(folders(i));
    file_dir = "";
    name = "";
    for j = 1:size(out(i).frame_list,1)
        file_dir = [file_dir; strcat(out(i).frame_list(j).folder,'\',out(i).frame_list(j).name)];
        name = [name; out(i).frame_list(j).name];
    end
    file_dir = [file_dir(2:end)]; % Remove initial empty string
    name = [name(2:end)]; % Remove initial empty string
    out(i).file_dir = file_dir;
    out(i).name = name;
    out(i).num_frames = size(file_dir,1);
    
end
end


% Determines information from last video frame
function x = last_frame_info(x)

for i = 1:size(x,2)
    i
    Imagestruc = dmread(x(i).file_dir(end)); %reads dm4 file
    x(i).x_pix = Imagestruc.ImageList.Unnamed0.ImageData.Dimensions.Unnamed0.Value; %extracts x-dimension of image
    x(i).y_pix =Imagestruc.ImageList.Unnamed0.ImageData.Dimensions.Unnamed1.Value; %extracts y-dimension of image
    x(i).scale_nm_pix =Imagestruc.ImageList.Unnamed0.ImageData.Calibrations.Dimension.Unnamed0.Scale.Value; %extracts scale of image (nm/pixel)
    if isequal(Imagestruc.ImageList.Unnamed0.ImageData.Calibrations.Dimension.Unnamed0.Scale.Value,1)
       x(i).scale_nm_pix = 0.14060315; %Change as needed 
    end
    x(i).intscale = 1/53.6 %Imagestruc.ImageList.Unnamed0.ImageData.Calibrations.Brightness.Scale.Value ; %extracts intensity scale for image (e-/value)
    x(i).fps = 1/Imagestruc.ImageList.Unnamed0.ImageTags.DataBar.ExposureTime0x28s0x29.Value; % FPS = 1/exposure
    x(i).final_dose =  x(i).fps*mean(Imagestruc.ImageList.Unnamed0.ImageData.Data.Value)*x(i).intscale/((100)*x(i).scale_nm_pix.^2); % image dose (e- A^-2 s^-1)
end
end


%reads dm4 file and outputs array in e-/(pixel)
function out = read_frame(x) %
    Imagestruc = dmread(x); 
    x_pix = Imagestruc.ImageList.Unnamed0.ImageData.Dimensions.Unnamed0.Value; %extracts x-dimension of image
    y_pix =Imagestruc.ImageList.Unnamed0.ImageData.Dimensions.Unnamed1.Value; %extracts y-dimension of image
    intscale =Imagestruc.ImageList.Unnamed0.ImageData.Calibrations.Brightness.Scale.Value ; %extracts intensity scale for image (e-/value)
    out = reshape(Imagestruc.ImageList.Unnamed0.ImageData.Data.Value,[x_pix,y_pix])'; %Normalize later
end


function out = read_frame_crop(x,crop_coord)
    I = read_frame(x);
    out = imcrop(I,[crop_coord(1), crop_coord(2), 0,0]+crop_coord(3).*[-.5, -.5, 1, 1]);
end


%Measures drift between two images using cross-correlation algorithm
%From: https://medium.com/@windless99/drift-correction-in-matlab-f9cc02860c09
function out = drift_coord(im1,im2)
[sy,sx]=size(im1);
fftim1=fft2(im1);
fftim2=fft2(im2);
cc=fftshift(ifft2(fftim1.*conj(fftim2)));
[shiftY,shiftX]=find(cc==max(cc(:)));
shiftY=shiftY-fix(sy/2)-1;
shiftX=shiftX-fix(sx/2)-1;
out = -[shiftX,shiftY]; % Correct for sign of drift
end


function out = flatten(I)
c_size = 32; %compression size, pixels
I_s = imresize(I,[c_size,c_size]); %Make small copy of image
[X,Y] = meshgrid(1:c_size,1:c_size); %Define coordinates
sf = fit([X(:), Y(:)],I_s(:),'poly23'); %Fit to polynomial
fit_points = feval(sf,[X(:), Y(:)]); %Get fit points
I_bg_s = reshape(fit_points,[c_size,c_size]); %Get background fit
I_bg_l = imresize(I_bg_s,size(I)); %Resize to normal size
out = I-I_bg_l; %Remove background
end


function view_frame(I)
imshow(imadjust(mat2gray(I)));
end


function x  = find_start(x,k)
parfor i = k
%Determine initial frame by comparing electrons detected to previous frames
count = size(x(i).file_dir,1);
trajectory_dose_frac =1;
%Fast backward scan
I_final_sum = sum(read_frame(x(i).file_dir(end)),'all');
while trajectory_dose_frac > 0.5; 
    Im = read_frame(x(i).file_dir(count)); % Read frame from count
    trajectory_dose_frac = sum(Im,'all')./I_final_sum;
 count = count-10%
 if count < 1
     count = 1;
     trajectory_dose_frac =0
 end
end


%Frame-by-frame forwards scan
while trajectory_dose_frac < 0.7;
if count < 1
     count = 1;
end
    Im = read_frame(x(i).file_dir(count)); % Read frame from count
    trajectory_dose_frac = sum(Im,'all')./I_final_sum;
 count = count+1 %
if count < 1
     count = 1;
end
end
x(i).start_frame = count;    
end
end


% GUI to select polygon region of interest from the video
% Inputs: (video_info, video, number of frames to analyze)
function x = roi_polygon(x,i)
    I = read_frame(x(i).file_dir(x(i).start_frame)); %Read frame for analysis
    I = I./max(I,[],'all');
    I = imadjust(I);
    imshow(I)
    roi = drawpolygon %Props user to draw rough polygon
    n_coord_initial = round(roi.Position);
    x(i).polygon_coord_initial = n_coord_initial;
end


% GUI to select region of interest from the video
% Inputs: (video_info, video, number of frames to analyze)
function x = roi_gui(x,i,n)
for j = i
framestart = x(j).start_frame;
frameend = x(j).num_frames;
frames = round(linspace(framestart,frameend,n));
for k = 1:n
I = read_frame(x(j).file_dir(frames(k)));
set(gcf,'name',x(j).name(frames(k)),'numbertitle','off')
imshow(imadjust(mat2gray(I)));
roi = drawrectangle('AspectRatio',1,'FixedAspectRatio',true,'InteractionsAllowed','all')
rectangle = customWait(roi);
% [xmin ymin width height]
center_x(k) = rectangle(1) + rectangle(3)/2;
center_y(k) = rectangle(2) + rectangle(4)/2;
width_frame(k) = ((rectangle(3)+rectangle(4))/2);
%[x(k),y(k)] = ginput(1)
end
width_frame_int = interp1(frames,width_frame,framestart:frameend);
x_int = interp1(frames,center_x,framestart:frameend);
y_int = interp1(frames,center_y,framestart:frameend);

% Place crop coordinates in convenient array
crop_coord_initial = [x_int(:), y_int(:), width_frame_int(:)]; 
% Fill in for all frames prior to the start frame 
crop_coord = uint32([repelem(crop_coord_initial(1,:),framestart-1,1); crop_coord_initial]);
crop_coord = single(crop_coord);
x(j).crop_coord = crop_coord;
end
close
end


% Wait function
function pos = customWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);

% Return the current position
pos = hROI.Position;

end


function clickCallback(~,evt)
if strcmp(evt.SelectionType,'double')
    uiresume;
end
end


function x = legacy_fix(x)
for i = 1:size(x,2)
test = x(i).frame_list(1).folder;
if isequal(test(end-8:end),'Second_00')

else
x(i).frame_list(1) =[]; %Remove the first file
x(i).file_dir(1) =[]; %Remove the first file
x(i).name(1) =[]; %Remove the first file
x(i).num_frames = x(i).num_frames - 1;
end
end
end