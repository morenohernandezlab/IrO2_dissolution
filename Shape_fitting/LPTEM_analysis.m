%% This script will generate dissolution trajectories for IrO2 nanocrystals starting from in situ TEM datasets
% Matteo Fratarcangeli
% Ivan A. Moreno-Hernandez
% August 2024
%
% for information on preparing the video_info file, see https://github.com/morenohernandezlab/IrO2_Dissolution

clear all
clc
[file,path] = uigetfile('*.mat'); % opens file selection dialog box. You may choose data from different folder
video_info = importdata(strcat(path,file)); % data will be imported from the specified path

%% Use this section to view videos
clf
v = 2;
sp = 120; %Speed factor, greater than 1
video_info = playback_video(video_info,v,sp);

%% Alline particle frames
% Define line of interest
v = 2;
video_info = loi_gui(video_info,v,15) %Draw line from top to bottom, third input is number of frames

%% Pull frames from video
%Show first and last frames of video
%Define image parameters
clf
I_level = 1; %Brightness level, 0-1
aspect = 1; %Aspect ratio of image
im_int = 5; %Frames used for image integration (integer, frames +- the frames listed)
end_factor = 1;
time_step = 10; %s
end_frame = round(end_factor*size(video_info(v).frame_list,1));
frames = video_info(v).start_frame:round(time_step*video_info(v).fps):end_frame;

%Pull line infro from_video_info and define aspect ratio
line_info_all = video_info(v).lines;
line_length = max(round(((line_info_all(:,1)-line_info_all(:,3)).^2+(line_info_all(:,2)-line_info_all(:,4)).^2).^(1/2)));
h = round(line_length*1.3); %Height of box in pixels
w = round(h/aspect);  %Width of box in pixels

I_start = crop_box_line2(video_info,v,frames(1),h,w,I_level,im_int);
I_end = crop_box_line2(video_info,v,frames(end),h,w,I_level,im_int);
imshowpair(I_start,I_end,'montage')

%% Create guesses for L_initial and L_final
%% Guess 1
L1 = 50; %thickness
L7 = 220; % length
f = [0.75,0.75,0.05,0.05]; % [top right, top left, bottom righ, bottom left]
angles = [45 90 90]; %Rotate in degrees with z plane being the screen, [center axis rotation, 45, 105]
L_projected = max(abs(sqrt(2)*L1*cos(angles(1)*pi/180)),abs(sqrt(2)*L1*cos((angles(1)+90)*pi/180)));
L1 = L1./L_projected*L1;
L2 = L1;
%L2 = L1./(sqrt(2)*cos(angles(1)*pi/180)) %(110) side 2
L3 = L7*f(1); %(111) tip 1 top
L4 = L7*f(2); %(111) tip 2 top
L5 = L7*f(3); %(111) tip 1 bottom
L6 = L7*f(4);  %(111) tip 2 bottom
L = [L1,L2,L3,L4,L5,L6,L7];
L_initial = L;
tic
clear report_fit

% Determine first two guesses
tic
frame_ends = [1,size(frames,2)];
pixel_scale_A = double(video_info(v).scale_nm_pix*10); %A/pixel

translation = [h/2,w/2]*pixel_scale_A; %position_tip+L5*.33*[sin((90-angle_tip)*pi/180),cos((90-angle_tip)*pi/180)] %Translate in A
[shp3d,shp2d] = nanocrystal_IrO2_shape(L,angles,translation);
Ibw = ibw_shape(I_start,pixel_scale_A,shp2d);
clf
imshowpair(I_start,Ibw)
%set(gcf,'units','pixel','position',[1500 25 800 975]);

%% Guess 2
L1_f = L1*0.4; %(110) side 1, thickness
L7_f = L7*0.5; %(001) tips, length
f_f = f.*[1,1,1,1];
angles = [45+0 90 90]; %Rotate in degrees %[90, 45, 105]

L2_f = L1_f; %(110) side 2
L3_f = L7_f*f_f(1); %(111) tip 1 top
L4_f = L7_f*f_f(2); %(111) tip 2 top
L5_f = L7_f*f_f(3); %(111) tip 1 bottom
L6_f = L7_f*f_f(4);  %(111) tip 2 bottom
L_final = [L1_f,L2_f,L3_f,L4_f,L5_f,L6_f,L7_f];
tic
clear report_fit

% Determine first two guesses
translation = [h/2,w/2]*pixel_scale_A; %position_tip+L5*.33*[sin((90-angle_tip)*pi/180),cos((90-angle_tip)*pi/180)] %Translate in A
[shp3d,shp2d] = nanocrystal_IrO2_shape(L_final,angles,translation);
Ibw = ibw_shape(I_end,pixel_scale_A,shp2d);
clf
imshowpair(I_end,Ibw)
%set(gcf,'units','pixel','position',[1500 100 800 800]);

L_trial = [L_initial; L_final];

L_guess = interp1(frames(frame_ends),L_trial,frames);

%% Save fit guesses (not angles)
%save all the info to video_info
video_info(v).L_guess = L_guess;
video_info(v).h = h;
video_info(v).w = w;
video_info(v).I_level = I_level;
video_info(v).im_int = im_int;
video_info(v).frames = frames;

save(strcat(path,file), 'video_info')

%% Draw outline
for v = v

 for i = 1:1:size(video_info(v).frames,2)
complete = i/size(video_info(v).frames,2)*100
     %if size(video_info(v).Ibw{video_info(v).frames(i)},1)> 2

     %else
     I = crop_box_line2(video_info,v,video_info(v).frames(i),video_info(v).h,video_info(v).w,.5,1); % Get image
     I = imadjust(I);
     imshow(I)
     %set(gcf,'units','pixel','position',[1500 100 800 800]);
    roi = drawpolygon
    Ibw = createMask(roi);
    imshow(Ibw)
    %set(gcf,'units','pixel','position',[1500 100 800 800]);
    video_info(v).Ibw{video_info(v).frames(i)} = Ibw;
    clf
     %end
     i
 end

end

%% save outlines
save(strcat(path,file), 'video_info')

%% Fit frames
tic
for v = 2%4:size(video_info,2)
    v
    %try put back on
frame_list = [];

%make a list of frames that have actual videos
parfor i = 1:size(video_info(v).frames,2) %make par

try
Ibw =  video_info(v).Ibw{video_info(v).frames(i)};
if size(Ibw,2) > 2
frame_list = [frame_list, i];
end
end
end 

% fit the actual frames
    clear fit_saved

parfor i = 1:size(frame_list,2) %make par later
try
i

Ibw =  video_info(v).Ibw{video_info(v).frames(frame_list(i))};

%Calculate picture limits
lim_x = [min(find(sum(Ibw,2) > 0))-10,max(find(sum(Ibw,2) > 0))+10];
lim_y = [min(find(sum(Ibw,1) > 0))-10,max(find(sum(Ibw,1) > 0))+10];

%Limiters
lim_x = max(lim_x,1);
lim_y = max(lim_y,1);
lim_x = min(lim_x,size(Ibw,1));
lim_y = min(lim_y,size(Ibw,2));

Ibw2 = Ibw(lim_x(1):lim_x(2),lim_y(1):lim_y(2));
imshow(Ibw2)
L_guess = video_info(v).L_guess(i,:);
pixel_scale_A = double(video_info(v).scale_nm_pix*10); %A/pixel
translation = [size(Ibw2,1)/2,size(Ibw2,2)/2]*pixel_scale_A;
%translation = [video_info(v).h/2,video_info(v).w/2]*pixel_scale_A;
angles = [45, 90, 90];
fit_scale_down = 2^0;
tic
[L_fit,angle_fit,translation_fit] = fit_IrO2_v3(L_guess,angles,translation,imresize(double(Ibw2),1/fit_scale_down),pixel_scale_A*fit_scale_down)

translation = [video_info(v).h/2,video_info(v).w/2]*pixel_scale_A;
[L_fit,angle_fit,translation_fit] = fit_IrO2_pos(L_fit,angle_fit,translation,imresize(double(Ibw),1/fit_scale_down),pixel_scale_A*fit_scale_down)

run_time = toc
fit_saved(i,:) = [video_info(v).frames(frame_list(i)), video_info(v).time(video_info(v).frames(frame_list(i))), L_fit,angle_fit,translation_fit];
%video_info(v).fit = fit_saved; %fit_saved;

%video_info(v).fit(i,:) = [video_info(v).frames(frame_list(i)), video_info(v).time(video_info(v).frames(frame_list(i))), L_fit,angle_fit,translation_fit]; %fit_saved;
end
end
video_info(v).fit = fit_saved;
   % end put back on
end

%% Save fits
save(strcat(path,file), 'video_info')

%% Watch video without fits

for v = 1
   for i = 1:size(video_info(v).Ibw,2)
       clf
       i
       clear Ibw
       try
              Ibw =  video_info(v).Ibw{i};
              if size(Ibw,2) > 2
       imshow(Ibw)
       pause(.5)
              end
       end
   end
end



%% Watch fits


for v = 1
    clear fit_saved
   %try
 for i = 1:size(video_info(v).fit,1)
     %clear Ibw
     %try
    Ibw =  video_info(v).Ibw{video_info(v).fit(i,1)};
     
    if size(Ibw,2) > 2
        %{
         L_guess = video_info(v).L_guess(:,i);

     pixel_scale_A = double(video_info(v).scale_nm_pix*10); %A/pixel
    translation = [video_info(v).h/2,video_info(v).w/2]*pixel_scale_A;
 angles = [45, 90, 90];
 translation = translation ;
 tic
 [~,shp2d] = nanocrystal_IrO2_shape(L_guess,angles,translation);
    Ibw_test = ibw_shape(Ibw,pixel_scale_A,shp2d);
run_time = toc
            %imshowpair(Ibw,Ibw_test)
            
           %fit_scale_down = 2^0;
    tic
run_time = toc
        %}
L_fit  = double(video_info(v).fit(i,3:9));
angle_fit = double(video_info(v).fit(i,10:12));
translation_fit = double(video_info(v).fit(i,13:14));
[~,shp2d_fit] = nanocrystal_IrO2_shape(L_fit,angle_fit,translation_fit);
    Ibw_test2 = ibw_shape(Ibw,pixel_scale_A,shp2d_fit);
        %video_info(v).fit(i,:) = [video_info(v).frames(i), video_info(v).time(video_info(v).frames(i)), L_fit,angle_fit,translation_fit];
%fit_saved(i,:) = [video_info(v).frames(i), video_info(v).time(video_info(v).frames(i)), L_fit,angle_fit,translation_fit];
       imshowpair(Ibw,Ibw_test2)
       %drawnow
       %pause(1)
%imshowpair(Ibw,Ibw_test2)
            %}
            pause(1)

    end
     %end      

 end
   %end
end

%% Compute properties from fit

v = 28;
clf
fit_trajectory = video_info(v).fit;
fit_time = video_info(v).fit(:,2);
fit_L110 = video_info(v).fit(:,3)*2
fit_L111_side1 = video_info(v).fit(:,5) + video_info(v).fit(:,8);
fit_L111_side2 = video_info(v).fit(:,6) + video_info(v).fit(:,7);
fit_L001 = video_info(v).fit(:,9);
hold on
%plot(fit_time,fit_L110/10,'o')
%plot(fit_time,fit_L001/10,'o')

%plot(fit_time,fit_L111_side1/10,'o')
%plot(fit_time,fit_L111_side2/10,'o')

plot(fit_L110(1)-fit_L110,fit_L001(1)-fit_L001,'o')

%plot(fit_L110(1)-fit_L110,fit_L111_side1(1)-fit_L111_side1,'o')
%xy_data = [fit_L110(1)-fit_L110,fit_L111_side1(1)-fit_L111_side1]

hold off
%%
% Functions start here
%}
function I_out = process_image(I_start,I_level,gain,bg_scale_down,s)

I_start = im_corr(I_start,I_level,gain,bg_scale_down);
I_start = ordfilt2(I_start,s.^2,ones(s,s));
I_out = I_start;

end


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


function Iout = crop_box_line2(x,v,f,h,w,m,s)
% x = video_info
% v = video numbe 
% f = frame number
% h = height
% w = width
% m = mean value for image intensity
% s = number of frames +- the centered frame
I_all = [];
for j = [f-s:f+s]
    try
if isempty(I_all) == 1
I_all = read_frame(x(v).file_dir(j));
else
    I_all = I_all +read_frame(x(v).file_dir(j));
end

    end
end
I = I_all; %read_frame(x(v).file_dir(f));
I = double(I);
I = I/mean(I,'all')*m;
[X,Y] = meshgrid(-(w-1)/2:(w-1)/2,-(h-1)/2:(h-1)/2);
XY = [X(:),Y(:)];
rot_mat = rotz(x(v).lines_angle(f));
XY_new = XY*rot_mat(1:2,1:2); %Rotate
line_info = x(v).lines(f,:);
line_center = [line_info(1)+line_info(3),line_info(2)+line_info(4)]/2;
XY_new = XY_new + [line_center(2),line_center(1)];
XY_new = round(XY_new);
size_I = size(I);
clear Intensity
for i= 1:size(XY_new,1)
    if XY_new(i,2) > size(I,2) | XY_new(i,2) < 1 | XY_new(i,1) > size(I,1) | XY_new(i,1) < 1 
    Intensity(i)  = m;
    else
    Intensity(i)  = I(XY_new(i,1),XY_new(i,2));
    end

end
Iout = flip(reshape(Intensity,size(X)),2);

end




function Iout = crop_box_line(x,v,f,h,w,m)
% x = video_info
% v = video numbe 
% f = frame number
% h = height
% w = width
% m = mean value for image intensity
I = read_frame(x(v).file_dir(f));
I = I/mean(I,'all')*m;
[X,Y] = meshgrid(-(w-1)/2:(w-1)/2,-(h-1)/2:(h-1)/2);
XY = [X(:),Y(:)];
rot_mat = rotz(x(v).lines_angle(f));
XY_new = XY*rot_mat(1:2,1:2); %Rotate
line_info = x(v).lines(f,:);
line_center = [line_info(1)+line_info(3),line_info(2)+line_info(4)]/2;
XY_new = XY_new + [line_center(2),line_center(1)];
XY_new = round(XY_new);
size_I = size(I);
clear Intensity
for i= 1:size(XY_new,1)
    if XY_new(i,2) > size(I,2) | XY_new(i,2) < 1 | XY_new(i,1) > size(I,1) | XY_new(i,1) < 1 
    Intensity(i)  = m;
    else
    Intensity(i)  = I(XY_new(i,1),XY_new(i,2));
    end

end
Iout = flip(reshape(Intensity,size(X)),2);

end


%reads dm4 file and outputs array in e-/(pixel)
function out = read_frame(x) %



    Imagestruc = dmread(x); 
    x_pix = Imagestruc.ImageList.Unnamed0.ImageData.Dimensions.Unnamed0.Value; %extracts x-dimension of image
    y_pix =Imagestruc.ImageList.Unnamed0.ImageData.Dimensions.Unnamed1.Value; %extracts y-dimension of image
    intscale =Imagestruc.ImageList.Unnamed0.ImageData.Calibrations.Brightness.Scale.Value ; %extracts intensity scale for image (e-/value)
    out = intscale*reshape(Imagestruc.ImageList.Unnamed0.ImageData.Data.Value,[x_pix,y_pix])'; %normalized image in e-;
end


function [sim_out] = im_similarity(Im_test,Im_actual)
% This function will calculate the similarity of an image, higher value is better
%Im_actual = imadjust(Im_actual);
%Im_actual = ordfilt2(Im_actual,1,ones(3,3));
Im_test  = imadjust(double(Im_test));
%imshowpair(Im_test,Im_actual)
%drawnow
try
    %im_diff = Im_test - Im_actual;
    
    im_diff = Im_test - Im_actual;
    t = im_diff >0; 
    im_diff(t) = 1*im_diff(t); %0.5 for different error
    %(im_diff>0) = im_diff;
    %im_diff(im_diff<0) = im_diff*1;
    im_sum = sum(abs(im_diff),'all');
    im_sum;
   % sim_out = corr2(Im_test,Im_actual);
    sim_out = -im_sum; %1./(1+im_sum);
catch
    sim_out = -10^10;
end


end


function [shp3d,shp2d] = nanocrystal_IrO2_shape(L,angles,translation)
%IrO2 parameters
spacings = [4.5051,4.5051,3.1586]; % In Angstroms %[4.4919, 4.4919, 3.1066];
unit_vectors = [1, 0, 0; 0, 1,0; 0,0,1];
b =  unit_vectors.*spacings';

%IrO2 Planes
L1 = L(1); %Side 1
L2 = L(2); %Side 2
L3 = L(3); %(111) tip 1 top
L4 = L(4); %(111) tip 2 top
L5 = L(5); %(111) tip 1 bottom
L6 = L(6); %(111) tip 2 bottom
L7 = L(7); %(001) tips
a = 1; %Integer values for for facets (aac)
c = 1; %Integer values for for facets (aac)
plane_info = [
    a, a, c, L3;
    -a, a, c, L3;
    a, -a, c, L4;
    -a, -a, c, L4;
    a, a, -c, L5;
    -a, a, -c, L5;
    a, -a, -c, L6;
    -a, -a, -c, L6;
    1, 1, 0, L1;
    -1, 1, 0, L2;
    -1, -1, 0, L1;
    1, -1, 0, L2;
    0, 0, 1, L7;
    0, 0, -1, 0];

[vertices_nanocrystal] = nanocrystal_shape2(plane_info,b);
vertices_nanocrystal = unique(vertices_nanocrystal,'rows'); %unique points
vertices_center = (max(vertices_nanocrystal)+min(vertices_nanocrystal))/2;
vertices_rotated = (vertices_nanocrystal-vertices_center)*rotz(angles(1))*roty(angles(3))*rotx(angles(2))+vertices_center;
vertices_rotated_translated = vertices_rotated + [translation,0]; %Translate particle
vertices_rotated_translated = double(vertices_rotated_translated); %Fix double error
shp3d = alphaShape(vertices_rotated_translated(:,1),vertices_rotated_translated(:,2),vertices_rotated_translated(:,3),Inf); %Output shape object
k = boundary(vertices_rotated_translated(:,1),vertices_rotated_translated(:,2),0); %indices for boundary
XY_boundary = [vertices_rotated_translated(k,1),vertices_rotated_translated(k,2)];
XY_boundary = unique(XY_boundary,'rows');
shp2d = alphaShape(XY_boundary(:,1),XY_boundary(:,2),Inf);
box_dim = max(shp2d.Points)-min(shp2d.Points);
if box_dim(1)*box_dim(2)./(L7*2*L1) < .5
%IrO2 Planes
L1 = L(1)+.0001; %Side 1
L2 = L(2)+.0001; %Side 2
L3 = L(3)+.0001; %(111) tip 1 top
L4 = L(4)+.0001; %(111) tip 2 top
L5 = L(5)+.0001; %(111) tip 1 bottom
L6 = L(6)+.0001; %(111) tip 2 bottom
L7 = L(7)+.0001; %(001) tips
a = 1; %Integer values for for facets (aac)
c = 1; %Integer values for for facets (aac)
plane_info = [
    a, a, c, L3;
    -a, a, c, L3;
    a, -a, c, L4;
    -a, -a, c, L4;
    a, a, -c, L5;
    -a, a, -c, L5;
    a, -a, -c, L6;
    -a, -a, -c, L6;
    1, 1, 0, L1;
    -1, 1, 0, L2;
    -1, -1, 0, L1;
    1, -1, 0, L2;
    0, 0, 1, L7;
    0, 0, -1, 0];

[vertices_nanocrystal] = nanocrystal_shape2(plane_info,b);
vertices_nanocrystal = unique(vertices_nanocrystal,'rows'); %unique points
vertices_center = (max(vertices_nanocrystal)+min(vertices_nanocrystal))/2;
vertices_rotated = (vertices_nanocrystal-vertices_center)*rotz(angles(1))*roty(angles(3))*rotx(angles(2))+vertices_center;
vertices_rotated_translated = vertices_rotated + [translation,0]; %Translate particle
vertices_rotated_translated = double(vertices_rotated_translated); %Fix double error
shp3d = alphaShape(vertices_rotated_translated(:,1),vertices_rotated_translated(:,2),vertices_rotated_translated(:,3),Inf); %Output shape object
k = boundary(vertices_rotated_translated(:,1),vertices_rotated_translated(:,2),0); %indices for boundary
XY_boundary = [vertices_rotated_translated(k,1),vertices_rotated_translated(k,2)];
XY_boundary = unique(XY_boundary,'rows');
shp2d = alphaShape(XY_boundary(:,1),XY_boundary(:,2),Inf);
end
end


function Ibw_out = ibw_shape(I,pixel_scale_A,shp2d)
[X,Y] = meshgrid(1:size(I,1),1:size(I,2));  %Define Coordinates
X = X.*pixel_scale_A; %In angstroms
Y = Y.*pixel_scale_A; %In angstroms
I_inshape = inShape(shp2d,X(:),Y(:));
Ibw_out = reshape(I_inshape,size(X))';

end



function [im_out,gauss_fit] = im_corr(I,multiplier,gain_scale,scaling)
Icorr = I./mean(I,'all')*multiplier;

% Compute background on small scale to save computing resources
Ismall = imresize(I,1/scaling); %Make image smaller
Ismall = Ismall./mean(Ismall,'all')*multiplier;
[X,Y] = meshgrid(1:size(Ismall,2),1:size(Ismall,1));
gauss_guess = [size(Ismall,2)/2,size(Ismall,1)/2,size(Ismall,2)/2,size(Ismall,1)/2,0]; %Angle
gauss_fit = bg_fit(Ismall,gauss_guess);
gauss = gauss_fit;
Ibg_small = rot_gauss(X,Y,gauss);
Ibg_small = Ibg_small./mean(Ibg_small,'all')*multiplier;
Ibg = imresize(Ibg_small,size(Icorr));

% Correct image
Icorr = Icorr-Ibg;
Icorr = Icorr-min(Icorr,[],'all');
% Get the mean value of the sides of the image
Icorr = Icorr./mean(Icorr,'all')*multiplier;
Icorr(Icorr>multiplier) =multiplier;
Icorr = Icorr.^gain_scale; %Change this to enhance contrast
Icorr = Icorr./mean(Icorr,'all')*multiplier;
Icorr(Icorr>multiplier) =multiplier;
im_out = Icorr;

end


function [gauss_fit] = bg_fit(I,gauss_guess)
I = I./mean(I,'all')*.5; %Normalize image
error_mask = ones(size(I)); %Weight the errors differently
error_mask(1:end,end*.3:end*.7) = 0;
[X_map,Y_map] = meshgrid(1:size(I,2),1:size(I,1));

[X,Y,Sx,Sy,angle] = ndgrid(-1:1,-1:1,-1:1,-1:1,-1:1); %[L1 and L2, L3, L4]
gauss_trial = gauss_guess; %Set up initial guess
step_size = 16;
fit_test = 0;
count_fits =  0;
while count_fits < 10
while fit_test == 0
gauss_range = [X(:),Y(:),Sx(:),Sy(:),angle(:)]*step_size;
for i = 1:243 %parfor slows down code
%Determine trial for background parameters
gauss_loop = gauss_trial +gauss_range(i,:);
Ib = rot_gauss(X_map,Y_map,gauss_loop);
Ib = Ib./mean(Ib,'all')*.5; %Normalize image

error(i) = sum(error_mask.*(I-Ib).^2,'all'); %Compute errror of fit

end
[~,minindx] = min(error); %find best guess
gauss_trial = gauss_trial + gauss_range(minindx,:);

fit_test = isequal(gauss_range(minindx,:),[0,0,0,0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = round(step_size/2);
if step_size <1
    step_size = 1; %Sets lower limit
end
end
gauss_fit = gauss_trial;

end


function [int_out] = rot_gauss(X,Y,param)
x0      = param(1);
y0      = param(2);
Sx      = param(3);
Sy      = param(4);
angle   = param(5)*pi/180;

a = cos(angle).^2./(2*Sx.^2)  + sin(angle).^2./(2*Sy.^2);
b = -sin(2*angle)./(4*Sx.^2)+ sin(2*angle)./(4*Sy.^2);
c = sin(angle).^2./(2*Sx.^2) + cos(angle).^2./(2*Sy.^2);
int_out = exp(-(a.*(X-x0).^2 + 2*b*(X-x0).*(Y-y0)+c.*(Y-y0).^2));

end

function [L_fit,angle_fit,translation_fit] = fit_IrO2(L,angles,translation,Icorr,pixel_scale_A)
[angle_fit,translation_fit] = spatial_fit(L,angles,translation,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_001_110(L,angle_fit,translation_fit,Icorr,pixel_scale_A);

for i = 1:2
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end

for i = 1:2 %:5
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_111_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_001_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end



%{
if tip_type == '111'
[L_fit] = dimension_fit_111(L,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_111(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
elseif tip_type == '001'
[L_fit] = dimension_fit_001(L,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_001(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end

for i = 1:2
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
if tip_type == '111'
[L_fit] = dimension_fit_111(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_111(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
elseif tip_type == '001'
[L_fit] = dimension_fit_001(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_001(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end
%}
%end

[~,shp2d] = nanocrystal_IrO2_shape(L_fit,angle_fit,translation_fit);
Ibw = ibw_shape(Icorr,pixel_scale_A,shp2d);
%imshowpair(Icorr,Ibw,'method','blend')


end

function [L_fit,angle_fit,translation_fit] = fit_IrO2_v2(L,angles,translation,Icorr,pixel_scale_A)

L_fit = L;
angle_fit = angles;
translation_fit = translation;
Icorr_shield = imcomplement(Icorr);
Icorr_shield(end*3/4:end,:) = 0;
Icorr_shield = imcomplement(Icorr_shield);

% Process image once so that it does not have to happen again
Icorr = imadjust(Icorr);
Icorr = imcomplement(Icorr);
%Icorr = ordfilt2(Icorr,1,ones(3,3)); %Get rid of backgroudn noise
%imshow(Icorr)

for i = 1 %:2 %:2
[L_fit,translation_fit] = dimension_fit_001_pos_v2(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,translation_fit] = dimension_fit_001_pos(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_111_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_001_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit,angle_fit] = rotation_fit_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);

%[L_fit] = dimension_fit_001_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit,translation_fit] = dimension_fit_001_pos(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit,angle_fit] = rotation_fit_global(L_fit,angle_fit,translation_fit,Icorr_shield,pixel_scale_A)
%[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr_shield,pixel_scale_A);
%[L_fit,angle_fit] = rotation_fit_global(L_fit,angle_fit,translation_fit,Icorr_shield,pixel_scale_A)
%[L_fit,translation_fit] = dimension_fit_001_pos(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr_shield,pixel_scale_A);
%[L_fit,angle_fit] = rotation_fit_global(L_fit,angle_fit,translation_fit,Icorr_shield,pixel_scale_A)

end

%{
for i = 1:3 %:5
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit] = rotation_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A)

[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end
%}



%{
if tip_type == '111'
[L_fit] = dimension_fit_111(L,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_111(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
elseif tip_type == '001'
[L_fit] = dimension_fit_001(L,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_001(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end

for i = 1:2
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
if tip_type == '111'
[L_fit] = dimension_fit_111(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_111(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
elseif tip_type == '001'
[L_fit] = dimension_fit_001(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_001(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end
%}
%end

%[~,shp2d] = nanocrystal_IrO2_shape(L_fit,angle_fit,translation_fit);
%Ibw = ibw_shape(Icorr,pixel_scale_A,shp2d);
%imshowpair(Icorr,Ibw,'method','blend')


end


function [L_fit,angle_fit,translation_fit] = fit_IrO2_v3(L,angles,translation,Icorr,pixel_scale_A)

L_fit = L;
angle_fit = angles;
translation_fit = translation;
% Process image once so that it does not have to happen again
%Icorr = imadjust(Icorr);
%Icorr = imcomplement(Icorr);
imshow(Icorr)

for i = 1:4 %:2 %:2
[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 
[L_fit] = dimension_fit_reshape_111_tip(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 

[L_fit] = dimension_fit_reshape_111_tip_part2(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 

[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit] = dimension_fit_110_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 

[L_fit] = dimension_fit_111_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 

[L_fit] = dimension_fit_001_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 

[L_fit] = dimension_fit_110_tip_global_refine(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 

[L_fit] = dimension_fit_111_tip_global_refine(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 

[L_fit] = dimension_fit_001_tip_global_refine(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 

[L_fit,translation_fit] = dimension_fit_001_pos(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 

[L_fit,translation_fit] = dimension_fit_001_pos_v2(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%L_fit = max(L_fit,0); 





%[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);

%[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit,angle_fit] = rotation_fit_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);

%[L_fit] = dimension_fit_001_tip_global(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit,translation_fit] = dimension_fit_001_pos(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[L_fit,angle_fit] = rotation_fit_global(L_fit,angle_fit,translation_fit,Icorr_shield,pixel_scale_A)
%[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr_shield,pixel_scale_A);
%[L_fit,angle_fit] = rotation_fit_global(L_fit,angle_fit,translation_fit,Icorr_shield,pixel_scale_A)
%[L_fit,translation_fit] = dimension_fit_001_pos(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
%[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr_shield,pixel_scale_A);
%[L_fit,angle_fit] = rotation_fit_global(L_fit,angle_fit,translation_fit,Icorr_shield,pixel_scale_A)

end

%{
for i = 1:3 %:5
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit] = rotation_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A)

[L_fit] = dimension_fit_001_110(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end
%}



%{
if tip_type == '111'
[L_fit] = dimension_fit_111(L,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_111(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
elseif tip_type == '001'
[L_fit] = dimension_fit_001(L,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_001(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end

for i = 1:2
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
if tip_type == '111'
[L_fit] = dimension_fit_111(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_111(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
elseif tip_type == '001'
[L_fit] = dimension_fit_001(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
[L_fit,angle_fit,translation_fit] = fit_three_001(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end
%}
%end

%[~,shp2d] = nanocrystal_IrO2_shape(L_fit,angle_fit,translation_fit);
%Ibw = ibw_shape(Icorr,pixel_scale_A,shp2d);
%imshowpair(Icorr,Ibw,'method','blend')


end



function [L_fit,angle_fit,translation_fit] = fit_IrO2_pos(L,angles,translation,Icorr,pixel_scale_A)

L_fit = L;
angle_fit = angles;
translation_fit = translation;
% Process image once so that it does not have to happen again
%Icorr = imadjust(Icorr);
%Icorr = imcomplement(Icorr);
%imshow(Icorr)

for i = 1:2 %:2 %:2
[angle_fit,translation_fit] = spatial_fit(L_fit,angle_fit,translation_fit,Icorr,pixel_scale_A);
end


end




%Fits (111) tip ratios
function [ L_out,angles_out] = rotation_fit_global(L,angles,translation,I,pixel_scale_A)
L_trial = L;
L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);
angles_test = linspace(45,90,46); %angles to test
for j = 1
for i = 1:size(angles_test,2)
L_loop =  L_trial;
angles_loop = angles_trial;
L_projected = max(abs(sqrt(2)*L_trial(1)*cos(angles_test(i)*pi/180)),abs(sqrt(2)*L_trial(1)*cos((angles_test(i)+90)*pi/180))); %Determine projected L
L_loop(1) = L_loop(1)./L_projected*L_loop(1);
L_loop(2) = L_loop(2)./L_projected*L_loop(2);
angles_loop(1) = angles_test(i); %replace with global search
try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_loop,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Ibw(end*3/4:end,:) =0; 
%imshowpair(I,Ibw)
%drawnow
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
    catch
corr_trial(i) = 0;
    end

end
[~,maxindx] = max(corr_trial); %find best guess
L_projected = max(abs(sqrt(2)*L_trial(1)*cos(angles_test(maxindx)*pi/180)),abs(sqrt(2)*L_trial(1)*cos((angles_test(maxindx)+90)*pi/180))); %Determine projected L
L_trial(1) = L_trial(1)./L_projected*L_trial(1);
L_trial(2) = L_trial(1)
angles_trial(1) = angles_test(maxindx);
end
L_out = L_trial;
angles_out = angles_trial;
end

function [angles_out,translation_out] = spatial_fit(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);


[X,Y,Z] = meshgrid(-1:1,-1:1,-1:1); %[X,Y,angle]

X = X;
Y = Y;
Z = Z/4; %Less deviation in angle
step_size = 4;
fit_test = 0;
count_fits =  0;
while count_fits < 6; %10
while fit_test == 0
XYZ_range = [X(:),Y(:),Z(:)]*step_size;
for i = 1:27 %parfor for speed
%Determine trial translation and angles
translation_loop = translation_trial +[XYZ_range(i,1),XYZ_range(i,2)];
angles_loop = angles_trial + [0,0,XYZ_range(i,3)];
%Compute shapes
try
[~,shp2d] = nanocrystal_IrO2_shape(L,angles_loop,translation_loop);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Ibw(end/2:end,:) =0; %
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
%drawnow
end
[~,maxindx] = max(corr_trial); %find best guess

translation_trial = translation_trial +[XYZ_range(maxindx,1),XYZ_range(maxindx,2)]; %Replaces old guess
angles_trial = angles_trial + [0,0,XYZ_range(maxindx,3)]; %Replaces old guess
fit_test = isequal(XYZ_range(maxindx,:),[0,0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.1
    step_size = .1; %Sets lower limit
end
end
angles_out = angles_trial;
translation_out = translation_trial;


end

%Fits (110) and (001) planes, keeps (111)/(001) ratios the same
function [L_out,translation_out] = dimension_fit_001_pos(L,angles,translation,I,pixel_scale_A)


L_trial = L;

L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);


[X,Y,L001] = meshgrid(-1:1,-1:1,-1:1); %[X,Y, 001]
step_size = 16;
fit_test = 0;
count_fits =  0;
while count_fits < 10; %10
while fit_test == 0
L_range = [X(:)*pixel_scale_A,Y(:)*pixel_scale_A,L001(:)]*step_size;
for i = 1:27 %parfor for speed
%Determine trial translation and angles
L_loop = L_trial;
translation_loop = translation_trial;
%L_loop(1) = L_loop(1)*(100+L_range(i,1))/100;
%L_loop(2) = L_loop(2)*(100+L_range(i,1))/100;
L_loop(7) = L_loop(7)*(100+L_range(i,3))/100;
translation_loop = translation_loop + [L_range(i,1),L_range(i,2)];
L_loop = [L_loop(1),L_loop(2),L_loop(7).*L_ratio(3),L_loop(7).*L_ratio(4),L_loop(7).*L_ratio(5),L_loop(7).*L_ratio(6),L_loop(7)]; %Scale (111) tips
%Compute shapes
try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_loop);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
drawnow
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial = L_trial;
%L_trial(1) = L_trial(1)*(100+L_range(maxindx,1))/100;
%L_trial(2) = L_trial(2)*(100+L_range(maxindx,1))/100;
translation_trial = translation_trial+[L_range(maxindx,1),L_range(maxindx,2)];
L_trial(7) = L_trial(7)*(100+L_range(maxindx,3))/100;
L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

fit_test = isequal(L_range(maxindx,:),[0,0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.1
    step_size = .1; %Sets lower limit
end
end
L_out = L_trial;
translation_out = translation_trial;

end

function [L_out,translation_out] = dimension_fit_001_pos_v2(L,angles,translation,I,pixel_scale_A)


L_trial = L;

L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);


[X,Y,L001] = meshgrid(-1:1,-1:1,-1:1); %[X,Y, 001]
step_size = 1;
fit_test = 0;
count_fits =  0;
while count_fits < 10; %10
while fit_test == 0
L_range = [X(:)*pixel_scale_A,Y(:)*pixel_scale_A,L001(:)*2]*step_size;
for i = 1:27 %parfor for speed
%Determine trial translation and angles
L_loop = L_trial;
translation_loop = translation_trial;
%L_loop(1) = L_loop(1)*(100+L_range(i,1))/100;
%L_loop(2) = L_loop(2)*(100+L_range(i,1))/100;
L_loop(7) = L_loop(7)*(100+L_range(i,3))/100;
translation_loop = translation_loop + [L_range(i,1),L_range(i,2)];
%L_loop = [L_loop(1),L_loop(2),L_loop(3),L_loop(4),L_loop(5),L_loop(6),L_loop(7)]; %Scale (111) tips
%Compute shapes
try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_loop);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
%drawnow
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial = L_trial;
%L_trial(1) = L_trial(1)*(100+L_range(maxindx,1))/100;
%L_trial(2) = L_trial(2)*(100+L_range(maxindx,1))/100;
translation_trial = translation_trial+[L_range(maxindx,1),L_range(maxindx,2)];
L_trial(7) = L_trial(7)*(100+L_range(maxindx,3))/100;
%L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

fit_test = isequal(L_range(maxindx,:),[0,0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.01
    step_size = .01; %Sets lower limit
end
end
L_out = L_trial;
translation_out = translation_trial;

end

function [L_out] = dimension_fit_001_110(L,angles,translation,I,pixel_scale_A)


L_trial = L;

L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; % imcomplement(I);
%Icomp = imadjust(Icomp);


[L110,L001] = meshgrid(-1:1,-1:1); %[110, 001]
step_size = 32;
fit_test = 0;
count_fits =  0;
while count_fits < 10; %10
while fit_test == 0
L_range = [L110(:),L001(:)]*step_size;
for i = 1:9 %parfor for speed
%Determine trial translation and angles
L_loop = L_trial;
L_loop(1) = L_loop(1)*(100+L_range(i,1))/100;
L_loop(2) = L_loop(2)*(100+L_range(i,1))/100;
L_loop(7) = L_loop(7)*(100+L_range(i,2))/100;
L_loop = [L_loop(1),L_loop(2),L_loop(7).*L_ratio(3),L_loop(7).*L_ratio(4),L_loop(7).*L_ratio(5),L_loop(7).*L_ratio(6),L_loop(7)]; %Scale (111) tips
%Compute shapes
try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_trial);

% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
%drawnow
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial = L_trial;
L_trial(1) = L_trial(1)*(100+L_range(maxindx,1))/100;
L_trial(2) = L_trial(2)*(100+L_range(maxindx,1))/100;
L_trial(7) = L_trial(7)*(100+L_range(maxindx,2))/100;
L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

fit_test = isequal(L_range(maxindx,:),[0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.1
    step_size = .1; %Sets lower limit
end
end
L_out = L_trial;

end
%Fits (001) tip



function [L_out] = dimension_fit_reshape_111_tip(L,angles,translation,I,pixel_scale_A)


L_trial = L;

L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; % imcomplement(I);
%Icomp = imadjust(Icomp);


[Ltip1,Ltip2] = meshgrid(-1:1,-1:1); %[110, 001]
step_size = 64;
fit_test = 0;
count_fits =  0;
while count_fits < 30; %10
while fit_test == 0
L_range = [Ltip1(:),Ltip2(:)]*step_size;
for i = 1:9 %parfor for speed
%Determine trial translation and angles
L_loop = L_trial;
%L_loop(3) = L_loop(3) + L_range(i,1);
L_loop(6) = L_loop(6) + L_range(i,1);
%L_loop(4) = L_loop(4) + L_range(i,2);
L_loop(5) = L_loop(5) + L_range(i,2);
%L_loop = [L_loop(1),L_loop(2),L_loop(7).*L_ratio(3),L_loop(7).*L_ratio(4),L_loop(7).*L_ratio(5),L_loop(7).*L_ratio(6),L_loop(7)]; %Scale (111) tips
%Compute shapes
try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_trial);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
%drawnow
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial = L_trial;
%L_trial(3) = L_trial(3) + L_range(maxindx,1);
L_trial(6) = L_trial(6) + L_range(maxindx,1);
%L_trial(4) = L_trial(4) + L_range(maxindx,2);
L_trial(5) = L_trial(5) + L_range(maxindx,2);


%L_trial(1) = L_trial(1)*(100+L_range(maxindx,1))/100;
%L_trial(2) = L_trial(2)*(100+L_range(maxindx,1))/100;
%L_trial(7) = L_trial(7)*(100+L_range(maxindx,2))/100;
%L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

fit_test = isequal(L_range(maxindx,:),[0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.1
    step_size = .1; %Sets lower limit
end
end
L_out = L_trial;

end
%Fits (001) tip


function [L_out] = dimension_fit_reshape_111_tip_part2(L,angles,translation,I,pixel_scale_A)


L_trial = L;

L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; % imcomplement(I);
%Icomp = imadjust(Icomp);


[Ltip1,Ltip2] = meshgrid(-1:1,-1:1); %[110, 001]
step_size = 64;
fit_test = 0;
count_fits =  0;
while count_fits < 30; %10
while fit_test == 0
L_range = [Ltip1(:),Ltip2(:)]*step_size;
for i = 1:9 %parfor for speed
%Determine trial translation and angles
L_loop = L_trial;
L_loop(3) = L_loop(3) - L_range(i,1);
L_loop(6) = L_loop(6) + L_range(i,1);
L_loop(4) = L_loop(4) - L_range(i,2);
L_loop(5) = L_loop(5) + L_range(i,2);
%L_loop = [L_loop(1),L_loop(2),L_loop(7).*L_ratio(3),L_loop(7).*L_ratio(4),L_loop(7).*L_ratio(5),L_loop(7).*L_ratio(6),L_loop(7)]; %Scale (111) tips
%Compute shapes
try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_trial);
% Compute BW image
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
%Compute correlation
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
%drawnow
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial = L_trial;
L_trial(3) = L_trial(3) - L_range(maxindx,1);
L_trial(6) = L_trial(6) + L_range(maxindx,1);
L_trial(4) = L_trial(4) - L_range(maxindx,2);
L_trial(5) = L_trial(5) + L_range(maxindx,2);


%L_trial(1) = L_trial(1)*(100+L_range(maxindx,1))/100;
%L_trial(2) = L_trial(2)*(100+L_range(maxindx,1))/100;
%L_trial(7) = L_trial(7)*(100+L_range(maxindx,2))/100;
%L_trial = [L_trial(1),L_trial(2),L_trial(7).*L_ratio(3),L_trial(7).*L_ratio(4),L_trial(7).*L_ratio(5),L_trial(7).*L_ratio(6),L_trial(7)]; %Scale (111) tips

fit_test = isequal(L_range(maxindx,:),[0,0]); %Checks to see if best fit leads to no change

end
count_fits = count_fits + 1;
fit_test = 0;
step_size = step_size/2;
if step_size <.1
    step_size = .1; %Sets lower limit
end
end
L_out = L_trial;

end
%Fits (001) tip


function [L_out] = dimension_fit_110_tip_global_refine(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_110_ratios = linspace(.98,1.02,51); %Ratios to use
for i = 1:size(Tip_110_ratios,2)
L_loop =  L_trial;
L_loop(1) =L_loop(1)*Tip_110_ratios(i); %Scale sides
L_loop(2) =L_loop(2)*Tip_110_ratios(i); %Scale sides

try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(1) = L_trial(1)*Tip_110_ratios(maxindx);
L_trial(2) = L_trial(2)*Tip_110_ratios(maxindx);

L_out = L_trial;

end




%Fits (111) tip ratios
function [L_out] = dimension_fit_111_tip_global_refine(L,angles,translation,I,pixel_scale_A)
L_trial = L;
L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; % imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_111_ratios = linspace(.96,1.04,51); %Ratios to use
for j = 3:6
    if j > 4
    Tip_111_ratios = linspace(.96,1.04,51); %Ratios to use

    end
for i = 1:size(Tip_111_ratios,2)
    
L_loop =  L_trial;
L_loop(j) =L_loop(j)*Tip_111_ratios(i); %Scale tip
    try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
    catch
corr_trial(i) = 0;
    end

%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(j) = L_trial(j)*Tip_111_ratios(maxindx);
end
L_out = L_trial;

end

%Fits (001) tip
function [L_out] = dimension_fit_001_tip_global_refine(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_001_ratios = linspace(.96,1.04,51); %Ratios to use
for i = 1:size(Tip_001_ratios,2)
L_loop =  L_trial;
L_loop(7) =L_loop(7)*Tip_001_ratios(i); %Scale tip
try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(7) = L_trial(7)*Tip_001_ratios(maxindx);

L_out = L_trial;

end









function [L_out] = dimension_fit_110_tip_global(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_110_ratios = linspace(.1,1.5,141); %Ratios to use
for i = 1:size(Tip_110_ratios,2)
L_loop =  L_trial;
L_loop(1) =L_trial(1)*Tip_110_ratios(i); %Scale sides
L_loop(2) =L_trial(2)*Tip_110_ratios(i); %Scale sides

try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(1) = L_trial(1)*Tip_110_ratios(maxindx);
L_trial(2) = L_trial(2)*Tip_110_ratios(maxindx);

L_out = L_trial;

end




%Fits (111) tip ratios
function [L_out] = dimension_fit_111_tip_global(L,angles,translation,I,pixel_scale_A)
L_trial = L;
L_ratio = L./L(7); %Compute all ratios compared to (001)
angles_trial = angles;
translation_trial = translation;
Icomp = I; % imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_111_ratios = linspace(.5,1.5,101); %Ratios to use
for j = 3:6
    if j > 4
    Tip_111_ratios = linspace(-.5,0.5,101); %Ratios to use

    end
for i = 1:size(Tip_111_ratios,2)
    
L_loop =  L_trial;
L_loop(j) =L_trial(7)*Tip_111_ratios(i); %Scale tip
    try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
    catch
corr_trial(i) = 0;
    end

%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(j) = L_trial(7)*Tip_111_ratios(maxindx);
end
L_out = L_trial;

end

%Fits (001) tip
function [L_out] = dimension_fit_001_tip_global(L,angles,translation,I,pixel_scale_A)
L_trial = L;
angles_trial = angles;
translation_trial = translation;
Icomp = I; %imcomplement(I);
%Icomp = imadjust(Icomp);
Tip_001_ratios = linspace(.1,1.5,141); %Ratios to use
for i = 1:size(Tip_001_ratios,2)
L_loop =  L_trial;
L_loop(7) =L_trial(7)*Tip_001_ratios(i); %Scale tip
try
[~,shp2d] = nanocrystal_IrO2_shape(L_loop,angles_trial,translation_trial);
Ibw = ibw_shape(I,pixel_scale_A,shp2d);
corr_trial(i) = im_similarity(Ibw,Icomp); %corr2(Icomp,Ibw);
catch
corr_trial(i) = 0;

end
%imshowpair(I,Ibw)
end
[~,maxindx] = max(corr_trial); %find best guess
L_trial(7) = L_trial(7)*Tip_001_ratios(maxindx);

L_out = L_trial;

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
line_pos = customWait(loi)
line_info(k,1:4) = [loi.Position(1,:),loi.Position(2,:)]; %[x1 y1 x2 y2]?
line_tip(k,1:2) = [line_info(k,1),line_info(k,2)];
line_angle(k,1) = atan2(line_info(k,4)-line_info(k,2),line_info(k,3)-line_info(k,1))*(180/pi);

end

line_info_interp = interp1(frames,line_info,framestart:frameend);
line_tip_interp = interp1(frames,line_tip,framestart:frameend);
line_angle_interp = interp1(frames,line_angle,framestart:frameend)';

% Fill in for all frames prior to the start frame 
line_info_interp_all = uint32([repelem(line_info_interp(1,:),framestart-1,1); line_info_interp]);
line_info_interp_all = single(line_info_interp_all);

line_tip_interp_all = uint32([repelem(line_tip_interp(1,:),framestart-1,1); line_tip_interp]);

line_angle_interp_all = ([repelem(line_angle_interp(1,:),framestart-1,1); line_angle_interp]);
x(j).lines = line_info_interp_all;
x(j).lines_angle = line_angle_interp_all;
x(j).lines_tip = line_tip_interp_all;
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
