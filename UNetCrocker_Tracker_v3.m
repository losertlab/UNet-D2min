% Required subfunctions:
%   - track.m
% Edited by Cathy Gu for Nucleus Tracking, Sep 2020
% Edited by Lenny Campanello for the UMD Matlab Boot Camp, June 6-10 2016
% Edited by Matt Harrington for the UMD Matlab Boot Camp, June 8-12 2015
% Created 2013/05/02 Rachel Lee
% Based on Lunhua_track_jpg (6/29/10, based on Philippe_Track), edited to 
% clean up formatting and increase commenting.

% Takes in segmented binary images for tracking
% Output is stored in variable tr, which has 8 columns: x, y, t, id,
% x_velocity, y_velocity, total_velocity, angle
% columns 5-8 are calculated based on columns 1-4, you can remove columns
% 5-8 if don't need

close all
clear
clc  

number_Frames =262 ;  % number of frames in the input image sequence
pathname = '/Users/cathygu/Documents/Unet_D2min_Github_package/example/segmented/'; % path to the segmented binary mask to perform tracking
savefolder = '/Users/cathygu/Documents/Unet_D2min_Github_package/example/'; % path to save tracks
file_short=strcat('sample');  % common name of the images
bwareaopen_size = 100;
xmin = 1;
xmax = 1024;
ymin = 1;
ymax = 1344;
param.dim=2;
param.good=8;  % how long does a track have to be
%lower_bpass=4;  
%higher_bpass=19;
%threshold=4;
%particle_size=19;
max_disp=12;      % how far are the cells allowed to move
param.mem=2;     % how long keep in memory when not found
%min_motion=2;  % minimum motion to be included in chemotactic index
param.quiet=0;
h=figure(1);
colormap('gray');
axis equal
axis off
number=1;
number_end=0;

for i=0:number_Frames
    
    filename=strcat(pathname,sprintf('%s%05u.tif',file_short,i));
    if isfile(filename)
        a =imread(filename);
        a1 = bwareaopen(a,bwareaopen_size);
        labeledDicty = bwlabel(a1);
        centroid = regionprops(labeledDicty,'Centroid');
        centroidloc = zeros(2,numel(centroid)); %centroid location(x,y) coord
        centroidx = zeros(1,numel(centroid));
        centroidy = zeros(1,numel(centroid));
        for j = 1: numel(centroid)
            centroidloc(:,j) = centroid(j).Centroid;
            centroidx(j) = centroidloc(1,j);
            centroidy(j) = centroidloc(2,j);
        end
        figure(h);
            imshow(a1);  
            axis equal
            axis off
            hold on;
            scatter(centroidx,centroidy,8,'r','o'); %plot centroids
            hold off;
            if (rem(i,10)==0)  %save every 10 frams
                savename=strcat(savefolder,sprintf('Positions_%04u',i));
                saveas(h,savename,'tif');
            end     
            number_end=number+size(centroidloc(1,:))-1;
            pos_lst(number:number_end(2),1)=centroidx.'; %x locations of particles in frame i
            pos_lst(number:number_end(2),2)=centroidy.'; %y loations of particles in frame i
            pos_lst(number:number_end(2),3)=i;  % record which frame i it is 
            number=number_end(2)+1;  
    
    end
    


end
tic
tr = track(pos_lst, max_disp,param);
toc
ntracks = max(tr(:,4));
col = repmat({'r', 'm', 'y', 'g', 'c', 'b'},1,ceil(ntracks/6));
% for i=10:10:number_Frames
%     filename=strcat(pathname,sprintf('%s%04u.tif',file_short,i));
%     a = imread(filename);
%     figure(h);
%     imagesc(a);
%     axis equal
%     axis off
%     hold on
%     tr_t = tr(tr(:,3) == i,:);
%     for j = 1:size(tr_t,1)
%         plot(tr_t(j,1),tr_t(j,2),'o','MarkerEdgeColor',col{tr_t(j,4)})
%     end
%     hold off
%     saveas(h,strcat(savefolder,sprintf('Tracks_%04u',i)),'tif');
% end
%% Determine the Velocity
tr(:,5)=NaN;   % velocity in the x direction
tr(:,6)=NaN;   % velocity in the y direction
tr(:,7)=NaN;   % total velocity
tr(:,8)=NaN;   % direction of velocity

% 2 step velocity calculation
for i=2:(length(tr(:,1))-1)
    if (tr(i,4) == tr(i-1,4)) && (tr(i,4) == tr(i+1,4)) % if same particle 
        if tr(i,3) == (tr(i-1,3)+1) && tr(i,3) == (tr(i+1,3)-1) % if no skip in time due to memory    
           tr(i,5) = (tr(i+1,1)-tr(i-1,1))/2; % x velocity is the difference in x coords
           tr(i,6) = (tr(i+1,2)-tr(i-1,2))/2; % y velocity is the difference in y coords
           tr(i,7) = (sqrt((tr(i+1,1)-tr(i-1,1))^2+(tr(i+1,2)-tr(i-1,2))^2))/2; %total velocity
           tr(i,8) = atan2(tr(i,5),tr(i,6)); %angle
        end 
    end
end

avel=nanmean(tr(:,7)); %mean velocity ignoring NaNs
figure;
vel = tr(:,7);
vel = vel(~isnan(vel));
n = histogram((vel),45,'BinLimits',[0,3],'Normalization','probability'); % create histogram for velocity distribution
xlim([0, 3]);
xlabel('Velocity')
ylabel('Frequency')
title('Histogram of Velocity', 'FontSize', 14);
h2=gcf;
saveas(h2,strcat(savefolder,'vel_figure_2step_bwareaopen_', num2str(bwareaopen_size)),'jpg');
saveas(h2,strcat(savefolder,'vel_figure_2step_bwareaopen_', num2str(bwareaopen_size)),'fig');

angle = rad2deg(tr(:,8));
angle = angle(~isnan(angle));
n2 = histogram(angle,40,'Normalization','probability');
h3=gcf;
xlabel('Direction of Motion')
ylabel('Frequency')
title('Histogram of Direction of Motion', 'FontSize', 14);
saveas(h3,strcat(savefolder,'Angle_figure_2step_bwareaopen_', num2str(bwareaopen_size)),'jpg');
saveas(h3,strcat(savefolder,'Angle_figure_2step_bwareaopen_', num2str(bwareaopen_size)),'fig');

close all;
clear h h2 h3

%% Save
save(strcat(savefolder,'tracks_2step_bwareaopen_', num2str(bwareaopen_size)));
%save(strcat(savefolder,'MSD.txt'),'MSD','-ASCII'); 
save(strcat(savefolder,'tracks_2step_bwareaopen_', num2str(bwareaopen_size), '.txt'),'tr','-ASCII'); 