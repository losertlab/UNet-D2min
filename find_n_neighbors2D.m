%  Edited by Cathy Gu for Nucleus Tracking, Sep 2020
%  RML 2018/03/23

clc
clear
close all

%must load tracks data before running this

%This script takes track files in the 4 column format: x,y,t,ID. Any extra
%columns are fine and can be ignored. It then identifies
%neighbors in each frame. 

%The output matrix 'neighbors' is large and somewhat redundant. It is the
%only thing that is saved, aside from the original track file. 
%Each row corresponds to a particular particle (i) and its neighbor (j) at a particular time. 
%C1 is time. Particle i and j are next IDed and located. 
%C8 gives the distance between them. 
%C9 gives the number of neighbors that particle i has in that frame. 

% VARIABLES:
r=6; %how many nearest neighbors to use
Nadd = r; % two names, useful to keep apart for now
loadname = 'tracks_2step_bwareaopen_100.mat'; %name of the tracks file
rootdir ='/Users/cathygu/Documents/Unet_D2min_Github_package/example/'; % path to the tracks file
fid = fullfile(rootdir, loadname);
load(fid);
particles = tr;
pathname='/Users/cathygu/Documents/Unet_D2min_Github_package/example/'; %path to save neighbor file

mintime=1;
maxtime=max(particles(:,3));
% frames=particles(:,[4 3 1 2]); %modify this row so it picks out the four relevant columns 
frames=particles(:,[3 4 1 2]); %modify this row so it picks out the four relevant columns 
%in this order: time, ID, x, y

frames(frames(:,1)>maxtime,:)=[]; %kill those frames
frames(frames(:,1)<mintime,:)=[]; %kill those frames

%preallocation
neighbors=zeros(100000000,9);
A=1;   %this will be used to count up neighbors
N_ref = 0;
averagenumber=0;
%THE MEAT: for each of the j particles in frame k, find all neighbors i within radius r
tic
for k=mintime:maxtime %which frame
    framesk=frames((frames(:,1)==k),:);
    TF=isempty(framesk);
    if TF == 0
        disp(k)
        for j=1:size(framesk,1) %which reference particle        
            deltad = ( (framesk(j,3)-framesk(:,3)).^2 + (framesk(j,4)-framesk(:,4)).^2 ).^.5;
            [~,indd] = sort(deltad);
%             % sanity check
%             if indd(1) ~= j
%                 error('Shortest distance should be to yourself -- check math')
%             end
            
%             if length(indd)-1 < r
%                 Nadd = length(indd) - 1;
%                 disp('Not enough neighbors found')
%             else
%                 Nadd = r;
%             end
            
            neighbors(A:A+Nadd-1,:) = [k*ones(Nadd,1) framesk(j,2)*ones(Nadd,1) framesk(j,3)*ones(Nadd,1) framesk(j,4)*ones(Nadd,1) ... %time id1 x1 y1
                framesk(indd(2:Nadd+1),2) framesk(indd(2:Nadd+1),3) framesk(indd(2:Nadd+1),4)... id2 x2 y2
                deltad(indd(2:Nadd+1)) Nadd*ones(Nadd,1)]; % distance N_neighbors ERK]; %%%%% ADDED ERK FOR THIS PROJECT
            
            A = A+Nadd;
            
        end
    end
end
toc
neighbors(neighbors(:,9)==0,:)=[];

% %% figures if you like!
% time=1;
% neighborplot=neighbors(neighbors(:,1)==time,:);
% figure(1)
% scatter(neighborplot(:,3),neighborplot(:,4),[],neighborplot(:,9),'filled')
% set(gca,'DataAspectRatio',[1 1 1],'FontSize',16)
% h = colorbar;
% set(get(h,'Label'),'String','Number of Neighbors','FontSize',16)
% axis([min(neighborplot(:,3))-5 max(neighborplot(:,3))+5 min(neighborplot(:,4))-5 max(neighborplot(:,4))+5])
% % saveas(gcf,[pathname 'NumberNeighbors_Frame' num2str(time,'%03u') '_r' num2str(r,'%03u') '.png'],'png')


%% SAVE

 clearvars -except particles neighbors pathname r bwareaopen_size averagenumber name
 save(strcat(pathname,'neighbors_r', num2str(r),'_bwareaopen100'));