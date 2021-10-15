%  Edited by Cathy Gu for Nucleus Tracking, Sep 2020
%  RML 2018/03/23 verion 4 corrects affine motion 
%neighbors = [time id1 x1 y1 id2 x2 y2 distance #neighbors];
%must load neighbors before running this
%output for all particles in all frames where D2min is well-defined
%[x y d2min strain dilation xaffine yaffine neigbors nborsused t dt ID actualx actualy E (x1 x2 x3 x4)]
%This is the D2min normalized by the number of neighbors. To get the
%unnormalized version: fval/npart -> fval when writing d2min matrix

% VARIABLES:

findneighborname = 'neighbors_r6_bwareaopen100.mat'; %Name of the neighbors to input
findneighbordir = '/Users/cathygu/Documents/Unet_D2min_Github_package/example/'; %Path to the neighbors file to input
neighbors = fullfile(findneighbordir, findneighborname);
load(neighbors)

savedfilename='D2min_r6'; %name of the D2min file to save
save_path='/Users/cathygu/Documents/Unet_D2min_Github_package/example/'; %path to save D2min file
m=1; %start frame, for impacts should be the 'start' of impact
dt=1; %delay time
total_frame=261; %you can analyze as many frames as you want
numbercut=2; % this is somewhat arbitrary, you may want to play with this. the point is to exclude particles that don't have enough neighbors, most importantly the ones on the edge

tic

D2min=NaN*zeros(size(neighbors,1),24);
B=0; %tracks the lost particles
C=1; % index counter

%the meat!
for t=m:dt:m+dt*(total_frame-1) %this is fine as long as no frames skipped
    disp(t)
    
    neighborst=neighbors((neighbors(:,1)==t),:);
    neighborsdt=neighbors((neighbors(:,1)==t+dt),:);
    [PIDlist,firstrow,~] = unique(neighborst(:,2),'first'); % unique particle IDs and the first time they occur in the variable neighbors

 
    for i=1:size(firstrow,1) %Loop through all particles

        clear rjt
        ri0=neighborst(firstrow(i),[3 4]); %xy for current particle        
        ri0=ri0(1,:);
        rit=neighborsdt((neighborsdt(:,2)==PIDlist(i)),[3 4]); %xy for current particle in next frame
        TF = isempty(rit);
        if TF == 0
            rit=rit(1,:); % All rows are idential, so only keep the first
            totaldisp=rit-ri0(1,:); % This is how much the current particle moved
            npart=neighborst((neighborst(:,2)==PIDlist(i)),5); %Ids of all the neighbors to this particle
            rj0=neighborst((neighborst(:,2)==PIDlist(i)),[6 7]); % xy for all the neighbors to this particle in current frame

            A=1; % Loop counter
            for n=1:size(npart,1) % for each neighbor find position in next frame
                rjtn=neighborsdt((neighborsdt(:,5)==npart(n,:)),[6 7]); % xy for neighbor in next frame
                TF = isempty(rjtn);
                if TF == 0
                    rjt(A,:)=rjtn(1,:); % rows should be identical, only keep first
                    A=A+1;
                else %Particle lost in next frame
                    rjt(A,:)=[0 0];
                    A=A+1;
                    B=B+1;%tracks the lost particles
                end
            end
 
            ind = find(rjt(:,1)); % Only those particles that were found are good
            rjt=rjt(ind,:);
            rj0=rj0(ind,:); % Get rid of bad initial positions as well
            npart=size(rj0,1); % number of neighbors?
 
            if npart > numbercut % this is somewhat arbitrary, you may want to play with this. the point is to exclude particles that don't have enough neighbors, most importantly the ones on the edge
                ri0 = repmat(ri0,npart,1); %current particle position
                rit = repmat(rit,npart,1);
                dri0=rj0-ri0; % Neighbor in current frame - current particle in current frame, d_ij(t)
                drjt=rjt-rit; % Neighbor in next frame - current particle in next framem, d_ij(t+dt)
                ddij = (drjt-dri0);
     
                d2array=@(x) sum((drjt(:,1)-(x(1)*dri0(:,1)+x(2)*dri0(:,2))).^2+(drjt(:,2)-(x(3)*dri0(:,1)+x(4)*dri0(:,2))).^2); %%%% THIS MIGHT BE WRONG?
%                d2array=@(x)
%                sum((ddij(:,1)-(x(1)*dri0(:,1)+x(2)*dri0(:,2))).^2+(ddij(:,2)-(x(3)*dri0(:,1)+x(4)*dri0(:,2))).^2);
%                %%Chen's version, Chen's E is Falk's E-I
 
                %define d2array as a function of x. minimize to get function (d2min) and x(i) values

                [x,fval] = fminsearch(d2array, [1 0 0 1]); % minimizes the d2min function and finds the location of the minimum (x) and the minimum values of D2min (fval)
                STR=x(2)+x(3); % "Local strain" as defined by Chen et al 2010
%                 DIL=((x(1)-1)*(x(4)-1)); % Dilation??  According to Chen et al dilation should actually be (x(1)+1)*(x(4)+1)-1
                DIL=((x(1)+1)*(x(4)+1))-1; % Dilation according to Chen et al 
                
                xyAFF = [x(1) x(2) ; x(3) x(4)]*dri0.'; % Calculated expected difference motion (difference between cell and neighbors in second frame)
                neighAFF = xyAFF + rit.';  % For pure affine motion, drjt is the same as the expected, drjt = neigh2-cell2new, so adding cell2new gives you the expected neighbor position in the lab frame
                
                distAFF = neighAFF - rit.' -drjt.'; % Expected minus actual positions in the lab frame for the neighbors
                xyNon = mean(distAFF,2); % What is the mean difference in the neighbor positions?
                nonaffine_angle = atan2(xyNon(2),xyNon(1));
                dxy = rit.' - ri0.'; % Actual displacement of the cell
                xyAffOnly = dxy - xyNon;
                affine_angle = atan2(xyAffOnly(2),xyAffOnly(1));
                
                
                XAFF=sum(x(1)*dri0(:,1)+x(2)*dri0(:,2))+totaldisp(1,1); % These two lines are E*dij (strain tensor times distance between neighbor and particle) -- expected motion
                YAFF=sum(x(3)*dri0(:,1)+x(4)*dri0(:,2))+totaldisp(1,2); % plus the particle's actual motion to be exact Xaff and Yaff
%                XAFF=x(1)*totaldisp(1,1)+x(2)*totaldisp(1,2); % This is the strain tensor times the particle's actual motion -- not sure what this is in this context...
%                YAFF=x(3)*totaldisp(1,1)+x(4)*totaldisp(1,2);
                
             
                D2min(C,:)=[ri0(1,1) ri0(1,2) ... % [x_current y_current
                    fval/npart STR DIL XAFF YAFF ... % D2min/n_neighbors local_strain dilation? x_affine? y_affine?
                    n npart t dt ... % n_all_neighbors n_good_neighbors current_frame dt_calc
                    PIDlist(i) totaldisp(1,1) totaldisp(1,2)... % particle_id x_disp_curr y_disp_curr]
                    x(1) x(2) x(3) x(4) xyNon(1) xyNon(2) nonaffine_angle affine_angle xyAffOnly(1) xyAffOnly(2) ]; 
                C=C+1; % advance counter
                
                
            
            end
        
        end 
    end
end

D2min(isnan(D2min(:,1)),:)=[]; % Get rid of empty/non-data rows
toc
save(strcat(save_path,savedfilename,'_bwareaopen', num2str(bwareaopen_size)),'D2min','particles','-v7.3'); % Save D2min matrix - described above for each column
      
%does not autosave these stats, save in txt file if interested 
daveraged2min=mean(D2min(:,3)) % mean d2min/n_neigh
dmediand2min=median(D2min(:,3)) % median d2min/n_neigh
dstdd2min=std(D2min(:,3)) % std of d2min/n_neigh
numbmeas=size(D2min,1) % number of data points

% %%  optional figure: D2min distribution
% %  play around with this (spacing and range), save histogram manually
% figure(1)
% %[f,x]=histogram(D2min(:,3),0:.2:10);
% histogram(D2min(:,3),60,'Normalization','probability');
% 
% %plot(x,f/trapz(x,f))
% saveas(gcf,[pathname 'D2min_distribution_2step_r' num2str(r,'%03u') '_bwareaopen' num2str(bwareaopen_size) '.png'],'png')
%  
% %% optional figure
% figure(2)
% q=100; % frame to plot
% D2minq=D2min(D2min(:,10)==q,:); % Find the rows in this frame
% scatter(D2minq(:,1),D2minq(:,2),15,(D2minq(:,3)),'filled')
% % caxis auto
% caxis ([0 20])
% daspect([1 1 1])
% h = colorbar;
% saveas(gcf,[pathname 'D2min_2step_frame' num2str(q) '_r' num2str(r,'%03u') '_bwareaopen' num2str(bwareaopen_size) '.png'],'png')  
