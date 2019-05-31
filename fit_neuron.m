%Arjun Balachandar 2017
%fit_neuron.m
%This program takes the observed firing-pattern of a given neuron at
%various stimulation intensities i_stim, and uses that in conjunction with the
%generated parameter-space slices to determine the possible underlying
%conductance values of gA and gsub characterizing the neuron.
%Speficially, this is done by loading 2 data sets containting 2 ranges of
%values of i_stim (and the FP of the neurons at each slice of i_stim for each combination of gA and gsub) and
%finding areas (conductance combinations) of the parameter-space where the neuron matches both sets of data.
%If no such area found, the algorithm returns no area.
%(See methods).

%***USER MODIFIABLE- refers to variables to be changed by user****

clc;
clear all;

%load simulation data
load('AutoSim_istim050_distim5_ioff0_dgA0.mat');

grid_size = 20;

xmax = grid_size;
xmin = 0;
ymax = grid_size;
ymin = 0;

del = 0.1;
dx = del;
dy = del;

num_x = (xmax-xmin)/dx + 1;
num_y = (ymax-ymin)/dy + 1;
FP_domain = zeros(num_istim,num_y,num_x);

best_matches1 = zeros(num_x*num_y,2);
best_matches2 = zeros(num_x*num_y,2);
num_best1 = 0; num_best2 = 0;
perc_match_best1 = 0; perc_match_best2 = 0;

%***USER MODIFIABLE: user inputs the FP of the neuron at each level of i_stim
%(see line below), which the algorithm uses to find areas of best match.
%I_stim order: [50,55,60,65,70,75,80]; from the cumulative data loaded from both data
%sets
%FP numbering system: [R, S, D, G, T] = [0, 1, 2, 3, 4] (i.e. [reluctant-,
%single, delayed-, gap- and tonic-spiking] )
target_tot = [0,0,2,3,3,3,3];
%Other example combinations of FPs for a neuron:
%[0,0,0,1,1,1,1]; ; %[0,2,2,2,3,3,3]; 
%[0,0,1,1,3,3,4]; [0,0,1,3,3,3,4]; [0,1,1,1,3,3,4]; [0,0,0,1,1,3,3];
%[0,2,2,3,3,3,3]; [0,2,3,3,4,4,4]; %[0,2,3,3,3,4,4]; [0,2,3,3,3,3,4]; [0,2,3,3,3,3,3];
%[0,0,0,0,1,3,3];[0,0,2,3,3,3,3];
target = target_tot(1:4); %take first part of target_tot initially (to compare to first data-set), then second part when comapring to second data-set

for t=1:2
    if t==2 %Load data for higher i_stims
        load('AutoSim_distim5_ioff0_dgA0'); %Contains i=70
        target = target_tot(5:7);
    end
    
    values = zeros(num_istim);
    values = values(:);
    
    display(num_istim);
    
    %value above which perc_match must be in order to be considered as a match:
    match_thres = (num_istim-1)/num_istim;
    if match_thres < 0.7
        match_thres = 0.7;
    end
    
    ind = 1;
    
    %location of best match (OUTPUT OF THE ALGORITHM):
    perc_match_best = 0.0;
    x_best = 0.0; y_best = 0.0;
    
    num_xTot = (max_gsub-min_gsub)/dx + 1; %for modified domain size
    num_x = (xmax-xmin)/dx + 1;
    num_y = (ymax-ymin)/dy + 1;
    x_domain = linspace(xmin,xmax,num_x);
    y_domain = linspace(ymin,ymax,num_y);
    
    [X, Y] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);
    
    FP_domain_slice = zeros(num_y,num_x);
    match_level = zeros(num_y,num_x);
    
    slice_size = num_x*num_y;
    param_array_slice = zeros(slice_size, 6 + nSpikeATPs);
    
    for k=1:num_istim
        %calculate each slice from i_stim = min_istim to i_stim max_istim,
        %and plot them for the user.
        value = min_istim + d_istim*(k-1);
        cur_ind = 1;
        
        for i=1:(num_istim)*(num_x)*(num_y)
            if param_array(i,ind) == value %current i_stim value
                for j=1:(6+nSpikeATPs)
                    param_array_slice(cur_ind,j) = param_array(i,j);
                end
                cur_ind = cur_ind + 1;
            end
        end
        
        for i=0:num_x-1
            for j=0:num_y-1
                FP_domain(k+4*(t-1),j+1,i+1) = param_array_slice(j*num_xTot + i + 1,4);
                FP_domain_slice(j+1,i+1) = param_array_slice(j*num_xTot + i + 1,4);
            end
        end
        
        x_size = num_x;
        y_size = num_y;
        gA_axis = linspace(ymin,ymax,num_y);
        gsub_axis = linspace(xmin,xmax,num_x);
        figure(k+4*(t-1));
        [aC aC] = contourf(FP_domain_slice.',5);
        
        title(['istim = ',int2str(value),', ioff = ',int2str(i_off)]);
        xlabel('gsub');
        ylabel('gA');
        set(aC,'LineStyle','none');
        set(gca, 'XTick', 1:x_size); % Change x-axis ticks
        set(gca, 'XTickLabel', gsub_axis); % Change x-axis ticks labels.
        set(gca, 'YTick', 1:y_size); % Change x-axis ticks
        set(gca, 'YTickLabel', gA_axis); % Change x-axis ticks labels.
        
        FP_domain_slice = [];
        FP_domain_slice = zeros(num_y,num_x);
    end
    
    num_best = 0;
    best_matches = [];
    %finds areas (combinations of cnductances) such that the neuron
    %firing-pattern at each i_stim (for that combination of conductances)
    %matches the FP seen in the parameter-space at that combination of
    %conductances:
    for i=1:num_x
        for j=1:num_y
            perc_match = 0;
            
            for k=1:num_istim
                if target(k) == FP_domain(k+4*(t-1),j,i)
                    perc_match = perc_match + 1;
                end
            end
            
            perc_match = perc_match/num_istim;
            if perc_match > perc_match_best && perc_match >= match_thres
                perc_match_best = perc_match;
                num_best = 1;
                x_best = i;
                y_best = j;
                best_matches = [];
                best_matches = zeros(num_x*num_y,2);
                best_matches(num_best,1) = x_best;
                best_matches(num_best,2) = y_best;
            elseif perc_match == perc_match_best && perc_match >= match_thres
                x_best = i;
                y_best = j;
                num_best = num_best + 1;
                best_matches(num_best,1) = x_best;
                best_matches(num_best,2) = y_best;
            end
            match_level(j,i) = perc_match;
        end
    end
    
    if t==1
        num_best1 = num_best;
        best_matches1 = best_matches;
        perc_match_best1 = perc_match_best;
    else
        num_best2 = num_best;
        best_matches2 = best_matches;
        perc_match_best2 = perc_match_best;
    end

end

%areas where the areas of best fit in both i_stim sections match
best_matches_final = intersect(best_matches1,best_matches2,'rows');
best_matches_final1 = zeros(num_best1,2);
best_matches_final2 = zeros(num_best2,2);

perc_match_best = (perc_match_best1*4 + perc_match_best2*3)/7.0;
display(perc_match_best);

num_most = max(num_best1,num_best2);
for i=1:num_best1
    best_matches_final1(i,1) = xmin + (best_matches1(i,1)-1)*dx;
    best_matches_final1(i,2) = ymin + (best_matches1(i,2)-1)*dy;
end
x_size = num_x;
y_size = num_y;
gA_axis = linspace(ymin,ymax,num_y);
gsub_axis = linspace(xmin,xmax,num_x);

%Display the areas of best fit in section 1, section 2, and both section 1
%and 2 (i.e. the real areas of best fit)

figure('name','Area of best fit1');
scat = scatter(best_matches_final1(:,2),best_matches_final1(:,1),10,'filled');

xlabel('gsub');
ylabel('gA');
xlim([xmin xmax]);
ylim([ymin ymax]);

for i=1:num_best2
    best_matches_final2(i,1) = xmin + (best_matches2(i,1)-1)*dx;
    best_matches_final2(i,2) = ymin + (best_matches2(i,2)-1)*dy;
end

figure('name','Area of best fit2');
scat = scatter(best_matches_final2(:,2),best_matches_final2(:,1),10,'filled');

xlabel('gsub');
ylabel('gA');
xlim([xmin xmax]);
ylim([ymin ymax]);

if size(best_matches_final) > 0
    perc_match_best = (perc_match_best1*4 + perc_match_best2*3)/7.0;
else
    perc_match_best = 0.0
end

for i=1:size(best_matches_final)
    best_matches_final(i,1) = xmin + (best_matches_final(i,1)-1)*dx;
    best_matches_final(i,2) = ymin + (best_matches_final(i,2)-1)*dy;
end

figure('name','Area of best fit - Total');
scat = scatter(best_matches_final(:,2),best_matches_final(:,1),10,'filled');

xlabel('gsub');
ylabel('gA');
xlim([xmin xmax]);
ylim([ymin ymax]);
            
