%Arjun Balachandar 2016
%fit_bivariate_sigma.m
%Algorithm that fits a Bivariate Normal Distribution (BND) to match
%observed proportions of neurons of each firing-type (i.e. percentage of
%observed neurons whose firing-pattern (FP) belongs to a given FP region)
%The observed proportions can be genrated by a 'target' BND, which the
%algorithm then attempts to re-create, or from inputted values of
%proportions in each region
%Note: this code ALSO FITS SIGMA (std-dev) of the distributions
%See methods for details.

%***USER MODIFIABLE- refers to variables to be changed by user****

clc;
clear all;

%****Note: see fit_bivariate.m for details, only difference is the fitting of
%sigma (sigma_x and sigma_y). This is done by taking various combinations of sigma
%values and running the same code as in fit_bivariate.m for each
%combination, yielding the sigma values that are optimal.

bivargauss = @bivariable_gaussian;


match_bivar = 1; %if 1, set target_volume to volume calculated from trial distribution, see if algorithm recreates this

loadVariables = 1;
%target_volumes = [R, SS, DO, GAP, RF];
target_volume = [0.12, 0.4, 0.15, 0.33, 0]; %set arbitrarily for now, set to proper values later

%Data-file (generated from AutoSim.m) used to create slice that is used for fitting
load('AutoSim_istim060_distim10_ioff20_dgA0_maxgsub15'); %i=60, i_off=0-20, 15x15

grid_size = 10;

if num_istim > 1
    value = 0;
else
    value = max_istim;
end
ind = 1;

xmax = grid_size;%max_gsub;
xmin = 0;%min_gsub;
ymax = grid_size;%max_gA;
ymin = 0;%min_gA;

target_volume = target_volume(:);
numRegions = length(target_volume);

x0 = (xmax-xmin)/2; y0 = (ymax-ymin)/2;
muX = x0; muY = y0;
mu = [muX; muY];
sigma_x = 1; %sigma_x and sigma_y values, set to specific values later in code
sigma_y = 1;

%***USER MODIFIABLE: target BND parameters (including target ('_match') sigma)
if match_bivar == 1
    p_match = 0;%-0.6;
    muX_match = 3.2;
    muY_match = 4.0;
    sigma_x_match = 1.2;%sigma_x;
    sigma_y_match = 0.8;%sigma_y;
end

%number of sigmas tested (USER MODIFIABLE) from sigma_min to sigma_max by
%intervals of ds (see below)
sigma_max = 1.4;
sigma_min = 0.6;
sigma_x_max = sigma_max;
sigma_x_min = sigma_min;
sigma_y_max = sigma_max;
sigma_y_min = sigma_min;

ds = 0.1; %resolution of change in sigma
dsx = ds;
dsy = ds;

num_sx = (sigma_x_max - sigma_x_min)/dsx + 1;
num_sy = (sigma_y_max - sigma_y_min)/dsy + 1;
num_sTot = num_sx*num_sy;

del = 0.1;
dx = del;
dy = del;

num_xTot = (max_gsub-min_gsub)/dx + 1; %for modified domain size
num_x = (xmax-xmin)/dx + 1;
num_y = (ymax-ymin)/dy + 1;
x_domain = linspace(xmin,xmax,num_x);
y_domain = linspace(ymin,ymax,num_y);

[X, Y] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);

FP_domain = zeros(num_x,num_y);

cur_ind = 1;
slice_size = num_x*num_y;
param_array_slice = zeros(slice_size, 6 + nSpikeATPs);
for i=1:(num_istim)*(num_x)*(num_y)
    if param_array(i,ind) == value
        for j=1:(6+nSpikeATPs)
            param_array_slice(cur_ind,j) = param_array(i,j);
        end
        cur_ind = cur_ind + 1;
    end
end
for i=0:num_x-1
    for j=0:num_y-1
        FP_domain(j+1,i+1) = param_array_slice(j*num_xTot + i + 1,4);
    end
end


integ = @trap_integ;

centroids = zeros(numRegions,2);
numPoints = zeros(numRegions);
for i=0:num_x-1
    for j=0:num_y-1
        next_point = FP_domain(i+1,j+1) + 1;
        centroids(next_point,1) = centroids(next_point,1) + i*dx;
        centroids(next_point,2) = centroids(next_point,2) + j*dy;
        numPoints(next_point) = numPoints(next_point) + 1;
    end
end

for i=1:numRegions
    centroids(i,1) = centroids(i,1)/numPoints(i);
    centroids(i,2) = centroids(i,2)/numPoints(i);
end

%list of values from all combinations of sigma_x, sigma_y:
uXs_calc = zeros(num_sTot);
uXs_calc = uXs_calc(:);
uYs_calc = zeros(num_sTot);
uYs_calc = uYs_calc(:);
ps_calc = zeros(num_sTot);
ps_calc = ps_calc(:);
sXs_calc = zeros(num_sTot);
sXs_calc = sXs_calc(:);
sYs_calc = zeros(num_sTot);
sYs_calc = sYs_calc(:);
min_errors_calc = zeros(num_sTot);
min_errors_calc = min_errors_calc(:);
min_minerrors_calc = zeros(num_sTot);
min_minerrors_calc = min_minerrors_calc(:);
aveerrors_calc = zeros(num_sTot);
aveerrors_calc = aveerrors_calc(:);
min_volumes_calc = zeros(num_sTot,numRegions);

error_thres = 0.001;
error = 1;
curAveComplete = 1;
dp = 0.1;
p_mini = -0.9;
p_max = abs(p_mini);
p_init = p_mini;
p = p_init;
p_min = p;
min_error = 20;

error_changeToAccu = 0.003;
p_mini_accu = -0.99;
p_max_accu = abs(p_mini_accu);
p_init_accu = p_mini_accu;
dp_accu = 0.01;
volumes = zeros(numRegions);
volumes = volumes(:);
min_volumes = zeros(numRegions);
min_volumes = min_volumes(:);
min_volumes_prev = zeros(numRegions,1);
min_volumes_prev = min_volumes_prev(:);

%after each iteration, store MaxError, used to show how error decreases
%over iterations:
errors = zeros(1000);
uX = zeros(1000);
uY = zeros(1000);
ps = zeros(1000);
errors = errors(:);
uX = uX(:);
uY = uY(:);
ps = ps(:);

iters = zeros(1000);
iters = iters(:);
iter = 0;

if match_bivar == 1
    volumes_match = zeros(numRegions); %volumes under target distribution
    
    fun_match = bivargauss(p_match,muX_match,muY_match,sigma_x_match,sigma_y_match);
    volumes_match = integ(fun_match,x_domain,y_domain,FP_domain);
    
    target_volume = volumes_match; %set target volume to this distributions volumes
    display(target_volume);
    
    Z = fun_match(X,Y);
    figure('name','Target Bivariate-fit Contour Plot');
    [hC hC] = contourf(X,Y,Z,1000);
    set(hC,'LineStyle','none');
    xlabel('gsub');
    ylabel('gA');
    title(['target vol=',mat2str(target_volume),' p_tar=',num2str(p_match)]);
end

k = 0; %index - keep track of number of combinations completed
min_error_ind = 1; %index (k) where least value of max_error occurs
min_error_best = 100; %value of least min_error
min_aveerror_ind = 1;
min_aveerror_best = 100;
min_minerror_ind = 1;%index (k) where least value of min_error occurs
min_minerror_best = 100;

for i=0:num_sx
    sigma_x = sigma_x_min + i*dsx;
    for j=0:num_sy
        sigma_y = sigma_y_min + j*dsy;
        k = k + 1;
        display(k);
        display(k/num_sTot);
        display(sigma_x);
        display(sigma_y);
        
        %reset initial conditions for each sigma_x, sigma_y combo
        error_thres = 0.001;
        error = 1;
        curAveComplete = 1;
        dp = 0.1;
        p_mini = -0.9;
        p_max = abs(p_mini);
        p_init = p_mini;
        p = p_init;
        p_min = p;
        min_error = 20;

        error_changeToAccu = 0.003;
        p_mini_accu = -0.99;
        p_max_accu = abs(p_mini_accu);
        p_init_accu = p_mini_accu;
        dp_accu = 0.01;
        volumes = zeros(numRegions);
        volumes = volumes(:);
        min_volumes = zeros(numRegions);
        min_volumes = min_volumes(:);
        min_volumes_prev = zeros(numRegions,1);
        min_volumes_prev = min_volumes_prev(:);
        
        %Run algorithm for current sigma_x, sigma_y
        while error~=0
            while curAveComplete~=0 | error~=0
                if min_error<error_changeToAccu & p_init~=p_init_accu
                    p_init = p_init_accu;
                    p_max = p_max_accu;
                    p_mini = p_mini_accu;
                    dp = dp_accu;
                end
                fun = bivargauss(p,muX,muY,sigma_x,sigma_y);
                volumes = integ(fun,x_domain,y_domain,FP_domain);
                diff = target_volume - volumes;
                max_diff = max(abs(diff));
                if max_diff<error_thres%abs(diff(1))<0.001 & abs(diff(3))<0.001
                    min_error = max_diff;
                    p_min = p;
                    min_volumes = volumes;
                    error = 0;
                    curAveComplete = 0;
                    break;
                else
                    if p < p_max
                        if min_error > max_diff %set new minimum error
                            min_error = max_diff;
                            p_min = p;
                            min_volumes = volumes;
                        end
                        p = p + dp;
                    else
                        curAveComplete = 0;
                        break;
                    end
                end
            end


            muX_prev = muX; muY_prev = muY;
            muX_new = muX;
            muY_new = muY;
            diff = target_volume - min_volumes;

            if error~=0
                %keep track of MaxError each iteration
                errors(iter+1) = max(abs(diff));
                uX(iter+1) = muX_new;
                uY(iter+1) = muY_new;
                ps(iter+1) = p_min;
                iters(iter+1) = iter;
                iter = iter + 1;

                mu = [muX; muY];
                for i=1:numRegions
                    ci = [centroids(i,1); centroids(i,2)];
                    muX_new = muX_new + ((centroids(i,1)-muX)/norm(ci-mu))*diff(i);
                    muY_new = muY_new + ((centroids(i,2)-muY)/norm(ci-mu))*diff(i);
                end
                muX = muX_new; muY = muY_new;
                mu = [muX; muY];
                p = p_init;
                if abs(muX-muX_prev)<=0.001 & abs(muY-muY_prev)<=0.001
                    if p_init~=p_init_accu
                        p_init = p_init_accu;
                        p_max = p_max_accu;
                        p_mini = p_mini_accu;
                        dp = dp_accu;
                    else
                        error = 0;
                        break;
                    end
                elseif max(abs(min_volumes_prev-min_volumes))<=0.001
                    if p_init~=p_init_accu
                        p_init = p_init_accu;
                        p_max = p_max_accu;
                        p_mini = p_mini_accu;
                        dp = dp_accu;
                    else
                        error = 0;
                        break;
                    end
                end

                min_volumes_prev = min_volumes;
            else
                break;
            end
        end
        %^^end of algorithm done for each s_x, s_y pair^^
        
        %error calculated using min_error:
        if min_error_best > min_error
            min_error_best = min_error;
            min_error_ind = k;
        end
        
        uXs_calc(k) = muX; uYs_calc(k) = muY;
        ps_calc(k) = p_min;
        sXs_calc(k) = sigma_x; sYs_calc(k) = sigma_y;
        min_errors_calc(k) = min_error;
        
        %alternate way to calculate error, using average error instead of min_error:
        diff = target_volume - min_volumes;
        min_minerror = min(abs(diff)); %smallest difference in volumes between regions
        aveerror = 0;
        for i=1:numRegions
            aveerror = aveerror + abs(diff(i));
        end
        aveerror = aveerror/numRegions;
        aveerrors_calc(k) = aveerror;
        if min_aveerror_best > aveerror
            min_aveerror_best = aveerror;
            min_aveerror_ind = k;
        end
        
        min_minerrors_calc(k) = min_minerror;
        if min_minerror_best > min_minerror
            min_minerror_best = min_minerror;
            min_minerror_ind = k;
        end
        
        for l=1:numRegions
            min_volumes_calc(k,l) = min_volumes(l);
        end
        
    end
end

muX_best = uXs_calc(min_error_ind); muY_best = uYs_calc(min_error_ind);
mu_best = [muX_best; muY_best];
sigma_x_best = sXs_calc(min_error_ind); sigma_y_best = sYs_calc(min_error_ind);
sigma_best = [sigma_x_best; sigma_y_best];
p_best = ps_calc(min_error_ind);
min_error_best = min_errors_calc(min_error_ind);
min_volumes_best = zeros(numRegions);
min_volumes_best = min_volumes_best(:);
for i=1:numRegions
    min_volumes_best(i) = min_volumes_calc(min_error_ind,i);
end

display(mu_best);
display(sigma_best);
display(p_best);
display(min_error_best);
display(min_error_ind);


%if instead of min_error, aveerror is used:
muX_best2 = uXs_calc(min_aveerror_ind); muY_best2 = uYs_calc(min_aveerror_ind);
mu_best2 = [muX_best2; muY_best2];
sigma_x_best2 = sXs_calc(min_aveerror_ind); sigma_y_best2 = sYs_calc(min_aveerror_ind);
sigma_best2 = [sigma_x_best2; sigma_y_best2];
p_best2 = ps_calc(min_aveerror_ind);
aveerror_best2 = min_errors_calc(min_aveerror_ind);
min_volumes_best2 = zeros(numRegions);
min_volumes_best2 = min_volumes_best(:);
for i=1:numRegions
    min_volumes_best2(i) = min_volumes_calc(min_aveerror_ind,i);
end

display(mu_best2);
display(sigma_best2);
display(p_best2);
display(aveerror_best2);
display(min_aveerror_ind);

%if instead of min_error, min_minerror is used:
muX_best3 = uXs_calc(min_minerror_ind); muY_best3 = uYs_calc(min_minerror_ind);
mu_best3 = [muX_best3; muY_best3];
sigma_x_best3 = sXs_calc(min_minerror_ind); sigma_y_best3 = sYs_calc(min_minerror_ind);
sigma_best3 = [sigma_x_best3; sigma_y_best3];
p_best3 = ps_calc(min_minerror_ind);
min_minerror_best3 = min_errors_calc(min_minerror_ind);
min_volumes_best3 = zeros(numRegions);
min_volumes_best3 = min_volumes_best(:);
for i=1:numRegions
    min_volumes_best3(i) = min_volumes_calc(min_minerror_ind,i);
end
display(mu_best3);
display(sigma_best3);
display(p_best3);
display(min_minerror_best3);
display(min_minerror_ind);

fun_min = bivargauss(p_best,muX_best,muY_best,sigma_x_best,sigma_y_best);
Z = fun_min(X,Y);

figure('name','Bivariate-fit Contour Plot');
[hC hC] = contourf(X,Y,Z,1000);
set(hC,'LineStyle','none');
xlabel('gsub');
ylabel('gA');
title(['max error = ',num2str(max(abs(diff))),'; target vol=',mat2str(target_volume),' p_tar=',num2str(p_match),' mux_tar=',num2str(muX_match),' muy_tar=',num2str(muY_match),'; p_calc=',num2str(p_best),' mux_calc=',num2str(muX_best),' muy_calc=',num2str(muY_best),' calc vol=',mat2str(min_volumes_best)]);

%%%%%%%%%%%%%%%%%%
%Plot FP-slice
% if loadVariables == 1
%     x_size = num_x;
%     y_size = num_y;
%     gA_axis = linspace(ymin,ymax,num_y);
%     gsub_axis = linspace(xmin,xmax,num_x);
%     figure('name','Contour Plot Slice - Firing Pattern');
%     [aC aC] = contourf(FP_domain.',5);
%     title(['istim = ',int2str(value),', ioff = ',int2str(i_off)]);
%     xlabel('gsub');
%     ylabel('gA');
%     set(aC,'LineStyle','none');
%     set(gca, 'XTick', 1:x_size); % Change x-axis ticks
%     set(gca, 'XTickLabel', gsub_axis); % Change x-axis ticks labels.
%     set(gca, 'YTick', 1:y_size); % Change x-axis ticks
%     set(gca, 'YTickLabel', gA_axis); % Change x-axis ticks labels.
% end

