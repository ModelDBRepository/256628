%Arjun Balachandar 2016
%fit_bivariate.m
%Algorithm that fits a Bivariate Normal Distribution (BND) to match
%observed proportions of neurons of each firing-type (i.e. percentage of
%observed neurons whose firing-pattern (FP) belongs to a given FP region)
%The observed proportions can be genrated by a 'target' BND, which the
%algorithm then attempts to re-create, or from inputted values of
%proportions in each region
%(Note: this code does not fit sigma (std-dev) of the distributions, instead done
%in fit_bivariate_sigma)
%See methods for details.

%***USER MODIFIABLE- refers to variables to be changed by user****

clc;
clear all;

bivargauss = @bivariable_gaussian; %used to generate BND equation

match_bivar = 1; %if 1, set target_volume to volume calculated from trial distribution, see if algorithm recreates this

loadVariables = 1;

%Data-file (generated from AutoSim.m) used to create slice that is used for fitting
load('AutoSim_istim060_distim10_ioff20_dgA0_maxgsub15'); %Data containing i_stim=60, i_off=0, 15x15 grid-size of gA vs. gsub


%size of grid (i.e. size of slice) used by algorithm - INPUT BY USER
grid_size = 10;
xmax = grid_size;
xmin = 0;
ymax = grid_size;
ymin = 0;

if num_istim > 1 %if many i_stim in the data-file, use value = inputted-by-user
    value = 60;
else
    value = max_istim; %if only one stimulation level, use that (i_stim = max_istim)
end
ind = 1;


target_volume = zeros(5);
%target_volumes = [R, SS, DO, GAP, RF] is the order of output of the firing-pattern type
target_volume = target_volume(:);
numRegions = length(target_volume);

x0 = (xmax-xmin)/2; y0 = (ymax-ymin)/2;
muX = x0; muY = y0;
mu = [muX; muY];
sigma_x = 1;
sigma_y = 1;

%***USER INPUT***
%The user inouts 'target' values of the BND, which the algorithm will try
%to determine simply from knowing the proportions of neurons in each
%FP-region that the target BND calculates.
if match_bivar == 1
    p_match = -0.6; %rho (correlation coefficient, i.e. rotation) of the distribution
    muX_match = 3; %centre (i.e. average) of the distribution
    muY_match = 4;
    sigma_x_match = 1;%sigma_x;
    sigma_y_match = 1;%sigma_y;
end

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

%create slice (of specified specified domain size 'grid_size')
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



integ = @trap_integ; %function to calculate area in a given region, by calculating double-integral using trapezoid method

%Centroid (centre of mass of each region) - calculated as shown below:
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


error_thres = 0.001; %if algorithm calculates poportions that are within error_thres of the observed proportions of the target BND, stop algorith and return optimized parameters
error = 1;
curAveComplete = 1;
dp = 0.1;
p_mini = -0.9; %initially, rho/p is varied in increments of 0.1 from -0.9 to 0.9 
p_max = abs(p_mini);
p_init = p_mini;
p = p_init;
p_min = p;
min_error = 20;

%After algorithm gets within error less than error_changetoAccu, algorithm
%changes rho/p by finer increments (0.01< from -0.99 to 0.99.
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

%after each iteration, store error (i.e. MaxError), used to show how error decreases
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
    
    %display target distribution
    Z = fun_match(X,Y);
    figure('name','Target Bivariate-fit Contour Plot');
    [hC hC] = contourf(X,Y,Z,1000);
    set(hC,'LineStyle','none');
    xlabel('gsub');
    ylabel('gA');
    title(['target vol=',mat2str(target_volume),' p_tar=',num2str(p_match)]);
end

%See methods for details
while error~=0
    %Step 1: vary rho given current uX and uY (centre of BND)
    while curAveComplete~=0 | error~=0
        if min_error<error_changeToAccu & p_init~=p_init_accu %conditions after which algorithm proceeds with more fine increments in rho, to fine-tune optimization
            p_init = p_init_accu;
            p_max = p_max_accu;
            p_mini = p_mini_accu;
            dp = dp_accu;
        end
        fun = bivargauss(p,muX,muY,sigma_x,sigma_y);
        volumes = integ(fun,x_domain,y_domain,FP_domain);
        diff = target_volume - volumes;
        max_diff = max(abs(diff)); %i.e. MaxError (maximum error in the proportions (max difference between calculated and target proportions))
        if max_diff<error_thres %if algorithm has reached target error_thres, break and return optimized parameters
            min_error = max_diff;
            p_min = p;
            min_volumes = volumes;
            error = 0;
            curAveComplete = 0;
            break;
        else
            if p < p_max %vary rho to find optimal value
                if min_error > max_diff %set new minimum error
                    min_error = max_diff;
                    p_min = p;
                    min_volumes = volumes;
                end
                p = p + dp;
            else
                curAveComplete = 0;
                display(min_volumes);
                break;
            end
        end
    end
    
    display(p_min);
    display([p_init; dp]);
    
    muX_prev = muX; muY_prev = muY;
    muX_new = muX;
    muY_new = muY;
    diff = target_volume - min_volumes;
    
    %Step 2: If after varying rho for given uX and uY, error still not
    %small enough, then move centre using centroids to calculate direction
    %vectors, which are scaled by the magnitude of error in each FP region
    %(hence move in large steps towards/from areas of large negative/positive
    %difference in calc and target proportions, and move in smaller steps
    %for regions where error is smaller).
    if error~=0
        %keep track of MaxError each iteration
        errors(iter+1) = max(abs(diff));
        uX(iter+1) = muX_new;
        uY(iter+1) = muY_new;
        ps(iter+1) = p_min;
        iters(iter+1) = iter;
        display(max(abs(diff)));
        display(iter);
        iter = iter + 1;
        
        mu = [muX; muY];
        for i=1:numRegions
            ci = [centroids(i,1); centroids(i,2)];
            muX_new = muX_new + ((centroids(i,1)-muX)/norm(ci-mu))*diff(i);
            muY_new = muY_new + ((centroids(i,2)-muY)/norm(ci-mu))*diff(i);
        end
        muX = muX_new; muY = muY_new;
        mu = [muX; muY];
        display(mu);
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

%***Display optimized parameters and plot of calculated BND***
display(min_volumes);
display(p_min);
display(mu);
fun_min = bivargauss(p_min,muX,muY,sigma_x,sigma_y);
Z = fun_min(X,Y);

figure('name','Bivariate-fit Contour Plot');
[hC hC] = contourf(X,Y,Z,1000);
set(hC,'LineStyle','none');
xlabel('gsub');
ylabel('gA');
title(['max error = ',num2str(max(abs(diff))),'; target vol=',mat2str(target_volume),' p_tar=',num2str(p_match),' mux_tar=',num2str(muX_match),' muy_tar=',num2str(muY_match),'; p_calc=',num2str(p_min),' mux_calc=',num2str(muX),' muy_calc=',num2str(muY),' calc vol=',mat2str(min_volumes)]);
%,'; calc vol=',mat2str(min_volumes)

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

