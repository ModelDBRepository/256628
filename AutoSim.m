%AutoSim.m
%Arjun Balachandar 2016
%***MAIN CODE TO RUN SIMULATION***
%This program automates simulation parameters i_stim, gA & g_sub.
%Using each combination of simulation parameters, a model neuron (using
%Simulate.m) is created, which outputs the firing-response due to the input
%conditions/parameters.
%Note: gA refers to A-type potassium channel, while g_sub/gK_lt corresponds to the
%subthreshold non-inactivating potassium channel


clear all; clc;
model = @Simulate;

%If custom = 1, use custom resting Vm by applying a pre-stimulus current I_off, but keep net current during stim constant, equal to I_stim
%If custon = 0, run as normal (no change in resting Vm)
custom = 0;

%Stimulation parameter ranges:
time = 0.500; %length of stimulation (in sec)
i_off = 0; %level of pre-stimulus current used to change Vm

%***i_stim between max_istim and min_istim, by intervals of d_istim (see below) are used for simulation
max_istim = 50;%maximum stimulus intensity
min_istim = 50;%minimum level of stimulus intensity

%***gA/gsub between max_g and min_g, by intervals of d_g (see below) are used for simulation
max_gA = 2;
min_gA = 2;
max_gsub = 8;
min_gsub = 8;

%intervals used for varying parameters i_stim, gA and gsub
d_istim = 10;
d_gA = 0.1;
d_gsub = 0.1;

%using above ranges of parameters set by user, determine number of
%combinations of each respective parameter
num_istim = (max_istim - min_istim)/d_istim + 1;
num_gA = (max_gA - min_gA)/d_gA + 1;
num_gsub = (max_gsub - min_gsub)/d_gsub + 1;

if custom == 1 %If in i_off mode, meaning before actual stimulation a current is applied to change resting Vm
    d_ioff = 10;
    max_ioff = 20;
    min_ioff = 20;
    num_ioff = (max_ioff - min_ioff)/d_ioff + 1;
    num_istim = num_ioff;
end

multipleNeurons = 1; %by deafult 1: if 1, multiple model-neuron simulations; if 0, only one simulation, program displays neuron-specific plots, such as V-t plot
if custom==1
    if num_ioff <= 1 & num_gA <=1 & num_gsub <=1
        multipleNeurons = 0;
    end
else
    if num_istim <= 1 & num_gA <=1 & num_gsub <=1
        multipleNeurons = 0;
    end
end

%3D-Array to store all pertinent simulation-specific input & output variables, including input gA and gsub, output firing-pattern and rate:
param_array = zeros((num_istim)*(num_gA)*(num_gsub),6 + nSpikeATPs); 

for i=0:(num_istim-1)
    i_stim = min_istim + d_istim*i;
    if custom == 1
        i_stim = min_istim - d_ioff*i
        i_off = min_ioff + d_ioff*i
    end
    display(i_stim); %Displays status of simulation to user (i.e. what percentage of simulations completed)
    display(i/num_istim);
    for j=0:(num_gsub-1)
        g_sub = min_gsub + d_gsub*j;
        for k=0:(num_gA-1)
            gA = min_gA + d_gA*k;
            firing_pattern = model(i_stim,i_off,g_sub,gA,time,multipleNeurons); %run model for given combination of input parameters
            
            %Parameter array for output to file
            if custom == 1
                param_array(i*(num_gsub*num_gA) + j*num_gA + k + 1, 1) = i_off;
            else
                param_array(i*(num_gsub*num_gA) + j*num_gA + k + 1, 1) = i_stim;
            end
            param_array(i*(num_gsub*num_gA) + j*num_gA + k + 1, 2) = g_sub;
            param_array(i*(num_gsub*num_gA) + j*num_gA + k + 1, 3) = gA;
            param_array(i*(num_gsub*num_gA) + j*num_gA + k + 1, 4) = firing_pattern(1); %Firing-pattern (i.e. T, S, D, G, R)
            param_array(i*(num_gsub*num_gA) + j*num_gA + k + 1, 5) = firing_pattern(2); %Firing-rate
        end
    end
end

%Save variables for later analysis with 'Sim_Analysis' code
if num_gA > 1 & num_gsub > 1
    if d_gA < 0.5
        if max_gsub >= 20
            save(['AutoSim_istim0',int2str(min_istim),'_distim',int2str(d_istim),'_ioff',int2str(i_off),'_dgA',int2str(d_gA),'.mat']);
        else
            save(['AutoSim_istim0',int2str(min_istim),'_distim',int2str(d_istim),'_ioff',int2str(i_off),'_dgA',int2str(d_gA),'_maxgsub',int2str(max_gsub),'.mat']);
        end
    else
        save(['AutoSim_distim',int2str(d_istim),'_ioff',int2str(i_off),'.mat']);
    end
    save('AutoSim.mat');
end




            