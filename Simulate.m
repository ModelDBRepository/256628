% Simulate.m
% Arjun Balachandar 2016
% Simulate neuron: inputs and monitor stimulation parameters

%AutoSim.m calls this file, inputting parameters i (i_stim), i_off
%(pre-stim stimulation to alter resting Vm), time of stimulation,
%(ATP-related input parameters not used here) and if this specific
%neuron-model is one of many, or is the only one being simulated by AutoSim
%(i.e. multipleNeurons: if 0, then display single-neuron-specific data,
%e.g. V-t traces)
function [firing_info] = Simulate(i,i_off,gsub,g_A, tim, multipleNeurons)
    check = @modified_morris_lecar; %uses modified_morris_lecar model from 

    %initialize variables
    type = 1;
    type_str = 'R'; %stores firing-pattern outputted by model neuron after simulation
    i_stim = i;
    gA = g_A;
    g_sub = gsub;
    time = tim;
    i_offset = i_off;
    
    %run modified morris lecar model using input parameters, and store
    %output variables (V-t, currents-t, conductances-t data etc..)
    [V,currents,conductances,spike,numAPs,t] = check(i_stim, i_offset, g_sub, gA, time);

    %I-t data after simulation (i.e. current through each channel)
    INa = currents(:,1);
    IK = currents(:,2);
    IgA = currents(:,3);
    Igsub = currents(:,4);
    Inet = currents(:,5);
    
    %instantaneous (i.e. current) value of g (conductance) at each
    %time-step:
    gA_current = conductances(:,1);
    gsub_current = conductances(:,2);
    

    %%Determine Firing Pattern
    spike_times = zeros(numAPs,1);
    curAP = 1; %index of AP currently being analyzed
    for i=1:length(spike) %spike stores spike-time of output response after simulation
        if spike(i) > 0
            spike_times(curAP) = t(i);
            curAP = curAP + 1;
        end
    end

    rate = numAPs/time; %calculate firing-rate
    
    %Classify neuron firing-pattern type
    %type corresponds to numrical classification system, while type_str is the
    %corresponding string-name
    if numAPs == 1
        type = 1;
        type_str = 'SS'; %if one spike, then single-spike neuron (SS)
    elseif numAPs < 1
        type_str = 'R'; %if no spikes, than reluctanct/non-spiking neuron (R)
        type = 0;
    elseif numAPs >= 3 %if more than 2 spikes, than can be of delayed-onset, gap- (Gap) or repetitive- (R - i.e. tonic in paper) spiking
        ISI2 = spike_times(3) - spike_times(2); %inter-spike interval between spikes 3 and 2
        ISI1 = spike_times(2) - spike_times(1); %ISI between t=0 and spike 1
        if 0.5*spike_times(2) > spike_times(1) %if the first spike does not occur late enough, then not considered a delayed-onset neuron
            %^specifically, if 0.5x the time for the second spike occurs after
            %the first spike occurs, then first spike did not occur late
            %enough to be 'delayed-onset'
            if ISI1 > 1.5*ISI2
                type_str = 'Gap';
                %If not delayed-onset, and if time between first and second
                %spike is more than 1.5x the time between 2nd and 3rd spikes,
                %then the initial gap between 1st and 2nd spikes is
                %considered long enough to classify neuron as 'gap-spiking'
                type = 3;
            else
                type_str = 'RF'; %if not a gap neuron, then considered RF
                type = 4;
            end
        else
            type_str = 'DO'; %delayed-onset'
            type = 2;
        end
    end
    
    %store all pertinent information in firing_info, for eventual use when
    %analyzing simulation data (Sim_Analysis.m)
    firing_info = [];
    firing_info(1) = type;
    firing_info(2) = rate;
    
    if multipleNeurons == 0
        
        %Plots --> if this neuron-model is the only one being simulated by
        %AutoSim.m, then display all pertinent output data to user
        
        %Conductance plots:
        figure('name','gA-t Plot');
        plot(t,gA_current);
        
        figure('name','gA-t (normalized) Plot');
        plot(t,gA_current*(1.0/gA));
        axis([0 tim*1000 0 1]);
        
        figure('name','gsub-t Plot');
        plot(t,gsub_current);
        
        figure('name','gsub-t (normalized) Plot');
        plot(t,gsub_current*(1.0/g_sub));
        axis([0 tim*1000 0 1]);
        
        figure('name','Conductances Plot');
        plot(t,gA_current,t,gsub_current);
        legend('gA','gsub');
        
        %Plot showing activation of each channel (normalized to the maximum
        %value of g)
        figure('name','Conductances (normalized) Plot');
        plot(t,gA_current*(1.0/gA),t,gsub_current*(1.0/g_sub));
        legend('gA/gA_max','gsub/gsub_max');
        axis([0 tim*1000 0 1]);
        
        figure('name','gA-V plot');
        plot(V,gA_current);
        
        %Plot showing (semi-real-time, to give impression of dynamic plot)
        %how the instantaneous conductance gA corresponds to V, over time:
%         figure('name','gA-V plot (real-time)');
%         hold on;
%         for k = 1:size(V)-1
%             plot(V(k:k+1),gA_current(k:k+1)*(1.0/gA));
%             pause(0.0001);
%         end
%         hold off;
        
        figure('name','gsub-V');
        plot(V,gsub_current);
        
        %V-t plots:
        
        figure('name','V-t Plot');
        plot(t,V);
        axis([0 tim*1000 -80 50]);
        file_type = 'V-t_';
        
    end
    
end



