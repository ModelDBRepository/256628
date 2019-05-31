% modified_morris_lecar.m
% Arjun Balachandar 2016
% Adapted from Husain Shakil and PLOS comp bio (Prescott & Sejnowski 2008)
%   This program simulates the Morris Lecar model, modified to include
%   Low-Threshold K+ and A-Type K+ currents.
%   The A-Type K+ current is modelled by a modified Connor-Stevens model

%This file is called by Simulate.m, which is in turn called by AutoSim.m
%for each combination of input parameters specified by the user in AutoSim.

%[V,currents,conductances,w,z,a,b,spike,numAPs,t,ATP,ATP_spike]
function [V,currents,conductances,spike,numAPs,t] = modified_morris_lecar(i_stim, i_off, g_sub, gA, time)
    tspan = time*1000; %time is in seconds, tspan is in milliseconds
    dt = 0.1;
    loop = ceil(tspan/dt); % no. of iterations of euler

    %Initialize constants:
    C = 2;
    numAPs = 0;
    %Fast Inward Current:
    gNa = 20;
    ENa = 50;
    %Delayed Rectifier Current (IKdr or recovery variable w)
    gK = 20;
    EK = -100;
    phi_w = 0.15;
    %Slow subthreshold Potassium Current
    phi_z = 0.15;
    %Leak Current (Il):
    gl = 2;
    El = -70;
    %AHP Current:
    g_ahp = 0;
    

    %Allocate memory for variable vectors:
    t = (1:loop)*dt;
    V = zeros(loop,1);
    w = zeros(loop,1);
    z = zeros(loop,1);
    a = zeros(loop,1);
    b = zeros(loop,1);
    q = zeros(loop,1);
    currents = zeros(loop,5);
    conductances = zeros(loop,2);
    INa = zeros(loop,1);
    IK = zeros(loop,1);
    IgA = zeros(loop,1);
    Igsub = zeros(loop,1);
    spike = zeros(loop,1);
    ref = zeros(loop,1);
    
    %gahp2*q2*(v-Vk)

    %Calculate steady-state initial conditions (current = 0) using Euler method:
    for i=1:loop-1 %loop/4
          %dV/dt = (i_stim-gna*minf(V)*(V-Vna)-gk*w*(V-VK)-gl*(V-Vl)-gsub*z*(V-VK)-gA*a^4*b*(V-VK))/c
          dV_dt = (i_off - gNa*m_inf(V(i))*(V(i)-ENa) - gK*w(i)*(V(i)-EK) - gl*(V(i)-El) - gA*((a(i))^4)*b(i)*(V(i)-EK) - g_sub*z(i)*(V(i) - EK) - g_ahp*q(i)*(V(i) - EK))/C;
          V(i+1) = V(i) + dt*dV_dt; %forward Euler equation

          dw_dt = phi_w*(w_inf(V(i))- w(i))/tau_w(V(i));
          w(i+1) = w(i) + dt*dw_dt; %forward Euler equation

          dz_dt = phi_z*(z_inf(V(i)) - z(i))/tau_z(V(i));
          z(i+1) = z(i) + dt*dz_dt; %forward Euler equation
          
          da_dt = (a_inf(V(i)) - a(i))/tau_a(V(i));
          a(i+1) = a(i) + dt*da_dt; %forward Euler equation
          
          db_dt = (b_inf(V(i)) - b(i))/tau_b(V(i));
          b(i+1) = b(i) + dt*db_dt; %forward Euler equation
          
          %AHP: (not used)
          dq_dt = (q_inf(V(i)) - q(i))/tau_q(V(i));
          q(i+1) = q(i) + dt*dq_dt;
          
          INa = gNa*m_inf(V(i))*(V(i)-ENa);
          IK = gK*w(i)*(V(i)-EK);
          IgA = gA*((a(i))^4)*b(i)*(V(i)-EK);
          Igsub = g_sub*z(i)*(V(i) - EK);
          Igahp = g_ahp*q(i)*(V(i) - EK);
          Inet = INa + IK + IgA + Igsub + Igahp + gl*(V(i)-El);
          currents(i,:) = [INa, IK, IgA, Igsub, Inet];
          
          gA_current = gA*((a(i))^4)*b(i);
          gsub_current = g_sub*z(i);
          conductances(i,:) = [gA_current, gsub_current];
          

          spike(i) = (V(i) > 0).*(~ref(i));
          ref(i+1) = (V(i) > 0);
          
    end
    
    % Set initial-conditions
    V(1) = V(loop-2); %V(loop/2 - 1)
    V_rest = V(1);
    display(V_rest);
    w(1) = w(loop-2);
    z(1) = z(loop-2);
    a(1) = a(loop-2);
    b(1) = b(loop-2);
    q(1) = q(loop-2);
    
    conductances(1,:) = conductances(loop-3,:);
    currents(1,:) = currents(loop-3,:);
    
    spike(1) = spike(loop-2); % =0 %loop -2
    ref(1) = ref(loop-2); % =0
    
    
    %Euler method:
    for i=1:loop-1
    %for i=(loop/6+1):(loop-1)
          if i>loop*0.8
              i_stim = 0.0;
          end
          %dV/dt = (i_stim-gna*minf(V)*(V-Vna)-gk*w*(V-VK)-gl*(V-Vl)-gsub*z*(V-VK)-gA*a^4*b*(V-VK))/c
          dV_dt = (i_off + i_stim - gNa*m_inf(V(i))*(V(i)-ENa) - gK*w(i)*(V(i)-EK) - gl*(V(i)-El) - gA*((a(i))^4)*b(i)*(V(i)-EK) - g_sub*z(i)*(V(i) - EK)  - g_ahp*q(i)*(V(i) - EK))/C;
          V(i+1) = V(i) + dt*dV_dt; %forward Euler equation

          dw_dt = phi_w*(w_inf(V(i))- w(i))/tau_w(V(i));
          w(i+1) = w(i) + dt*dw_dt; %forward Euler equation

          dz_dt = phi_z*(z_inf(V(i)) - z(i))/tau_z(V(i));
          z(i+1) = z(i) + dt*dz_dt; %forward Euler equation
          
          da_dt = (a_inf(V(i)) - a(i))/tau_a(V(i));
          a(i+1) = a(i) + dt*da_dt; %forward Euler equation
          
          db_dt = (b_inf(V(i)) - b(i))/tau_b(V(i));
          b(i+1) = b(i) + dt*db_dt; %forward Euler equation
          
          %AHP:
          dq_dt = (q_inf(V(i)) - q(i))/tau_q(V(i));
          q(i+1) = q(i) + dt*dq_dt;
          
          INa = gNa*m_inf(V(i))*(V(i)-ENa);
          IK = gK*w(i)*(V(i)-EK);
          IgA = gA*((a(i))^4)*b(i)*(V(i)-EK);
          Igahp = g_ahp*q(i)*(V(i) - EK);
          Inet = INa + IK + IgA + Igsub + Igahp + gl*(V(i)-El);
          
          gA_current = gA*((a(i))^4)*b(i);
          gsub_current = g_sub*z(i);
          
          conductances(i,:) = [gA_current, gsub_current];
          currents(i,:) = [INa, IK, IgA, Igsub, Inet];
              

          spike(i) = (V(i) > 0).*(~ref(i));
          ref(i+1) = (V(i) > 0);
          
          if spike(i) > 0
              numAPs = numAPs + 1;
          end
    end
    
end


%Steady-state and Decay functions for gating variables:
%% Fast Sodium Channel:
function [minf] = m_inf(V)
    beta_m = -1.2;
    gamma_m = 18;
    minf = 0.5*(1+tanh((V-beta_m)/gamma_m));
end

%% Delayed Rectifier Potassium Channel:
function [winf] = w_inf(V)
    beta_w = -10;
    gamma_w = 10;
    winf = 0.5*(1+tanh((V-beta_w)/gamma_w));
end

function [tauw] = tau_w(V)
    beta_w = -10;
    gamma_w = 10;
    tauw = 1/cosh((V-beta_w)/(2*gamma_w));
end

%% AHP Potassium Channel
function [qinf] = q_inf(V)
    vhalfq = 0;
    vslopeq = -5;
    qinf = 1/(1+exp((V-vhalfq)/vslopeq));
end

function [tauq] = tau_q(V)
    tauq = 20;
end
%           q2_inf = 1/(1+exp((v-vhalfq2)/vslopeq2))
%             param gahp2=50
%          param vhalfq2=0
%             param vslopeq2=-5
%             param tau_q2=200
%             q2(0)=0

%% Slow Threshold Potassium Current (I_sub):
function [zinf] = z_inf(V)
    beta_z = -21;
    gamma_z = 15;
    zinf = 0.5*(1+tanh((V-beta_z)/gamma_z));
end

function [tauz] = tau_z(V)
    beta_z = -21;
    gamma_z = 15;
    tauz = 1/cosh((V-beta_z)/(2*gamma_z));
end

%% A-Type Potassium Current (I_A) - Modified Connor-Stevens model:
function [ainf] = a_inf(V)
    ainf = 1/(1+exp(-(V+60)/8.5));
end

function [binf] = b_inf(V)
    binf = 1/(1 + exp((V+78)/6));
end

function [taua] = tau_a(V)
    taua = 1/(exp((V + 35.82)/19.69) + exp(-(V + 79.69)/12.7)) + 0.37;
end

function [taub] = tau_b(V)
    r = V+63;
    s = -V-63;
    taub = 19*heav(r)+(1/((exp((V+46.05)/ 5)+exp(-(V+238.4)/37.45))))*heav(s);
end

%%
function [hs] = heav(x)
    if x > 0
        hs = 1;
    elseif x < 0
        hs = 0;
    else
        hs = 0.5;
    end
end
    










