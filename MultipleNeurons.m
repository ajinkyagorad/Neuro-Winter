close all;
clear;
clc;
% Implementation of single spiking neuron with multiple inputs 
inputNeurons = 10;
outputNeurons =2;
spikeEvents = 10;
% Neuron Spike Matrix ( Each Row represents a set of spikes)
Iin = randi([0 1],inputNeurons,spikeEvents);
Tspike = randperm(200,spikeEvents)*1E-3;
Iout = zeros(outputNeurons,1);
wts = rand(outputNeurons,inputNeurons);
Tstep = 1E-3;
Tsim = 200E-3;
tau = 5E-3;
IoutPlot = [];
TimePlot = [];
spikedPlot = [];
spike =zeros(outputNeurons,1);
tLastSpike =0;
Iinput = 0;
figHandle = figure ('Position',[100,100,1049,895]);
for t=Tstep:Tstep:Tsim
    index = find(Tspike<=t);
    msg = sprintf('t = %d',t); disp(msg);
    if ~isempty(index)
        Iinput = wts*Iin(:,index);
        Iout = Iout+wts*Iin(:,index);
        Tspike(index) = []; % remove selected point
        Iin(:,index) = [];  % remove
        tLastSpike = t;
    end
    
   
    % update current
    %dI_dT = -Iout/tau; % Linear
    tp = t-tLastSpike;
    dI_dT = -Iout/tau+Iinput/tau/(1+tau).*exp(-tp/tau).*(1-tp/tau);
    Iout = Iout+dI_dT*Tstep;
    
    % Code Simulation Error Checking
%     err = find(Iout<0);
%     if(~isempty(err))
%         Iout(err) = 0;
%     end
    
    % check for spiking
    spikedNeurons = find(Iout>0.3*inputNeurons);
    if(~isempty(spikedNeurons))
        msg = sprintf('spiked %s at t = %d\n',48+spikedNeurons,t); disp(msg);
        spike(spikedNeurons) = Iout(spikedNeurons);
        Iout(spikedNeurons) = 0;
    else
        spike = zeros(outputNeurons,1);
    end
    
    
    
    
    IoutPlot = [IoutPlot Iout];
    spikedPlot = [spikedPlot spike];
    TimePlot = [TimePlot t];
    
    subplot(1,2,1);    plot(TimePlot,IoutPlot);
    subplot(1,2,2);    plot(TimePlot,spikedPlot);

    pause(0.001);
end

surf(wts);
title('weights')