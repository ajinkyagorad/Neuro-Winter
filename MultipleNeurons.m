close all;
clear;
clc;
% Implementation of single spiking neuron with multiple inputs 
inputNeurons = 10;
outputNeurons = 10;
spikeEvents = 20;
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

for t=Tstep:Tstep:Tsim
    index = find(Tspike<=t);
    msg = sprintf('t = %d',t); disp(msg);
    
    if ~isempty(index)
        Iout = Iout+wts*Iin(:,index);
        Tspike(index) = []; % remove selected point
        Iin(:,index) = [];  % remove
    end
    
   
    % update current
    dI_dT = -Iout/tau;
    Iout = Iout+dI_dT*Tstep;
    spikedNeurons = find(Iout>0.5*inputNeurons);
    if(~spikedNeurons)
        spike(spikedNeurons) = Iout;
        Iout(spikedNeurons) = 0;
    else
        spike = 0;
    end
    
    
    
    IoutPlot = [IoutPlot Iout];
    spikedPlot = [spikedPlot spike];
    TimePlot = [TimePlot t];
    figure (1)
    plot(TimePlot,IoutPlot);
    figure(2)
    plot(TimePlot,spikedPlot);
    
    pause(0.001);
end

surf(wts)
title('weights')