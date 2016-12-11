close all;
clear;
clc;
% Implementation of single spiking neuron with multiple inputs 
inputNeurons = 10;
spikeEvents = 20;
Iin = randi([0 1],spikeEvents,inputNeurons);
% Neuron Spike Matrix ( Each Row represents a set of spikes)
Tspike = randperm(200,spikeEvents)*1E-3;
Iout = 0;
wts = rand(inputNeurons,Iout);
Tstep = 1E-3;
Tsim = 200E-3;
tau = 5E-3;
IoutPlot = [];
TimePlot = [];
spikedPlot = [];
spike = 0;
for t=Tstep:Tstep:Tsim
    index = find(Tspike<=t);
    msg = sprintf('t = %d',t); disp(msg);
    IoutPlot = [IoutPlot Iout];
    spikedPlot = [spikedPlot spike];
    TimePlot = [TimePlot t];
    plot(TimePlot,IoutPlot);
    hold on;
    plot(TimePlot,spikedPlot);
    hold off;
    pause(0.001);
    if ~isempty(index)
        Iout = Iout+sum(Iin(index,:));
        Tspike(index) = []; % remove selected point
        Iin(index,:) = [];  % remove
    end
    
   
    % update current
    dI_dT = -Iout/tau;
    Iout = Iout+dI_dT*Tstep;
    
    if(Iout>0.5*inputNeurons)
        spike = Iout;
        Iout = 0;
    else
        spike = 0;
    end
end