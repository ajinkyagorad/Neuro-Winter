close all;
clear;
clc;
% Implementation of single spiking neuron with multiple inputs 
inputNeurons = 12;
outputNeurons =2;
spikeEvents = 20;
% Neuron Spike Matrix ( Each Row represents a set of spikes)
% Iin = randi([0 1],inputNeurons,spikeEvents);      % Random Input
Iin = repmat( [ 1 1 1 1 0 0 0 0 1 1 1 1]',[1 spikeEvents]);
Tspike = sort(randperm(200,spikeEvents))*1E-3;    % Times for input Iin spikes (must be unique)
msg = sprintf('Time-Input Spikes'); disp(msg);
for i = 1:length(Tspike)
    disp([ Tspike(i) Iin(:,i)']);
    
end
Iout = zeros(outputNeurons,1);              % Output Current as time
wts = rand(outputNeurons,inputNeurons);     % Weights matrix for input and output neuron
Tstep = 1E-3;
Tsim = 200E-3;
tau = 5E-3;
IoutPlot = [];
TimePlot = [];
spikedPlot = [];

PostSpike =zeros(outputNeurons,1);      % spikes from output Neuron ( To analyse output spikes )
tPostSpike = zeros(outputNeurons,1);    % time of output Neuron Spike
tPreSpike = zeros(inputNeurons,1);      % time of spikes of input neurons ( For  STDP)
Iinput = zeros(inputNeurons,1);         % Input Current variable
figHandle = figure ('Position',[10,10,1049,565]);   % For  customizing display
%% Simulation Starts here
for t=Tstep:Tstep:Tsim
    index = find(Tspike<=t);        % Find the spike at current time ( the event which has spiked is later removed from the list)
    msg = sprintf('t = %d',t); disp(msg);
    if ~isempty(index)              % If spiked (ie we have non zero input)
        Iinput = Iin(:,index);
        
        Iout = Iout+wts*Iinput;
        Tspike(index) = []; % remove selected point
        Iin(:,index) = [];  % remove
        tPreSpike(Iinput>0) = t;   % save time for input neuron spike
        %apply LTD here
        
    end
    
   
    % update current
    dI_dT = -Iout/tau; % Leaky integrate and fire model
   % tp = t-tLastSpike; % alpha function model ( not working as required)
   % dI_dT = Iout/tau/(1+tau).*exp(-tp/tau).*(1-tp/tau);
    Iout = Iout+dI_dT*Tstep;        
    
   
    % check for spiking accoring to threshold
    spikedNeurons = find(Iout>0.3*inputNeurons);    
    if(~isempty(spikedNeurons))
        msg = sprintf('spiked %s at t = %d\n',48+spikedNeurons,t); disp(msg);
        PostSpike(spikedNeurons) = Iout(spikedNeurons);
        Iout(spikedNeurons) = 0;    % reset the spike  neuron current to zero
        tPostSpike(spikedNeurons) = t;     % save the spiking time
        
        %apply LTP here
        
        
    else
        PostSpike = zeros(outputNeurons,1);
    end
    
    %Weights Update
    tPostSpikeM = repmat( tPostSpike, [1 length(tPreSpike)]);
    tPreSpikeM = repmat( tPreSpike', [length(tPostSpike) 1]);
    delta_wts = heaviside(tPostSpikeM-tPreSpikeM).*exp(-(tPostSpikeM-tPreSpikeM)/50E-3)-...
                heaviside(-(tPostSpikeM-tPreSpikeM)).*exp((tPostSpikeM-tPreSpikeM)/50E-3);
    %idZero = find(tPreSpike(spikedNeurons) == t);
    %delta_wts(idZero) =1;
    wts= wts+delta_wts;   % update weights
    
    
    
    IoutPlot = [IoutPlot Iout];
    spikedPlot = [spikedPlot PostSpike];
    TimePlot = [TimePlot t];
    N = 2;M = 2;
    subplot(N,M,1);    plot(TimePlot,IoutPlot); title('Iout vs Time');
    subplot(N,M,2);    plot(TimePlot,spikedPlot); title('Output Spikes vs Time');
    subplot(N,M,3);    surf(wts);    title('Weights');
    subplot(N,M,4);    surf(Iin);    title('InputCurrent');
    pause(0.001);
end


title('weights')