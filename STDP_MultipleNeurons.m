close all;
clear;
clc;
%% NOTES
% those inputs with less input spike intensity weakly modify weights !!
%% Implementation of single spiking neuron with multiple inputs 
inputNeurons = 12;
outputNeurons =5;
spikeEvents = 50;

% Neuron Spike Matrix ( Each Row represents a set of spikes)
% Iin = randi([0 1],inputNeurons,spikeEvents);      % Random Input
Iin = repmat( [ 1 1 0 0 1 1 1 0 0 0 1 1]',[1 spikeEvents]);
%Iin = 1*xor(Iin,randi([0,1],size(Iin)));
%Iin = 0.9*Iin+0.1;
Tspike = sort(randperm(300,spikeEvents))*1E-3;    % Times for input Iin spikes (must be unique)
msg = sprintf('Time-Input Spikes'); disp(msg);
for i = 1:length(Tspike)
    disp([ Tspike(i) Iin(:,i)']);
    
end
Iout = zeros(outputNeurons,1);              % Output Current as time
Wmax = 1
Wmin = 0.2
wts = (Wmax-Wmin)*rand(outputNeurons,inputNeurons)+Wmin;     % Weights matrix for input and output neuron
Tstep = 1E-3;
Tsim = 300E-3;
tau = 5E-3;
tauLTP = 5E-3
tauLTD = 10E-3;
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
    
    % calc error correction
    err  = find(Iout<0);
    Iout(err) = 0;
   
    % check for spiking accoring to threshold
    spikedNeurons = find(Iout>0.3*inputNeurons);    
    if(~isempty(spikedNeurons))
        msg = sprintf('spiked %s at t = %d\n',48+spikedNeurons,t); disp(msg);
        PostSpike(spikedNeurons) = Iout(spikedNeurons);
        Iout(spikedNeurons) = 0;    % reset the spike  neuron current to zero
        tPostSpike(spikedNeurons) = t+Tstep;     % save the spiking time
        
        
        
        
    else
        PostSpike = zeros(outputNeurons,1);
    end
    
    %Weights Update
    tPostSpikeM = repmat( tPostSpike, [1 length(tPreSpike)]);
    tPreSpikeM = repmat( tPreSpike', [length(tPostSpike) 1]);
    deltaT = tPostSpikeM-tPreSpikeM
    delta_wts = heaviside(deltaT).*exp(-(deltaT)/tauLTP)-...
                heaviside(-(deltaT)).*exp((deltaT)/tauLTD);
    %idZero = find(tPreSpike(spikedNeurons) == t);
    delta_wts(deltaT==0) = -1;
    %Add weighted contribution of current into weight change
    %delta_wts = (delta_wts.*repmat(Iinput',[outputNeurons 1]));
   
    wts= wts+0.1*delta_wts.*(Wmax-wts).*(wts-Wmin).*heaviside((wts-Wmin).*(Wmax-wts));   % update weights
%     pos = find(delta_wts>0);
%     neg = find(delta_wts<0);
%     wts(pos) = wts(pos)+0.1*delta_wts(pos).*(1-wts(pos));
%     wts(neg) = wts(neg)+0.1*delta_wts(neg).* wts(neg);
    
    
    
    IoutPlot = [IoutPlot Iout];
    spikedPlot = [spikedPlot PostSpike];
    TimePlot = [TimePlot t];
    
    
    N = 2;M = 2;
    subplot(N,M,1);    plot(TimePlot,IoutPlot); title('Iout vs Time');
    subplot(N,M,2);    plotSpike(TimePlot,spikedPlot); title('Output Spikes vs Time');
    subplot(N,M,3);    surf(wts);    title('Weights');
    subplot(N,M,4);    surf(Iin);    title('InputCurrent');
    pause(0.001);
    
end


title('weights')