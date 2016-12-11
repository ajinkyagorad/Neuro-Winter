%Initiallize Network => Network Parameters and Weight parameters
%Prepare(get) => Training/Learning Data (output of AER sensor)
%Define Simulation & relative learning paramters
%simulation(learn)
%  - check Spiking Pixels of the input Data
%  - Spiking Pixels-> Inhibit and apply STDP learning Rule
%  - No spiking ->Depretiation of Weights

clc;clear;close all;
% Time for Simulation
%Tsim  = 0.4;
Tstep = 1E-3;
% Data Generation and Network Parameters
Nout = 6;
N = 16;
M = 4; 
msg = sprintf('Size of Network = [%dx%d][%d]',M,N,Nout);
msg = sprintf('Total Neurons = %d',M*N*Nout);disp(msg);
msg = sprintf('Total Synapses = %d',M*N*2*Nout);disp(msg);
[AERout,EventTime,DataSmooth]= generateData(M,N);
EventTime = EventTime*Tstep*1;
print('AER output Data Generated');


%Simulation and learning Parameters
w_max = 2000;
w_min = 100;
w_mean = 800;
weight = w_mean+100.*randn(N*M,2,Nout);
weight((weight>w_max))= w_max; % Bound the weights 
weight((weight<w_min))= w_min;
% Learning Rate for Neurons
% LTP Rule(Long Term Potentiation)
alpha_plus =  300+ 10.*randn(N*M,2,Nout);
beta_plus = 0.5;
%LTD Rule (Long Term DepreciIation)
alpha_minus = 300+ 5.*randn(N*M,2,Nout);
beta_minus = 0.5;
% Neuron Paramters
Iout = zeros(Nout,1);
Trefrac = 10E-3;
T_LTP  = 2E-3;
Tinhibit = 1.5E-3;
tau_leak = 5E-3;
I_threshold = 0.5*N*M*2*w_mean;
msg = sprintf('Threshold Current =%d',I_threshold); disp(msg);
spikeInputPrev = zeros(Nout,1);
%Flags for conditional processing
%isLateralInhibit = zeros(Nout,1);  %use spikedNeurons Instead
TimeRefractory = zeros(Nout,1);
TimeInhibit = zeros(Nout,1);

% loop variables
tEventLast = 0; %%****
%simulation in loop
msg = sprintf('Data Duplicate'); disp(msg);
msg = sprintf('Initial Size of EventTime, %d', length(EventTime)); disp(msg);
EventTime_ = [];
for i = 1:1:5
    EventTime_ = [EventTime 0.2+EventTime_];
end
close all;
pause(0.01);
EventTime = EventTime_;
msg = sprintf('Final Size of EventTime, %d', length(EventTime)); disp(msg);
% Plots 
IoutPlot = [];
TimeInhibitPlot = [];
for i = 1:1:length(EventTime)
    
    tEvent = EventTime(i);
    AERsize = size(AERout);
    AERdata = AERout(:,:,1+mod(i,AERsize(3))); % repeat the data accordingly
    spikesP = find(AERdata==1);
    spikesN = find(AERdata==-1);
    spikes = [spikesP' spikesN'];
    i = i+1;
    msg = sprintf('Iter:%d EventTime:%d',i,tEvent);disp(msg);
    TimeInhibit(TimeInhibit>0) = TimeInhibit(TimeInhibit>0)-(tEvent-tEventLast);
    TimeInhibit(TimeInhibit<0)=0;
    TimeRefractory(TimeRefractory>0) = TimeRefractory(TimeRefractory>0)-(tEvent-tEventLast);
    TimeRefractory(TimeRefractory<0)=0;
    TimeInhibitPlot = [TimeInhibitPlot TimeInhibit];
    activeNeurons = intersect(find(TimeInhibit==0),find(TimeRefractory==0));
    % COMPUTE CURRENT
    I = Iout(activeNeurons);
    frac = exp(-(tEvent-spikeInputPrev(activeNeurons))/tau_leak);
    spikeInputPrev(activeNeurons) = tEvent;
    Iout_(activeNeurons) = I.*frac;
    Iout_(activeNeurons) = Iout_(activeNeurons)+sum(reshape(sum(weight(spikes,:,activeNeurons)),2,length(activeNeurons)));
    
    % check spiking Neurons 
    spikingNeurons = find(Iout_>I_threshold);
    Iout_(spikingNeurons) = 0;  % reset current of spikes neurons
    % update Flags
    TimeRefractory(spikingNeurons) = Trefrac;
    % update weights for spiked Neurons
    
    
     
    
   
    figure(1)
    Iout = Iout_';
    plot(Iout)
    
    IoutPlot = [IoutPlot Iout];
    pause(0.001)
end  
figure(2)
plot(EventTime,IoutPlot)
xlabel('Time(sec)')
ylabel('I(arb)')
%surf(IoutPlot)