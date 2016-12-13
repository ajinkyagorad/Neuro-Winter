function[] = plotSpike(time,spikes)
    [Nout,len] = size(spikes);
    spikes(find(spikes)) = (1+mod(find(spikes),Nout))/5;
    plot(time,spikes,'.');
end
