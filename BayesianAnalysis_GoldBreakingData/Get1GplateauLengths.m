function [lifetimes,L_b,breakingspeeds] = ...
    Get1GplateauLengths(TraceStruct,PiezoSpeed,AdjustingWindow)

    %Function Description: Will compute the lifetimes, breaking lengths
    %(L_b), and breaking speeds from a trace structure based on a given
    %conductance window around 1G0, say from 0.8 to 1.2 which is a standard
    %choice.
    
    %Inputs:TraceStruct should be processed with the full attenuation
    %distribution before running through the following code. The data
    %should be logged. Piezo speed should match how the data was collected.
    %adjusting window is not generally useful, typically set to 0.

    %Output: lifetimes (s), breaking lengths/L_b (nm), breaking speeds (nm/s).

    %if different speeds are ran in the future, new # of DAQ samples will
    %have to be added. This can be found in the main program for running
    %conductance traces, inside of intialize piezo ramp.
    if PiezoSpeed <= 40
        NDAQsamples = 16384;
    end
    if PiezoSpeed == 50
        NDAQsamples = 16000;
    end
    if PiezoSpeed == 60
        NDAQsamples = 13333;
    end
    if PiezoSpeed == 80
        NDAQsamples = 10000;
    end
    if PiezoSpeed == 100
        NDAQsamples = 8000;
    end
    
    %Calculate the amount of time between each data point
    PiezoRange = 40; %microns
    ActivePiezoTime = PiezoRange/PiezoSpeed; %micron/s
    SecondsPerDAQSample = ActivePiezoTime/NDAQsamples; %second/DAQSample
        
    if AdjustingWindow == 1
        %define adjusting window intervals in linear scale
        lifetimes = zeros(TraceStruct.Ntraces,1);
        L_b = zeros(TraceStruct.Ntraces,1);
        upperWindow = 1.99:-.01:1.01;
        lowerWindow = 0.01:0.01:0.99;
        medianL_b = zeros(length(upperWindow),1);
        breakingspeeds = zeros(TraceStruct.Ntraces,1);

        for w = 1:length(upperWindow)

            for i = 1:TraceStruct.Ntraces
                trace = TraceStruct.Traces{i};
                tracecut1 = trace(trace(:,2)>log10(lowerWindow(w)),:);
                tracecut2 = tracecut1(tracecut1(:,2)<log10(upperWindow(w)),:);

                if length(tracecut2(:,1)) <= 2
                    lifetimes(i,1) = 0;
                    L_b(i,1) = 0;
                else
                    L_btosum = zeros(length(tracecut2)-1,1);
                    lifetimestosum = zeros(length(tracecut2)-1,1);
                    %get ∆d between each point in the trace
                    for l = 1:length(tracecut2)-1
                        diff = tracecut2(l+1,1) - tracecut2(l,1);
                        L_btosum(l,1) = diff;
                        lifetimestosum(l,1) = SecondsPerDAQSample;
                    end
                    lifetimes(i,1) = sum(lifetimestosum);
                    L_b(i,1) = sum(L_btosum); 
                    breakingspeeds(i,1) = L_b(i,1)/lifetimes(i,1);     
                end
            end

            %remove zeros from plateaulengths
            keep = lifetimes>0;
            lifetimes = lifetimes(keep);
            L_b = L_b(keep);
            breakingspeeds = breakingspeeds(keep);
            
            %remove negatives from breaking lengths
            keep = L_b>0;
            lifetimes = lifetimes(keep);
            L_b = L_b(keep);
            breakingspeeds = breakingspeeds(keep);
            
            %compute median values
            medianL_b(w,1) = median(L_b);

        end

        figure();
        plot(upperWindow,medianL_b);
        xlabel('Upper Conductance Window Value (G_0)')
        ylabel('median L_b (nm)')
    end
    
    %initialize empty arrays
    lifetimes = zeros(TraceStruct.Ntraces,1);
    L_b = zeros(TraceStruct.Ntraces,1);
    breakingspeeds = zeros(TraceStruct.Ntraces,1);
    
    for i = 1:TraceStruct.Ntraces
        trace = TraceStruct.Traces{i};
        tracecut1 = trace(trace(:,2)>log10(.8),:);
        tracecut2 = tracecut1(tracecut1(:,2)<log10(1.2),:);   
        
        if length(tracecut2(:,1)) <= 2
            lifetimes(i,1) = 0;
            L_b(i,1) = 0;
        else
            L_btosum = zeros(length(tracecut2)-1,1);
            lifetimestosum = zeros(length(tracecut2)-1,1);
            %get ∆d between each point in the trace
            for l = 1:length(tracecut2)-1
                diff = tracecut2(l+1,1) - tracecut2(l,1);
                L_btosum(l,1) = diff;
                lifetimestosum(l,1) = SecondsPerDAQSample;
            end
            lifetimes(i,1) = sum(lifetimestosum);
            L_b(i,1) = sum(L_btosum); 
            breakingspeeds(i,1) = L_b(i,1)/lifetimes(i,1);
        end
    end
    
    %remove zeros from plateaulengths
    keep = lifetimes>0;
    lifetimes = lifetimes(keep);
    L_b = L_b(keep);
    breakingspeeds = breakingspeeds(keep);

    %remove zeros/negatives from breaking lengths
    keep = L_b>0;
    lifetimes = lifetimes(keep);
    L_b = L_b(keep);
    breakingspeeds = breakingspeeds(keep);
   
    %histogram the length, force and lifetime values
    figure();
    histogram(log10(lifetimes));
    xlabel('Log lifetime (s)');
    ylabel('Counts');
    figure();
    histogram(log10(L_b));
    xlabel('L_b (nm)');
    ylabel('Counts');
    figure();
    histogram(log10(breakingspeeds))
    xlabel('breakingspeeds (nm/s)');
    ylabel('Counts');
end