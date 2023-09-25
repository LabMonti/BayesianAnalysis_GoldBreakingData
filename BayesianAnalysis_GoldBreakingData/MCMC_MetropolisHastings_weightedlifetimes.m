function [weightedxdata] = MCMC_MetropolisHastings_weightedlifetimes(x_dag,k_0,G_dag,Nu,K,breakingspeeds,breakingspeedprobability,nsamples)
    
    %Function Description: Will sample from the probability density function
    %from  many different breaking speeds and return the array lifetimes.
    
    %notes: be careful about the number of requested samples, you dont want
    %1000's of samples from each speed if you have 1000's of speeds. more 
    %like 10-100 is reasonable.
    
    %Inputs:TraceStruct should be processed with the full attenuation
    %distribution before running through the following code. The data
    %should be logged. Breakingspeeds and their probabilities are calculated 
    %when computing the prior distribtuion. 

    %Typical values for Au junctions...
    % G_dag = 30 < G_dag < 100 pN*nm 
    % x_dag = 0.02 < x_dag < 0.30 nm 
    % k_0 = .061 1/s
    % Nu = 1/2 < Nu < 2/3 

    %Output: Sampled data.

    KbT = 4.114; %pN*nm
    originalnsamples = nsamples;
    weightedxdata = zeros(nsamples*length(breakingspeeds),1);
    loopoffset = 0;
    
    %The first thing to do is generate equally spaced
    %points along the lifetime axis of the governing
    %function. This is needed because there is a
    %cut-off where the lifetimes become imaginary which
    %depends on the parameter values chosen for
    %Nu,G_dag,x_dag, and k0.

    %These boundary values initializing tvector should 
    %always be smaller/larger than the lifetimes you 
    %are measuring experimentally. Otherwise your
    %integrations when computing fisher information
    %may be partial.
    tvector = [1.0e-4 1000];

    Nmags_tvector = order(tvector(2))-order(tvector(1));
    equalspacedtvectorArray = zeros(Nmags_tvector*10,1);
    magnitude = order(tvector(1));
    counter1 = 0;
    for i = 1:Nmags_tvector
        startvalue = 1.0*10^(magnitude+i-1);
        endvalue = 1.0*10^(magnitude+i);
        stepsize = (endvalue-startvalue)/(10-1);
        counter2 = 0;
        for l = 1:10
            equalspacedtvectorArray(l+counter1,1) = ...
                round(startvalue+(stepsize*counter2),abs(order(tvector(1)))+1);
            counter2 = counter2+1;
        end
        counter1 = counter1+10;
    end
    tArray = unique(equalspacedtvectorArray);
                    
    for i = 1:length(breakingspeeds)
        FullReqSampleArray = zeros(1,1);
        previoussamplelength = 0;
        nsamples = originalnsamples;
        [xdatarange] = DetermineSamplingStepsize_lifetimes(x_dag,k_0,G_dag,Nu,K,breakingspeeds(i),KbT,tArray);

        magnitude = order(xdatarange(2));
        xinterval = (1.0 * 10^(magnitude))/2;

        while length(FullReqSampleArray) < originalnsamples

            pdf = @(x)((k_0*(1-(((Nu*K*breakingspeeds(i)*x*x_dag)/(KbT))/(G_dag/KbT)))^((1/Nu)-1))...
            *exp((G_dag/KbT)*(1-(1-((Nu*K*breakingspeeds(i)*x*x_dag)/(KbT))/(G_dag/KbT))^(1/Nu)))*exp((k_0*KbT)/(x_dag*K*breakingspeeds(i)))...
            *exp((-(((k_0*(1-(((Nu*K*breakingspeeds(i)*x*x_dag)/(KbT))/(G_dag/KbT)))^((1/Nu)-1))...
            *exp((G_dag/KbT)*(1-(1-((Nu*K*breakingspeeds(i)*x*x_dag)/(KbT))/(G_dag/KbT))^(1/Nu))))*KbT)/(x_dag*K*breakingspeeds(i)))*...
            (1-(((Nu*K*breakingspeeds(i)*x*x_dag)/(KbT))/(G_dag/KbT)))^(1-(1/Nu))))*breakingspeedprobability(i); % Target distribution

            proppdf = @(x,y) normpdf(y-x,0,1); %proposal distribution density
            proprnd = @(x) x + rand*2*xinterval-xinterval; %random number generator for proposal dist
            smpl = mhsample(0,nsamples,'pdf',pdf,'proprnd',proprnd,'proppdf',proppdf,'symmetric',1);

            %remove any negative lifetimes
            smpl = smpl(smpl>0);
            %remove forces above the imaginary cutoff
            smpltemp = (k_0^2)-((k_0^2*Nu.*K.*breakingspeeds(i).*smpl.*x_dag)./(G_dag));
            smpl = smpl(smpltemp<.95);

            %loop over the samples and put them into a vector
            for l = 1:length(smpl)
                FullReqSampleArray(l+previoussamplelength,1) = smpl(l);
            end

            if length(FullReqSampleArray) >= originalnsamples
                FullReqSampleArray = FullReqSampleArray(1:originalnsamples);
                break
            end

            if ~isempty(smpl)
                samplelength = length(smpl)+previoussamplelength;
                previoussamplelength = samplelength;
                nsamples = originalnsamples - previoussamplelength;
            end
        end
        weightedxdata(1+loopoffset:originalnsamples*i,1) = FullReqSampleArray;
        loopoffset = originalnsamples*i;
    end
end