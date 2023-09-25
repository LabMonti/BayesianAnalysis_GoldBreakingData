function [PriorStruct] = ComputePrior_weighted(TraceStruct,PiezoSpeed,Gdagvector,xdagvector,k0,Nuvector,uniformprior)
    
    %Function Description: Will compute a uniform or reference prior
    %distribution depending on user input. 
    
    %notes: be careful if changing code to include k0 as an estimator...
    %the for loops in the uniform prior, reference prior, and when
    %computing the posterior need to match.
    
    %Inputs:TraceStruct should be processed with the full attenuation
    %distribution before running through the following code. The data
    %should be logged.

    %Typical values for Au junctions...
    % G_dag = 30 < G_dag < 100 pN*nm 
    % x_dag = 0.02 < x_dag < 0.30 nm 
    % k_0 = .061 1/s
    % Nu = 1/2 < Nu < 2/3 

    %Output: Structure holding all relavent information of the prior.
    
    %set constants
    KbT = 4.114; %pN*nm
    K = 3800; %pN/nm
    PiezoRange = 40; %microns
    
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
    ActivePiezoTime = PiezoRange/PiezoSpeed;
    SecondsPerDAQSample = ActivePiezoTime/NDAQsamples;
    lifetimes = zeros(TraceStruct.Ntraces,1);
    L_b = zeros(TraceStruct.Ntraces,1);
    breakingspeeds = zeros(TraceStruct.Ntraces,1);
    
    for i = 1:TraceStruct.Ntraces
        
        %draw a conductance window
        trace = TraceStruct.Traces{i};
        tracecut1 = trace(trace(:,2)>log10(0.8),:);
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
    
    %remove any zeros from the lifetimes and breakingspeeds
    keep = lifetimes>0;
    lifetimes = lifetimes(keep);
    L_b = L_b(keep);
    breakingspeeds = breakingspeeds(keep);
    keep = breakingspeeds>0;
    lifetimes = lifetimes(keep);
    L_b = L_b(keep);
    breakingspeeds = breakingspeeds(keep);  
    
    %get probabilities of breakingspeeds
    [~,breakingspeedsindex] = sort(breakingspeeds);
    h = histogram(log10(breakingspeeds));
    total_loopcounter = 0;
    sortedbreakingspeedprobability = zeros(length(breakingspeeds),1);
    for i = 1:length(h.Values)
        for l = 1:h.Values(i)
            total_loopcounter = total_loopcounter+1;
            sortedbreakingspeedprobability(total_loopcounter,1) = h.Values(i)/sum(h.Values);
        end
    end
    
    %reorder so lifetimes,breakingspeeds and their probabilities match
    breakingspeedprobability = zeros(length(breakingspeeds),1);
    for i = 1:length(breakingspeeds)
        breakingspeedprobability(breakingspeedsindex(i)) = sortedbreakingspeedprobability(i);
    end
    
    %allocate empty matrices and vectors
    total_loopcount = 0;
    priorvector = zeros(length(Nuvector)*length(Gdagvector)*length(xdagvector),1);
    xvector = zeros(length(Nuvector)*length(Gdagvector)*length(xdagvector),1);
    yvector = zeros(length(Nuvector)*length(Gdagvector)*length(xdagvector),1);
    zvector = zeros(length(Nuvector)*length(Gdagvector)*length(xdagvector),1); 
      
%~~~~Uniform prior code~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if uniformprior == 1 
        for n = 1:length(Nuvector)
            for g = 1:length(Gdagvector)
                for x = 1:length(xdagvector)
                    total_loopcount = total_loopcount+1;
                    xvector(total_loopcount,1) = xdagvector(x);
                    yvector(total_loopcount,1) = Gdagvector(g);
                    zvector(total_loopcount,1) = Nuvector(n);
                    priorvector(total_loopcount,1) = 1;
                end
            end
        end
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    else

%~~~~Reference Prior Code~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        %Set stepsizes
        Nustepsize = Nuvector(2) - Nuvector(1);
        Gdagstepsize = Gdagvector(2) - Gdagvector(1);
        
        %set the triple for loop for each combination
        for n = 1:length(Nuvector)
            for g = 1:length(Gdagvector)
                for x = 1:length(xdagvector)
                    
                    %set the indices for the about to be calculated reference
                    %prior point
                    total_loopcount = total_loopcount+1;
                    xvector(total_loopcount,1) = xdagvector(x);
                    yvector(total_loopcount,1) = Gdagvector(g);
                    zvector(total_loopcount,1) = Nuvector(n);

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
                    tvector = unique(equalspacedtvectorArray);
                    finalmaxtVector = 0;

                    %compute largest value of xdata before going imaginary
                    for i = 1:length(breakingspeeds)
                        tArraytemp = (k0^2)-((k0^2*...
                            Nuvector(n).*K.*breakingspeeds(i).*tvector.*xdagvector(x))./(Gdagvector(g)));
                        tVector_cut = tvector(tArraytemp > 0);
                        if max(tVector_cut) >= finalmaxtVector
                            finalmaxtVector = max(tVector_cut);
                        end
                    end

                    %Now that we know the largest lifetime value given
                    %the specific Nu, G_dag, x_dag, and k0 values, we
                    %need to adjust the largest lifetime of our synthetic 
                    %data to match exactly.
                    adjustedvalue = finalmaxtVector;
                    tomultiply = order(finalmaxtVector);
                    if tomultiply < 0
                        for i = 1:-1*tomultiply
                            adjustedvalue = adjustedvalue*10;
                        end
                    else
                        for i = 1:tomultiply
                            adjustedvalue = adjustedvalue/10;
                        end
                    end

                    %now generate the synthetic data which will be
                    %integrating over when we calculate our fisher
                    %information
                    Npointsperdecade_xdata = 100;
                    Nmags = order(finalmaxtVector)-order(min(tvector));
                    equalspacedVector = zeros(Nmags*Npointsperdecade_xdata,1);
                    magnitude = order(min(tvector));
                    counter1 = 0;
                    for i = 1:Nmags
                        startvalue = 1.0*10^(magnitude+i-1);
                        if i == Nmags
                            endvalue = adjustedvalue*10^(magnitude+i);
                            stepsize = (endvalue-startvalue)/(Npointsperdecade_xdata-1);
                            counter2 = 0;
                            for l = 1:Npointsperdecade_xdata
                                equalspacedVector(l+counter1,1) = ...
                                    round(startvalue+(stepsize*counter2),abs(order(min(tvector)))+1);
                                counter2 = counter2+1;
                            end
                            counter1 = counter1+Npointsperdecade_xdata;
                        else
                            endvalue = 1.0*10^(magnitude+i);
                            stepsize = (endvalue-startvalue)/(Npointsperdecade_xdata-1);
                            counter2 = 0;
                            for l = 1:Npointsperdecade_xdata
                                equalspacedVector(l+counter1,1) = ...
                                    round(startvalue+(stepsize*counter2),abs(order(min(tvector)))+1);
                                counter2 = counter2+1;
                            end
                            counter1 = counter1+Npointsperdecade_xdata;
                        end

                    end
                    xdata = unique(equalspacedVector);
                    xdatastepsizeVector = zeros(length(xdata),1);
                    for i = 2:length(xdata)
                        xdatastepsizeVector(i) = xdata(i) - xdata(i-1);
                    end                       

%compute the fisher information for Nu for the triplet prior point~~~~~~~~~

                    finaltermstosum = zeros(length(xdata),1);
                    for i = 1:length(xdata)

                        %compute the inner sum of the derivative terms
                        innerderivativetermstosum = zeros(length(breakingspeeds),1);
                        inneroriginaltermstosum = zeros(length(breakingspeeds),1);
                        for l = 1:length(breakingspeeds)

                            %compute the log partial derivative w.r.t. Nu
                            %of our main function
                            term1 = -(1/(Nuvector(n)^2))-(1/Nuvector(n));
                            term2 = -(Gdagvector(g)/KbT)*(-((((K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))...
                                    /(Nuvector(n)*(((-(Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))+1)))...
                                    -(log((((-(Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))+1))/Nuvector(n)^2))*...
                                    (((-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))+1))^(1/Nuvector(n));
                            term3 = (k0*exp((Gdagvector(g)/KbT)*(-(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))...
                                    /(Gdagvector(g)/KbT)))^(1/Nuvector(n))+1))*(Gdagvector(g)/KbT)*KbT*(log(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))...
                                    /(Gdagvector(g)/KbT)))-((((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))...
                                    /(Gdagvector(g)/KbT))*log(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))...
                                    /(Gdagvector(g)/KbT))))+(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))...
                                    /(Gdagvector(g)/KbT)))*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))...
                                    /(Gdagvector(g)/KbT)))^((1-Nuvector(n))/Nuvector(n)))/(Nuvector(n)^2*xdagvector(x)*K*breakingspeeds(l));

                            finalterm = term1+term2-term3;

                            innerderivativetermstosum(l,1) = finalterm;

                            %compute the main function value
                            originalfuncvalue = ((k0*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^((1/Nuvector(n))-1))...
                            *exp((Gdagvector(g)/KbT)*(1-(1-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))^(1/Nuvector(n))))*exp((k0*KbT)/(xdagvector(x)*K*breakingspeeds(l)))...
                            *exp((-(((k0*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^((1/Nuvector(n))-1))...
                            *exp((Gdagvector(g)/KbT)*(1-(1-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))^(1/Nuvector(n)))))*KbT)/(xdagvector(x)*K*breakingspeeds(l)))*...
                            (1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^(1-(1/Nuvector(n)))))*breakingspeedprobability(l);

                            inneroriginaltermstosum(l,1) = originalfuncvalue;
                        end

                        %remove any speeds that may produce imaginary
                        %or negative numbers
                        keep = imag(innerderivativetermstosum)==0;
                        innerderivativetermstosum = innerderivativetermstosum(keep);
                        inneroriginaltermstosum = inneroriginaltermstosum(keep);
                        keep = imag(inneroriginaltermstosum)==0;
                        innerderivativetermstosum = innerderivativetermstosum(keep);
                        inneroriginaltermstosum = inneroriginaltermstosum(keep);
                        keep = inneroriginaltermstosum>0;
                        innerderivativetermstosum = innerderivativetermstosum(keep);
                        inneroriginaltermstosum = inneroriginaltermstosum(keep);

                        summedderivativeterm = sum(innerderivativetermstosum,'omitnan');                    
                        summedoriginalterms = sum(inneroriginaltermstosum,'omitnan');

                        finaltermstosum(i,1) = summedderivativeterm^2*summedoriginalterms*xdatastepsizeVector(i);
                    end                        

                    FisherInformation = sum(finaltermstosum,'omitnan');

%now compute the conditional prior point for (Nu|Gdag)~~~~~~~~~~~~~~~~~~~~~                        

                    %compute the integral of the exp{(1/2)*log(fisherinformation(theta))}
                    termstosum = zeros(n,1);
                    for l = 1:n
                        termstosum(l,1) = exp((1/2)*log(FisherInformation))*Nustepsize;
                    end

                    c_Nu = sum(termstosum,'omitnan');

                    Nu_Gdag_prior = (exp((1/2)*log(FisherInformation))/c_Nu);

%compute the fisher information for Gdag for the triplet prior point~~~~~~~                      

                    %compute the integral via reimann sum
                    finaltermstosum = zeros(length(xdata),1);
                    for i = 1:length(xdata)

                        %compute the inner sum of the derivative terms
                        innerderivativetermstosum = zeros(length(breakingspeeds),1);
                        inneroriginaltermstosum = zeros(length(breakingspeeds),1);
                        for l = 1:length(breakingspeeds)

                            %compute the log partial derivative
                            term1 = (1-(1-((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*xdagvector(x)))^(1/Nuvector(n)))/KbT;
                            term2 = -((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)*xdagvector(x))*(1-((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*xdagvector(x)))^((1/Nuvector(n))-1))...
                                    /((Gdagvector(g)/KbT)*KbT*Nuvector(n));
                            term3 = term1+term2;
                            term4 = ((k0*KbT*exp((Gdagvector(g)/KbT)*(1-(1-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))^(1/Nuvector(n)))))*term3)...
                                    /(K*breakingspeeds(l)*xdagvector(x));

                            finalterm = term3-term4;

                            innerderivativetermstosum(l,1) = finalterm;

                            %compute the orginal function value
                            originalfuncvalue = ((k0*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^((1/Nuvector(n))-1))...
                            *exp((Gdagvector(g)/KbT)*(1-(1-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))^(1/Nuvector(n))))*exp((k0*KbT)/(xdagvector(x)*K*breakingspeeds(l)))...
                            *exp((-(((k0*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^((1/Nuvector(n))-1))...
                            *exp((Gdagvector(g)/KbT)*(1-(1-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))^(1/Nuvector(n)))))*KbT)/(xdagvector(x)*K*breakingspeeds(l)))*...
                            (1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^(1-(1/Nuvector(n)))))*breakingspeedprobability(l);

                            inneroriginaltermstosum(l,1) = originalfuncvalue;
                        end

                        %remove any speeds that may produce imaginary
                        %or negative numbers
                        keep = imag(innerderivativetermstosum)==0;
                        innerderivativetermstosum = innerderivativetermstosum(keep);
                        inneroriginaltermstosum = inneroriginaltermstosum(keep);
                        keep = imag(inneroriginaltermstosum)==0;
                        innerderivativetermstosum = innerderivativetermstosum(keep);
                        inneroriginaltermstosum = inneroriginaltermstosum(keep);
                        keep = inneroriginaltermstosum>0;
                        innerderivativetermstosum = innerderivativetermstosum(keep);
                        inneroriginaltermstosum = inneroriginaltermstosum(keep);

                        summedderivativeterm = sum(innerderivativetermstosum,'omitnan');                    
                        summedoriginalterms = sum(inneroriginaltermstosum,'omitnan');

                        finaltermstosum(i,1) = summedderivativeterm^2*summedoriginalterms*xdatastepsizeVector(i);
                    end

                    FisherInformation = sum(finaltermstosum,'omitnan');

%now compute the conditional prior point for (Gdag|xdag)~~~~~~~~~~~~~~~~~~~

                    %compute the integral of the exp{(1/2)*log(fisherinformation(theta))}
                    termstosum = zeros(g,1);
                    for l = 1:g
                        termstosum(l,1) = exp((1/2)*log(FisherInformation))*Gdagstepsize;
                    end

                    c_Gdag = sum(termstosum,'omitnan');

                    Gdag_xdag_prior = (exp((1/2)*log(FisherInformation))/c_Gdag);

%now compute the fisher information for xdag for the triplet prior point~~~

                    %compute the integral via reimann sum
                    finaltermstosum = zeros(length(xdata),1);
                    for i = 1:length(xdata)

                        %compute the inner sum of the derivative terms
                        innerderivativetermstosum = zeros(length(breakingspeeds),1);
                        inneroriginaltermstosum = zeros(length(breakingspeeds),1);
                        for l = 1:length(breakingspeeds)

                            %compute the log partial derivative
                            term1 = -((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*(Nuvector(n)-1))/(Nuvector(n)*(((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*xdagvector(x))-1));
                            term2 = ((Gdagvector(g)/KbT)*((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*xdagvector(x))^((1/Nuvector(n))-1)))/...
                                    Nuvector(n);
                            term3 = -(KbT*k0)/(K*breakingspeeds(l)*xdagvector(x)^2);
                            term4 = -(KbT*k0*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*xdagvector(x))^(2-(2/Nuvector(n)))...
                                    *exp((Gdagvector(g)/KbT)*(1-(1-((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*xdagvector(x)))^(1/Nuvector(n)))))/(K*breakingspeeds(l)*xdagvector(x)^2);
                            term5 = -((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*KbT*k0*(2-(2/Nuvector(n)))*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*xdagvector(x))^(1-(2/Nuvector(n)))...
                                    *exp((Gdagvector(g)/KbT)*(1-(1-((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*xdagvector(x)))^(1/Nuvector(n)))))/(K*breakingspeeds(l)*xdagvector(x));
                            term6 = ((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*(Gdagvector(g)/KbT)*KbT*k0*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*xdagvector(x))^(1-(1/Nuvector(n)))...
                                    *exp((Gdagvector(g)/KbT)*(1-(1-((((Nuvector(n)*K*breakingspeeds(l)*xdata(i))/KbT)/(Gdagvector(g)/KbT))*xdagvector(x)))^(1/Nuvector(n)))))/(Nuvector(n)*K*breakingspeeds(l)*xdagvector(x));

                            finalterm = term1+term2+term3-term4-term5-term6;

                            innerderivativetermstosum(l) = finalterm;

                            %compute the orginal function value
                            originalfuncvalue = ((k0*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^((1/Nuvector(n))-1))...
                            *exp((Gdagvector(g)/KbT)*(1-(1-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))^(1/Nuvector(n))))*exp((k0*KbT)/(xdagvector(x)*K*breakingspeeds(l)))...
                            *exp((-(((k0*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^((1/Nuvector(n))-1))...
                            *exp((Gdagvector(g)/KbT)*(1-(1-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))^(1/Nuvector(n)))))*KbT)/(xdagvector(x)*K*breakingspeeds(l)))*...
                            (1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^(1-(1/Nuvector(n)))))*breakingspeedprobability(l);

                            inneroriginaltermstosum(l,1) = originalfuncvalue;
                        end

                        %remove any speeds that may produce imaginary
                        %or negative numbers
                        keep = imag(innerderivativetermstosum)==0;
                        innerderivativetermstosum = innerderivativetermstosum(keep);
                        inneroriginaltermstosum = inneroriginaltermstosum(keep);
                        keep = imag(inneroriginaltermstosum)==0;
                        innerderivativetermstosum = innerderivativetermstosum(keep);
                        inneroriginaltermstosum = inneroriginaltermstosum(keep);
                        keep = inneroriginaltermstosum>0;
                        innerderivativetermstosum = innerderivativetermstosum(keep);
                        inneroriginaltermstosum = inneroriginaltermstosum(keep);

                        summedderivativeterm = sum(innerderivativetermstosum,'omitnan');
                        summedoriginalterms = sum(inneroriginaltermstosum,'omitnan');
                        
                        finaltermstosum(i,1) = summedderivativeterm^2*summedoriginalterms*xdatastepsizeVector(i);
                    end

                    FisherInformation = sum(finaltermstosum,'omitnan');

%now compute the final joint reference prior point~~~~~~~~~~~~~~~~~~~~~~~~~

                    %integrate out the first two parameters
                    termstosum_n = zeros(n,1);
                    termstosum_g = zeros(g,1);  
                    for q = 1:n
                        for s = 1:g
                            termstosum_g(s) =  Nu_Gdag_prior*Gdag_xdag_prior*...
                                    log(FisherInformation)*Gdagstepsize;
                        end
                        termstosum_n(q) = sum(termstosum_g,'omitnan')*Nustepsize;
                    end

                    %compute the x_dag prior point
                    xdag_prior = exp(-(1/2)*sum(termstosum_n,'omitnan'));  

                    %compute the final rference prior point
                    priorvector(total_loopcount,1) = Nu_Gdag_prior*Gdag_xdag_prior*...
                        xdag_prior;
                end
            end
        end
    end

    %plot the reference prior
    figure();
    scatter3(xvector,yvector,zvector,80,priorvector,'filled')
    xlabel('$x^\ddagger\:(nm)$','Interpreter','latex','FontSize',18)
    ylabel('$\Delta G^\ddagger\:(pN nm)$','Interpreter','latex','FontSize',18)
    zlabel('\nu','FontSize',18)
    title('Prior Distribution','FontSize',18)
    colormap(copper)
    colorbar
    a = colorbar;
    a.Label.Interpreter = 'latex';
    a.Label.String = '$p(\:x^\ddagger,\Delta G^\ddagger,\nu)$';
    a.FontSize = 18;
    
    %always make outputs row vectors
    if length(Nuvector(:,1))==1
        Nuvector = Nuvector';
    end
    if length(Gdagvector(:,1))==1
        Gdagvector = Gdagvector';
    end
    if length(xdagvector(:,1))==1
        xdagvector = xdagvector';
    end
    
    %assign structure fields
    PriorStruct.prior = priorvector;
    PriorStruct.Nuvector = Nuvector;
    PriorStruct.Gdagvector = Gdagvector;
    PriorStruct.xdagvector = xdagvector;
    PriorStruct.k0vector = k0;
    PriorStruct.xvector = xvector;
    PriorStruct.yvector = yvector;
    PriorStruct.zvector = zvector;
    PriorStruct.K = K;
    PriorStruct.lifetimes = lifetimes;
    PriorStruct.L_b = L_b;
    PriorStruct.breakingspeeds = breakingspeeds;
    PriorStruct.breakingspeedprobability = breakingspeedprobability;
    PriorStruct.KbT = KbT; 
    
%spare code~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     %generate equally spaced points for Gdag if varying over more than 1
%     %order of magnitude
%     if length(Gdagvector) > 1
%         Gdagstepsize = Gdagvector(2)-Gdagvector(1);
%         if order(Gdagvector(1)) < order(Gdagvector(2)) && length(Gdagvector) == 2
%             Nmags_g = order(Gdagvector(2))-order(Gdagvector(1));
%             equalspacedGdagvectorArray = zeros(Nmags_g*Npointsperdecade,1);
%             magnitude = order(Gdagvector(1));
%             counter1 = 0;
% 
%             for i = 1:Nmags_g
%                 startvalue = 1.0*10^(magnitude+i-1);
%                 endvalue = 1.0*10^(magnitude+i);
%                 stepsize = (endvalue-startvalue)/(Npointsperdecade-1);
%                 counter2 = 0;
%                 for l = 1:Npointsperdecade
%                     equalspacedGdagvectorArray(l+counter1,1) = ...
%                         round(startvalue+(stepsize*counter2),abs(order(Gdagvector(1)))+1);
%                     counter2 = counter2+1;
%                 end
%                 counter1 = counter1+Npointsperdecade;
%             end
%             Gdagvector = unique(equalspacedGdagvectorArray);
%         end
%     end
%     
%     %generate equally spaced points for k0 if varying over more than 1
%     %order of magnitude
%     if length(k0vector) > 1
%         if order(k0vector(1)) < order(k0vector(2)) && length(k0vector) == 2
%             Nmags_k = order(k0vector(2))-order(k0vector(1));
%             equalspacedk0Array = zeros(Nmags_k*Npointsperdecade,1);
%             magnitude = order(k0vector(1));
%             counter1 = 0;
% 
%             for i = 1:Nmags_k
%                 startvalue = 1.0*10^(magnitude+i-1);
%                 endvalue = 1.0*10^(magnitude+i);
%                 stepsize = (endvalue-startvalue)/(Npointsperdecade-1);
%                 counter2 = 0;
%                 for l = 1:Npointsperdecade
%                     equalspacedk0Array(l+counter1,1) = ...
%                         round(startvalue+(stepsize*counter2),abs(order(k0vector(1)))+1);
%                     counter2 = counter2+1;
%                 end
%                 counter1 = counter1+Npointsperdecade;
%             end
%             k0vector = unique(equalspacedk0Array);
%         end
%     end
%     
%     %generate equally spaced points for xdag if varying over more than 1
%     %order of magnitude
%     if length(xdagvector) > 1
%         if order(xdagvector(1)) < order(xdagvector(2)) && length(xdagvector) == 2
%             multimags = 1;
%             Nmags_x = order(xdagvector(2))-order(xdagvector(1));
%             equalspacedxdagArray = zeros(Nmags_x*Npointsperdecade,1);
%             magnitude = order(xdagvector(1));
%             xdagstepsizeArray = zeros(Nmags_x,1);
%             counter1 = 0;
% 
%             for i = 1:Nmags_x
%                 startvalue = 1.0*10^(magnitude+i-1);
%                 endvalue = 1.0*10^(magnitude+i);
%                 stepsize = (endvalue-startvalue)/(Npointsperdecade-1);
%                 counter2 = 0;
%                 for l = 1:Npointsperdecade
%                     equalspacedxdagArray(l+counter1,1) = ...
%                         round(startvalue+(stepsize*counter2),abs(order(xdagvector(1)))+1);
%                     counter2 = counter2+1;
%                 end
%                 counter1 = counter1+Npointsperdecade;
%                 xdagstepsizeArray(i) = stepsize;
%             end
%             xdagvector = unique(equalspacedxdagArray);
% 
%             %set up arrays for keepig track of which stepsize to use
%             rangeArray_x = zeros(Nmags_x,1);
%             rangeArray_x(1) = Npointsperdecade;
%             for i = 1:Nmags_x-1
%                 rangeArray_x(i+1) = Npointsperdecade+((Npointsperdecade-1)*i);
%             end
%         else
%             multimags = 0;
%             xdagstepsize = xdagvector(2) - xdagvector(1); 
%         end
%     end
%     
%                         %compute the integral of the exp{(1/2)*log(fisherinformation(theta))}
%                         termstosum = zeros(x,1);
%                         rangecounter_x = 1;
%                         for l = 1:x
%                             if multimags == 1
%                                  if l < rangeArray_x(rangecounter_x)
%                                      xdagstepsize = xdagstepsizeArray(rangecounter_x);
%                                  elseif l > rangeArray_x(rangecounter_x)
%                                      xdagstepsize = xdagstepsizeArray(rangecounter_x+1);
%                                      rangecounter_x = rangecounter_x+1;
%                                  end
%                             end
%                             termstosum(l,1) = exp((1/2)*log(FisherInformation))*xdagstepsize;
%                         end
% 
%                         c_xdag = sum(termstosum,'omitnan');
% 
%                         %compute the prior for the given theta value and store it
%                         xdag_k0_prior = (exp((1/2)*log(FisherInformation))/c_xdag);

%                         finaltermstosum = zeros(length(xdata),1);
%                         for i = 1:length(xdata)
% 
%                             %compute the inner sum of the derivative terms
%                             innerderivativetermstosum = zeros(length(breakingspeeds),1);
%                             inneroriginaltermstosum = zeros(length(breakingspeeds),1);
%                             for l = 1:length(breakingspeeds)
% 
%                                 %compute the log partial derivative
%                                 term1 = 1/(k0vector(k));
%                                 term2 = KbT/(xdagvector(x)*K*breakingspeeds(l));
%                                 term3 = (-(((k0vector(k)*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^((1/Nuvector(n))-1))...
%                                 *exp((Gdagvector(g)/KbT)*(1-(1-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))^(1/Nuvector(n)))))*KbT)/(xdagvector(x)*K*breakingspeeds(l)))*...
%                                 (1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^(1-(1/Nuvector(n)));
% 
%                                 finalterm = term1+term2+term3;
% 
%                                 innerderivativetermstosum(l) = finalterm;
% 
%                                 %compute the orginal function value
%                                 originalfuncvalue = ((k0vector(k)*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^((1/Nuvector(n))-1))...
%                                 *exp((Gdagvector(g)/KbT)*(1-(1-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))^(1/Nuvector(n))))*exp((k0vector(k)*KbT)/(xdagvector(x)*K*breakingspeeds(l)))...
%                                 *exp((-(((k0vector(k)*(1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^((1/Nuvector(n))-1))...
%                                 *exp((Gdagvector(g)/KbT)*(1-(1-((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT))^(1/Nuvector(n)))))*KbT)/(xdagvector(x)*K*breakingspeeds(l)))*...
%                                 (1-(((Nuvector(n)*K*breakingspeeds(l)*xdata(i)*xdagvector(x))/(KbT))/(Gdagvector(g)/KbT)))^(1-(1/Nuvector(n)))))*breakingspeedprobability(l);
% 
%                                 inneroriginaltermstosum(l,1) = originalfuncvalue;
%                             end
% 
%                             %remove any speeds that may produce imaginary
%                             %numbers
%                             keep = imag(innerderivativetermstosum)==0;
%                             innerderivativetermstosum = innerderivativetermstosum(keep);
%                             inneroriginaltermstosum = inneroriginaltermstosum(keep);
%                             keep = imag(inneroriginaltermstosum)==0;
%                             innerderivativetermstosum = innerderivativetermstosum(keep);
%                             inneroriginaltermstosum = inneroriginaltermstosum(keep);
% 
%                             summedderivativeterm = sum(innerderivativetermstosum,'omitnan');
%                             summedoriginalterms = sum(inneroriginaltermstosum,'omitnan');
% 
%                             %compute the outer sum 
%                             if isreal(summedderivativeterm^2*summedoriginalterms*xdatastepsizeVector(i)) 
%                                 finaltermstosum(i,1) = summedderivativeterm^2*summedoriginalterms*xdatastepsizeVector(i);
%                             else
%                                 finaltermstosum(i,1) = 0;
%                             end
%                         end
% 
%                         FisherInformation = sum(finaltermstosum,'omitnan');
%spare code~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end