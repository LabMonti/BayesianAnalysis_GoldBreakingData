function DetermineSamplingStepsize_lifetimes_weighted(breakingspeeds,breakingspeedprobability,x_dag,k_0,G_dag,Nu,K)

    %Function Description: Will plot a weighted probability density function 
    %of lifetimes. 
    
    %Inputs:Breakingspeeds and their probabilities are calculated when
    %computing the prior distribtuion. 

    %Typical values for Au junctions...
    % G_dag = 30 < G_dag < 100 pN*nm 
    % x_dag = 0.02 < x_dag < 0.30 nm 
    % k_0 = .061 1/s
    % Nu = 1/2 < Nu < 2/3
    % K = 3800 pN/nm

    %Output: Plot of the weighted lifetime probability density function
    %given some parameter choices.

    KbT = 4.114; %pN/nm
    
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
        tArraytemp = (k_0^2)-((k_0^2*...
            Nu.*K.*breakingspeeds(i).*tvector.*x_dag)./(G_dag));
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

    %now generate the synthetic data
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
    
    %compute the ydata points at each breaking speed then sum them
    ydatareal = zeros(length(xdata),1);
    ydataArray = zeros(length(xdata),length(breakingspeeds));
    for l = 1:length(breakingspeeds)
        for i = 1:length(xdata)
            ydata = ((k_0*(1-(((Nu*K*breakingspeeds(l)*xdata(i)*x_dag)/(KbT))/(G_dag/KbT)))^((1/Nu)-1))...
                *exp((G_dag/KbT)*(1-(1-((Nu*K*breakingspeeds(l)*xdata(i)*x_dag)/(KbT))/(G_dag/KbT))^(1/Nu)))*exp((k_0*KbT)/(x_dag*K*breakingspeeds(l)))...
                *exp((-(((k_0*(1-(((Nu*K*breakingspeeds(l)*xdata(i)*x_dag)/(KbT))/(G_dag/KbT)))^((1/Nu)-1))...
                *exp((G_dag/KbT)*(1-(1-((Nu*K*breakingspeeds(l)*xdata(i)*x_dag)/(KbT))/(G_dag/KbT))^(1/Nu))))*KbT)/(x_dag*K*breakingspeeds(l)))*...
                (1-(((Nu*K*breakingspeeds(l)*xdata(i)*x_dag)/(KbT))/(G_dag/KbT)))^(1-(1/Nu))))*breakingspeedprobability(l);

            if isreal(ydata) && ydata>0
                ydatareal(i,1) = ydata;
            else
                ydatareal(i,1) = 0;
            end
        end
        ydataArray(:,l) = ydatareal; 
    end
    
    summedydata = sum(ydataArray,2,'omitnan');
    normalization = sum(summedydata,'omitnan');
    data = summedydata./normalization;%renormalize everything
    hold on
    plot(log10(xdata),data,'linewidth',3)
    xlabel('time (s)')
    ylabel('p(t)')
end