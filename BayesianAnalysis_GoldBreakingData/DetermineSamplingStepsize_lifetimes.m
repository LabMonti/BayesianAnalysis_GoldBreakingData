function [xdatarange] = DetermineSamplingStepsize_lifetimes(x_dag,k_0,G_dag,pot,K,V,KbT,tArray)
    
    %Function Description: Can be used to return the data range given some 
    %set of paramaters or to  visualize the probability density function
    %given some set of parameters.
    
    %notes: be careful if plotting, when used inside of
    %MCMC_Metropolis_hastings function, many plots could potentially be
    %produced.
    
    %Inputs:

    %Typical values for Au junctions...
    % G_dag = 30 < G_dag < 100 pN*nm 
    % x_dag = 0.02 < x_dag < 0.30 nm 
    % k_0 = .061 1/s
    % Nu = 1/2 < Nu < 2/3 
    %tArray can be [1.0e-4:stepsize:1000] but needs to include the lifetimes you
    %would expect from the parameter choices and stepsize should be
    %resonable (10 points per decade).
    
    %Output: xdata range of lifetimes given some parameter values.

    %check to make sure we stay in the real data range only
    tArraytemp = (k_0^2)-((k_0^2*pot.*K.*V.*tArray.*x_dag)./(G_dag));
    tArray = tArray(tArraytemp > 0);
    
    %set the stepsize
    magnitude = order(max(tArray));
    stepsize = 1.0*10^(magnitude-3);
    xdata = min(tArray):stepsize:max(tArray);
    ydata = zeros(length(xdata),1);

    for i = 1:length(xdata)
            ydata(i,1) = (k_0*(1-(((pot*K*V*xdata(i)*x_dag)/(KbT))/(G_dag/KbT)))^((1/pot)-1))...
                *exp((G_dag/KbT)*(1-(1-((pot*K*V*xdata(i)*x_dag)/(KbT))/(G_dag/KbT))^(1/pot)))*exp((k_0*KbT)/(x_dag*K*V))...
                *exp((-(((k_0*(1-(((pot*K*V*xdata(i)*x_dag)/(KbT))/(G_dag/KbT)))^((1/pot)-1))...
                *exp((G_dag/KbT)*(1-(1-((pot*K*V*xdata(i)*x_dag)/(KbT))/(G_dag/KbT))^(1/pot))))*KbT)/(x_dag*K*V))*...
                (1-(((pot*K*V*xdata(i)*x_dag)/(KbT))/(G_dag/KbT)))^(1-(1/pot)));
    end
    
    %cut all values before the first zero going from left to right
    leftindex = 1;
    for i = 2:length(ydata)
        if real(ydata(i)) == 0
            leftindex = i;
        elseif real(ydata(i-1)) < real(ydata(i))
            leftindex = 1;
            break
        end
    end
    
    xdata = xdata(leftindex:length(ydata));
    ydata = ydata(leftindex:length(ydata));
    
    %cut all the values after the last zero going right ot left
    keep = real(ydata) == 0;
    
    if sum(keep) == 0
        rightindex = length(ydata);
    else
        for i = 1:length(ydata)
            if keep(i) == 1
                rightindex = i;
                break
            end 
        end
    end
    
    xdata = xdata(leftindex:rightindex);
%     ydata = real(ydata(leftindex:rightindex));
%     hold on
%     scatter(xdata,ydata,'filled')
%     xlabel('t (s)')
%     ylabel('p(t|V)')
    xdatarange(1,1) = min(xdata);
    xdatarange(1,2) = max(xdata);
end