function [PosteriorStruct] = ComputePosterior_weighted(PriorStruct)

    %Function Description: Will compute a posterior and 95% credible interval 
    %distribution. 
    
    %notes: be careful if changing code to include k0 as an estimator...
    %the for loops in the uniform prior, reference prior, and when
    %computing the posterior need to match.
    
    %Inputs:The PriorStruct should come from ComputePrior_weighted.

    %Typical values for Au junctions...
    % G_dag = 30 < G_dag < 100 pN*nm 
    % x_dag = 0.02 < x_dag < 0.30 nm 
    % k_0 = .061 1/s
    % Nu = 1/2 < Nu < 2/3 

    %Output: Structure holding all relavent information of the posterior.

    %allocate empty matrices and vectors
    tempposteriorvector = zeros(length(PriorStruct.prior),1);

    for i = 1:length(PriorStruct.lifetimes)
        total_loopcount = 0;
        for n = 1:length(PriorStruct.Nuvector)
            for g = 1:length(PriorStruct.Gdagvector)
                for x = 1:length(PriorStruct.xdagvector)
                    
                    %compute the posterior
                    total_loopcount = total_loopcount+1;

                    likelihood = ((PriorStruct.k0vector.*...
                        (1-(((PriorStruct.Nuvector(n).*PriorStruct.K.*PriorStruct.breakingspeeds*PriorStruct.lifetimes(i).*PriorStruct.xdagvector(x))./(PriorStruct.KbT))./(PriorStruct.Gdagvector(g)./PriorStruct.KbT)))...
                        .^((1./PriorStruct.Nuvector(n))-1)).*exp((PriorStruct.Gdagvector(g)./PriorStruct.KbT).*(1-(1-((PriorStruct.Nuvector(n).*PriorStruct.K.*PriorStruct.breakingspeeds.*PriorStruct.lifetimes(i).*...
                        PriorStruct.xdagvector(x))./(PriorStruct.KbT))./(PriorStruct.Gdagvector(g)./PriorStruct.KbT)).^(1/PriorStruct.Nuvector(n))))...
                        .*exp((PriorStruct.k0vector*PriorStruct.KbT)./(PriorStruct.xdagvector(x).*PriorStruct.K.*PriorStruct.breakingspeeds))...
                        .*exp((-(((PriorStruct.k0vector.*(1-(((PriorStruct.Nuvector(n).*PriorStruct.K.*PriorStruct.breakingspeeds.*PriorStruct.lifetimes(i).*PriorStruct.xdagvector(x))...
                        ./(PriorStruct.KbT))./(PriorStruct.Gdagvector(g)./PriorStruct.KbT))).^((1./PriorStruct.Nuvector(n))-1))...
                        .*exp((PriorStruct.Gdagvector(g)./PriorStruct.KbT).*(1-(1-((PriorStruct.Nuvector(n).*PriorStruct.K.*PriorStruct.breakingspeeds.*PriorStruct.lifetimes(i).*PriorStruct.xdagvector(x))...
                        ./(PriorStruct.KbT))./(PriorStruct.Gdagvector(g)./PriorStruct.KbT)).^(1./PriorStruct.Nuvector(n))))).*PriorStruct.KbT)./(PriorStruct.xdagvector(x).*PriorStruct.K.*PriorStruct.breakingspeeds)).*...
                        (1-(((PriorStruct.Nuvector(n).*PriorStruct.K.*PriorStruct.breakingspeeds.*PriorStruct.lifetimes(i).*...
                        PriorStruct.xdagvector(x))./(PriorStruct.KbT))./(PriorStruct.Gdagvector(g)./PriorStruct.KbT))).^(1-(1/PriorStruct.Nuvector(n))))).*PriorStruct.breakingspeedprobability;

                    keep = imag(likelihood)==0;
                    likelihood = likelihood(keep);
                    keep = likelihood>0;
                    likelihood = likelihood(keep);
                    summedlikelihood = sum(likelihood,'omitnan');
                    if summedlikelihood ~= Inf
                        tempposteriorvector(total_loopcount) = summedlikelihood*PriorStruct.prior(total_loopcount);
                    else
                        tempposteriorvector(total_loopcount) = 0;
                    end
                end
            end
        end
        %update the prior
        normalizationterm = sum(tempposteriorvector,'omitnan');
        if normalizationterm ~= 0
            posteriorvector = tempposteriorvector;
            PriorStruct.prior = posteriorvector./normalizationterm;
        end 
    end
    
    %assign posterior vector
    posterior = PriorStruct.prior;

    %plot full posterior
    figure();
    scatter3(PriorStruct.xvector,PriorStruct.yvector,PriorStruct.zvector,80,posterior,'filled')
    xlabel('$x^\ddagger\:(nm)$','Interpreter','latex','FontSize',18)
    ylabel('$\Delta G^\ddagger\:(pN nm)$','Interpreter','latex','FontSize',18)
    zlabel('\nu','FontSize',18)
    title('Full Reference Posterior Distribution','FontSize',18)
    colormap(copper)
    colorbar
    a = colorbar;
    a.Label.Interpreter = 'latex';
    a.Label.String = '$p(\:x^\ddagger,\Delta G^\ddagger,\nu|\tau)$';
    a.FontSize = 18;
    
    [sorted_posterior,sorted_index] = sort(posterior,'descend');
    if ~isempty(sorted_posterior)
        for i = 1:length(sorted_posterior)
            if sum(sorted_posterior(1:i)) >= .95
                totalindex = i;
                break
            end
            totalindex = i;
        end
        
        cut_posterior = zeros(totalindex,1);
        cut_xvector = zeros(totalindex,1);
        cut_yvector = zeros(totalindex,1);
        cut_zvector = zeros(totalindex,1);

        for i = 1:totalindex
            cut_posterior(i) = posterior(sorted_index(i));
            cut_xvector(i) = PriorStruct.xvector(sorted_index(i));
            cut_yvector(i) = PriorStruct.yvector(sorted_index(i));
            cut_zvector(i) = PriorStruct.zvector(sorted_index(i));
        end
   
        %plot 95% credible interval
        figure();
        scatter3(cut_xvector,cut_yvector,cut_zvector,80,cut_posterior,'filled')
        xlabel('$x^\ddagger\:(nm)$','Interpreter','latex','FontSize',18)
        ylabel('$\Delta G^\ddagger\:(pN nm)$','Interpreter','latex','FontSize',18)
        zlabel('\nu','FontSize',18)
        title('95% Credible Interval','FontSize',18)
        colormap(copper)
        colorbar
        a = colorbar;
        a.Label.Interpreter = 'latex';
        a.Label.String = '$p(\:x^\ddagger,\Delta G^\ddagger,\nu|\tau)$';
        a.FontSize = 18;
        
        PosteriorStruct.posterior = posterior;
        PosteriorStruct.Nuvector = PriorStruct.Nuvector;
        PosteriorStruct.Gdagvector = PriorStruct.Gdagvector;
        PosteriorStruct.xdagvector = PriorStruct.xdagvector;
        PosteriorStruct.k0vector = PriorStruct.k0vector;
        PosteriorStruct.xvector = PriorStruct.xvector;
        PosteriorStruct.yvector = PriorStruct.yvector;
        PosteriorStruct.zvector = PriorStruct.zvector;
        PosteriorStruct.CredibleRegion = cut_posterior;
        PosteriorStruct.CredibleRegionxvector = cut_xvector;
        PosteriorStruct.CredibleRegionyvector = cut_yvector;
        PosteriorStruct.CredibleRegionzvector = cut_zvector; 
        
    else
        cut_posterior = [];
        cut_xvector = [];
        cut_yvector = [];
        cut_zvector = [];

        PosteriorStruct.posterior = posterior;
        PosteriorStruct.Nuvector = PriorStruct.Nuvector;
        PosteriorStruct.Gdagvector = PriorStruct.Gdagvector;
        PosteriorStruct.xdagvector = PriorStruct.xdagvector;
        PosteriorStruct.k0vector = PriorStruct.k0vector;
        PosteriorStruct.xvector = PriorStruct.xvector;
        PosteriorStruct.yvector = PriorStruct.yvector;
        PosteriorStruct.zvector = PriorStruct.zvector;
        PosteriorStruct.CredibleRegion = cut_posterior;
        PosteriorStruct.CredibleRegionxvector = cut_xvector;
        PosteriorStruct.CredibleRegionyvector = cut_yvector;
        PosteriorStruct.CredibleRegionzvector = cut_zvector;    
    end
end