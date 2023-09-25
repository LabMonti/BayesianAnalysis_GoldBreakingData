function PlotPosteriorDistribution(PosteriorStruct)
    
        %plot full posterior
        figure();
        scatter3(PosteriorStruct.xvector,PosteriorStruct.yvector,PosteriorStruct.zvector,150,PosteriorStruct.posterior,'filled')
        ylabel('$\Delta G^\ddagger\:(pN nm)$','Interpreter','latex','FontSize',20)
        xlabel('$x^\ddagger\:(nm)$','Interpreter','latex','FontSize',20)
        zlabel('\nu','FontSize',20)
        title('Full Posterior Distribution','FontSize',20)
        colormap(copper)
        a = colorbar('Location','northoutside');
        set(gca,'ColorScale','log')
        a.Label.Interpreter = 'latex';
        a.Label.String = '$p(\:x^\ddagger,\Delta G^\ddagger,\nu|\tau)$';
        a.FontSize = 20;
        ax = gca;
        ax.FontSize = 20;
        set(gcf,'color','white');

        %plot 95% credible region
        figure();
        scatter3(PosteriorStruct.CredibleRegionyvector,PosteriorStruct.CredibleRegionxvector,PosteriorStruct.CredibleRegionzvector,150,PosteriorStruct.CredibleRegion,'filled')
        ylabel('$\Delta G^\ddagger\:(pN nm)$','Interpreter','latex','FontSize',20)
        xlabel('$x^\ddagger\:(nm)$','Interpreter','latex','FontSize',20)
        zlabel('\nu','FontSize',20)
        title('95% Credible Interval','FontSize',20)
        colormap(copper)
        a = colorbar('Location','northoutside');
        set(gca,'ColorScale','log')
        a.Label.Interpreter = 'latex';
        a.Label.String = '$p(\:x^\ddagger,\Delta G^\ddagger,\nu|\tau)$';
        a.FontSize = 20;
        ax = gca;
        ax.FontSize = 20;
        set(gcf,'color','white');
end