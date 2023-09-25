function PlotPriorDistribution(PriorStruct)

        %plot the reference prior
        figure();
        scatter3(PriorStruct.xvector,PriorStruct.yvector,PriorStruct.zvector,150,PriorStruct.prior,'filled')
        xlabel('$x^\ddagger\:(nm)$','Interpreter','latex','FontSize',20)
        ylabel('$\Delta G^\ddagger\:(pN nm)$','Interpreter','latex','FontSize',20)
        zlabel('\nu','FontSize',20)
        title('Prior Distribution','FontSize',20)
        colormap(copper)
        a = colorbar('Location','northoutside');
        set(gca,'ColorScale','log')
        a.Label.Interpreter = 'latex';
        a.Label.String = '$p(\:x^\ddagger,\Delta G^\ddagger,\nu)$';
        a.FontSize = 20;
        ax = gca;
        ax.FontSize = 20;
        set(gcf,'color','white');
end