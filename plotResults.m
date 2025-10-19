function f = plotResults(tt, SB, t, SF, tt2, relDif)
    fig1 = figure();
    pF1  = newplot(fig1);
    yyaxis left
    set(gca,'ColorOrderIndex',1);
    linestyleorder('default', 'beforecolor');
    plot(tt, SB, 'Color', "b", 'linewidth', 2); 
    hold(pF1, 'on');
    plot(t, SF, 'Color', "r", 'linewidth', 2); 
    xH = xlabel(pF1,'t');
    yH = ylabel(pF1,'$S_B(t)$', Interpreter='latex');
    xT = title(pF1,'The exercise boundary computed by two methods', Interpreter='latex');
    set(xT,'FontSize',12);
    set(xH,'FontSize',12);
    set(yH,'FontSize',12);
    set(gca,'FontSize',12);
    grid(pF1,'on');
    axis(pF1, 'tight');

    hold(pF1, 'on');
    set(gca,'ColorOrderIndex',1)
    L.AutoUpdate = 'on'; 
    yyaxis right
    plot(tt2, relDif, '--', 'Color', 'black', 'linewidth', 2); 
    yH = ylabel(pF1,'rel dif', Interpreter='latex');
    L = legend(pF1, {'$S_B(t)$', '$S_F(t)$', 'Rel. dif'}, 'Location', 'best', Interpreter='latex');
end