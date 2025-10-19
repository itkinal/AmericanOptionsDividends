function f = plotData(tt, r, q, sig)
    fig1 = figure();
    pF1  = newplot(fig1);
    yyaxis left
    set(gca,'ColorOrderIndex',1);
    linestyleorder('default', 'beforecolor');
    plot(tt, r, 'Color', "b", 'linewidth', 2, 'DisplayName', '$r(t)$'); 
    hold(pF1, 'on');
    plot(tt, q, 'Color', "r", 'linewidth', 2, 'DisplayName', '$q(t)$'); 
    xH = xlabel(pF1,'t');
    yH = ylabel(pF1,'$r(t), q(t)$', Interpreter='latex');
    xT = title(pF1,'$r(t)$, $q(t)$ and $\sigma(t)$ as functions of the time $t$', Interpreter='latex');
    % xT = title(pF1,'');
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
    plot(tt, sig, 'Color', 'g', 'linewidth', 2, 'DisplayName','$\theta(t)$'); 
    yH = ylabel(pF1,'$\sigma(t)$', Interpreter='latex');
    L = legend(pF1, {'r(t)', 'q(t)', '$\sigma(t)$'}, 'Location', 'best', Interpreter='latex');
end
