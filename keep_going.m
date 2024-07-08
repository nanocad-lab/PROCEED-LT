save(matname,'J','J_DB','DB','txtname','Xtxtname');
if(flag_plot)
    switch mode
        case 1 
            plot_x = log10(delay_weight)+J(:,1);
            plot_y = J(:,2);
        case 2
            plot_x = log10(delay_weight)+J(:,1);
            plot_y = J(:,5) ;
        case 3
            plot_x = J(:,5);
            plot_y = J(:,2);
    end
    h=plot(plot_x,plot_y,'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
    title(vtSelect + " " + vdSelect);
    drawnow;
end
sort_clear_J;
save(matname,'J','J_DB','DB','txtname','Xtxtname');
if(flag_plot)
    switch mode
        case 1 
            plot_x = log10(delay_weight)+J(:,1);
            plot_y = J(:,2);
        case 2
            plot_x = log10(delay_weight)+J(:,1);
            plot_y = J(:,5) ;
        case 3
            plot_x = J(:,5);
            plot_y = J(:,2);
    end
    h=plot(plot_x,plot_y,'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
    title(vtSelect + " " + vdSelect);
    drawnow;
end