function  plot_history(fig)
%
% Plot and save the history of the optimization objective and constraint.
% The user can costomize the figure sizes here.
global OPT


markertype = 'none';
markersize = 6;
linestyle = '-';

f = OPT.history.fval;

figure(fig); cla;
if OPT.stress_needed
    nsubplots = 3;
else
    nsubplots = 2;
end
subplot(nsubplots,1,1);
a = semilogy(f);
a.LineStyle = linestyle;
a.Marker = markertype;
a.MarkerSize = markersize;
title('objective history')
xlabel('iteration')
legend(OPT.functions.f{1}.name)

if isfield(OPT.history,'fconsval')
    g = reshape(OPT.history.fconsval,OPT.functions.n_func-1,[]).' ...
        + OPT.functions.constraint_limit;
    label = cell(1,OPT.functions.n_func-1);

    for i = 2:OPT.functions.n_func
       label{i-1} = OPT.functions.f{i}.name;
    end
    
    subplot(nsubplots,1,2); hold on
    plot(g,'Linestyle', linestyle, 'Marker', markertype,...
        'MarkerSize', markersize);
    title('constraint history');
    xlabel('iteration');

    hold on;
    cons_lim = 0*g + OPT.functions.constraint_limit;
    plot(cons_lim,'Color','r','LineStyle',':','LineWidth',1.5)
    hold off;
    legend(label);
    grange = abs(max(g)-min(g));
    ylim( [ min(min(g),OPT.functions.constraint_limit) - 0.2*grange, ...
        max(max(g),OPT.functions.constraint_limit) + 0.2*grange ]);
%     ylim([-.05,1.8*max(OPT.functions.constraint_limit)]);
end

% % Plot the true maximum stress history if there are stress constraints
% if OPT.stress_needed
%     subplot(nsubplots, 1, 3); hold on;
%     plot(OPT.history.true_stress_max,'Linestyle', linestyle, 'Marker', markertype,...
%         'MarkerSize', markersize);    
%     title('true maximum stress history');
%     xlabel('iteration');  
%     
%     hold on
%     cons_lim = OPT.parameters.slimit*ones(size(OPT.history.true_stress_max));
%     plot(cons_lim,'Color','r','LineStyle',':','LineWidth',1.5)
%     hold off
%     legend('true maximum stress'); 
%     srange = abs(max(OPT.history.true_stress_max)-min(OPT.history.true_stress_max));
%     ylim( [ min(min(OPT.history.true_stress_max),OPT.parameters.slimit)- 0.2*srange, ...
%         max(max(OPT.history.true_stress_max),OPT.parameters.slimit) + 0.2*srange ]);    
%     % ylim([-.05,1.2*max(OPT.history.true_stress_max)]);    
% end
