
f1 = gcf();

%set font for all subplots
for ax = get(f1,'Children');
    set(ax,'FontSize',14);
     set(ax,'LineWidth',2);
    lines = get(ax,'Children');
    for kk = 1:length(lines);
        try
        set(lines(kk),'LineWidth',2);
        catch
            fprintf('error\n')
        end
    end
end
