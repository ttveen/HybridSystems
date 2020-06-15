function  saveaspdf(fig,FILENAME)
% Save the figure as a pdf with correct dimensions

fig = gcf;

% make units equal
set(fig,'PaperUnits','centimeters');
set(fig,'Units','centimeters');

% get figure sizes 
pos = get(fig,'Position');
figWidth = pos(3);
figHeight = pos(4);
figSize = [figWidth, figHeight];

% set figure sizes
set(fig, 'PaperPosition', [0 0 figWidth figHeight]); %Position plot at left hand corner with width 5 and height 5.
set(fig, 'PaperSize', figSize); %Set the paper to have width 5 and height 5.

% save figure
print(fig,'-dpdf','-painters','-r600','-bestfit',FILENAME);
end



