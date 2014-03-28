%% Utilitiy functions

function fixTitle(hT)
	%Multiline title
	%trick:<cite>http://www.mathworks.com/matlabcentral/answers/93295</cite>
	%prevents title from being cut-off in matlab
	axpos = get(gca,'pos');
	set(hT,'units','normalized');
	extent = get(hT,'extent');
	set(gca,'pos',[axpos(1) axpos(2) axpos(3) axpos(4)-0.33*extent(4)]);
end

