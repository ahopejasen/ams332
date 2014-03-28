%%plotting function
function plotSimulation(timeVector,X,t_ode_cells,y_ode_cells,theParms,doLog)
if (~ exist('doLog','var'))
	doLog=0; %do secondary y-axis plot by default instead of semilog
end
IC=theParms.IC;
vs=theParms.vs;

numVars=size(IC,1);
numRuns=size(IC,2);
for theRun=1:numRuns
	figure();
	plot(timeVector(:,theRun).',X(:,:,theRun).'); %plot() wants column vectors
	hold on;
	if (theParms.N.doOde45)
		%plot ode45 as crosses
		%plot(t_ode(:,theRun),y_ode(:,:,theRun).','+')
		plot(t_ode_cells{theRun},y_ode_cells{theRun}(:,:),'+')
		%	t_ode(:,theRun),y_ode(:,n_ip,theRun).','g+', ...
		%	t_ode(:,theRun),y_ode(:,n_or,theRun).','r+', ...
		%	t_ode(:,theRun),y_ode(:,n_op,theRun).','y+');
	end % if doOde
	legTxt={'cI mRNA','cI protein','cro mRNA','cro protein'};
	legend(legTxt,'Location','East');
	% trick to do multiline title *with* variable values:
	% <cite>http://mechatronics.me.wisc.edu/labresources/MatlabTipsNTricks.htm</cite>
	% <cite>http://www.mathworks.com/matlabcentral/answers/93295</cite>
	titleTxt1={'Numerical solution to lysis gene model'};
	titleTxt2={	['Solid lines show forward Euler with time-step: ', ...
			num2str(theParms.N.timeStep),'s']}; 
    if (theParms.N.doOde45)
        titleTxt3={'crosses show MatLab ode45() solver'};
    else
        titleTxt3={''};
    end
	titleTxt4={['ICs: cI_{r}=',num2str(IC(vs.ir,theRun)),'  cI_{p}=',num2str(IC(vs.ip,theRun)), ...
			'  cro_{r}=',num2str(IC(vs.o_r,theRun)),  '  cro_{p}=',num2str(IC(vs.op,theRun))]; ...
		['\mu=',num2str(theParms.P.muCI),'  \omega=', num2str(theParms.P.wCI),'  \chi_{cI}=',num2str(theParms.P.xCIRna), ...
			'  \chi_{cro}=',  num2str(theParms.P.xCroRna), '  k=',num2str(theParms.P.kCI)]};
	hT=title([titleTxt1;titleTxt2;titleTxt3;titleTxt4]);
	fixTitle(hT);
    xlabel('time (s)');
	ylabel('concentration (molecules/cell)')
    hold off;

	%do plot to see smaller vars (semilogy or secondary axis)
	figure();
	if doLog
		%deal with zeros
		X1=X+1;
		semilogy(timeVector(:,theRun).',X1(:,:,theRun).'); %plot() wants column vectors
        hold  on;
		if (theParms.N.doOde45)
			%plot ode45 as crosses
			%plot(t_ode(:,theRun),y_ode(:,:,theRun).','+')
			y1=cellfun(@(x) x+1,y_ode_cells,'UniformOutput',false);
			semilogy(t_ode_cells{theRun},y1{theRun}(:,:),'+')
			%	t_ode(:,theRun),y_ode(:,n_ip,theRun).','g+', ...
			%	t_ode(:,theRun),y_ode(:,n_or,theRun).','r+', ...
			%	t_ode(:,theRun),y_ode(:,n_op,theRun).','y+');
		end % if doOde
		hold off;
		legend(legTxt,'Location','East');
		hT=title([titleTxt1;titleTxt2;titleTxt3;titleTxt4]);
		fixTitle(hT);
		xlabel('time (s)');
		ylabel('concentration (molecules/cell)')

	else %do plotyy to plot variable with largest values on a seperate axis
		figure();
		thisX=X(:,:,theRun); %seperating this out allows the nice max(thisX(:)) trick below
		%see <cite>http://stackoverflow.com/questions/2635120/how-can-i-find-the-maximum-or-minimum-of-a-multi-dimensional-matrix-in-matlab</cite>
		[~,xPos]=max(thisX(:)); %xPos uses 1-D (collum-major) indexing
		[maxVarRow,~]=ind2sub(size(thisX),xPos); %convert to 2-d indexing
		%extract column to plot it on secondary axis
		maxVar=thisX(maxVarRow,:);
		otherRows=~ismember(1:numVars,maxVarRow);
		otherVars=thisX(otherRows,:);
		[hax,hlines,hmax]=plotyy(timeVector(:,theRun).',otherVars.', ...
			timeVector(:,theRun).',maxVar.','plot','plot');
		titleTxt2={['Note second Y-scaling for ',num2str(maxVarRow),'th var']};
		%give secondary axis var a distinctive color
		set(hmax,'color',[1,0,1]); %magenta
		set(hmax,'marker','none');
		set(hmax,'lineStyle',':');
		set(hax(2),'ycolor',[1,0,1]);%magenta
		ht=title([titleTxt1;titleTxt2;titleTxt3]);
		fixTitle(hT);
		legend([hlines;hmax],[legTxt(otherRows),legTxt(maxVarRow)],'Location','East');
		xlabel('time (s)');
		ylabel(hax(1),'concentration (molecules/cell)')
		ylabel(hax(2),'concentration (molecules/cell)')
	end %if doLog
end %looping over runs to make plots
end %function plot


