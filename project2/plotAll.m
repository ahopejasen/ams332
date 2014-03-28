%% plot all runs of different ICs on same plot
function plotAll(timeVector,X,t_ode_cells,y_ode_cells, theParms)
	numRuns=size(X,3) ; % timeVector is timeSteps x runs
	icStart=theParms.IC(theParms.vs.ir,1);
	icEnd=theParms.IC(theParms.vs.ir,end);
	figure();
	for theRun=1:numRuns
		plot(timeVector(:,theRun),X(:,:,theRun).');
		hold on;
		if theParms.N.doOde45>0 %ode45 requested
			plot(t_ode_cells{theRun},y_ode_cells{theRun}(:,:),'+')
		end
	end %runs
	legTxt={'cI mRNA','cI protein','cro mRNA','cro protein'};
	legend(legTxt,'Location','East');
	titleTxt1={'Numerical solution to lysis gene model'};
	titleTxt2={ ['Solid lines show forward Euler with time-step: ', ...
				num2str(theParms.N.timeStep),'s'];};
	titleTxt3={['Different runs of same color show different ICs']; ...
				['ICs for RNA range from ',num2str(icStart),'-',num2str(icEnd),' with no initial protein']; ...
				['\mu=',num2str(theParms.P.muCI),'  \omega=', num2str(theParms.P.wCI),'  \chi_{cI}=',num2str(theParms.P.xCIRna), ...
					'  \chi_{cro}=',  num2str(theParms.P.xCroRna), '  k=',num2str(theParms.P.kCI)]};
	hT=title([titleTxt1;titleTxt2;titleTxt3]);
	xlabel('time (s)');
	ylabel('concentration (molecules/cell)')
	hold off;

	%plot all vars log-log
	X1=X+1; %to get rid of zeros. data already guarateed non-negative
	figure();
	for theRun=1:numRuns
		semilogy(timeVector(:,theRun),X1(:,:,theRun).');
		hold on;
		if theParms.N.doOde45>0 %ode45 requested
			y1=cellfun(@(x) x+1,y_ode_cells,'UniformOutput',false);
			semilogy(t_ode_cells{theRun},y1{theRun}(:,:),'+')
		end	
	end %runs
	clear X1;
	legend(legTxt,'Location','East');
	hT=title([titleTxt1;titleTxt2;titleTxt3]);
	xlabel('time (s)');
	ylabel('concentration (molecules/cell)')
	hold off;


