%% computes and plots abs relative differences over time in X array
% function [XrelDif]=relDif(X,diffTime)
function [XrelDiff,XDiff]=relDiff(X,diffTime,timeStep,titleCells)
% X is vars x time x runs
% diffTime is the time (in seconds) that differences should be computed
% timeStep is the step size of the simulation

%% check bounds
[numVars,dataPoints,numRuns]=size(X);
boundsErr=0;
if  (numVars<1 || dataPoints<2 || numRuns<1)
	boundsErr=1;
%else if (diffTime > (dataPoints./timeStep))
%%% diffTime greater  than length of simulation
%	boundsErr=1;
%else if 	(diffTime < timeStep)
%%% diffTime greater than/equal to timeStep
%	boundsErr=1;
end
%
if (~ boundsErr)
	steps=diffTime/timeStep; %num of steps for time lag

	titleTxtR={['Relative (percent) change over ',num2str(diffTime),' s.']};
	titleTxtA={['Absolute change over ',num2str(diffTime),' s.']};
	if exist('titleCells','var');
		titleTxtR=[titleTxtR;titleCells];
		titleTxtA=[titleTxtA;titleCells];
	end
	legTxt={'cI mRNA','cI protein','cro mRNA','cro protein'};

	XDiff=abs((X(:,1:end-steps,:)-X(:,steps+1:end,:)));
	XDenom=X(:,steps+1:end,:);
	XDenom(0==XDenom)=eps; %avoid divide by zero eps is smallest floating point
	XrelDiff=100*abs(XDiff./XDenom); %percentages
	timeVect=[1:size(XrelDiff,2)].*timeStep;
	timeVect=timeVect.';



	%absolute diffs
	figure
	for theRun=1:numRuns
		semilogy(timeVect,XDiff(:,:,theRun).');
		hold on;
	end
	hold off;
	title(titleTxtA);
	legend(legTxt);
	xlabel('time(s)');
	ylabel('relative change'); 



	%relative diffs
	figure
	for theRun=1:numRuns
		semilogy(timeVect,XrelDiff(:,:,theRun).');
		hold on;
	end
	hold off;
	title(titleTxtR);
	legend(legTxt);
	xlabel('time(s)');
	ylabel('relative change'); 

else  % bounds not OK
	disp('Bad bounds passed to XrelDiff()');
end
end %function
