%% stochastic simulation of lysis/lysogeny using gilespie algorithm



%%% Note on my choice of data structures:
% Even though the algorithm only requires the current and next states,
% I'm keeping all time steps in memory to allow doing
% stats (mean, variance) across replicates, (and across time for stable states)
%

SHAREDIR='../../shared/';
STARTPATH=addpath(SHAREDIR);
% TODO: exit code/cleanup should restore STARTPATH


%% switches to control script behavior
if ~exist('fp','var') %output file handle
	fp=0; %stdout
else
	if fp<0 %invalid file
		fp=0;
	end
end
if ~exist('runsPerIC','var') %allows me to overide this from the workspace
    runsPerIC=20; %replicates for each IC 
end

if ~exist('doPhasePlots','var') %allows me to overide this from the workspace
    doPhasePlots=true; % plot log(_cro_{pro}_) vs 
                        % log(_cI_{pro}_) phase plot
end

if ~exist('doXYPlots','var') %allows me to overide this from the workspace
    doXYPlots=true; % plot data shows all four vars on linear scale
                    % in practice only the equilibrium protein
                    % will be visible above the X axis
end
if ~exist('doLogPlots','var') %allows me to overide this from the workspace
    doLogPlots=false; % plot data on semilogy scale
                    % this is good for seeing low concentration data
                    % as well as showing rna and protein on same plot
end

if ~exist('doSimulation','var') %allows me to overide this from the workspace
    doSimulation=true; %otherwise just (re)-compute stats and/or plot
			%existing data from workspace
end

if ~exist('maxStep','var')
	maxStep=50000; % maximum steps
end



%% initializations
if doSimulation
	%%% model variables


	% 5 vars x r replicate array x n steps
	numVars=6; %four chemicals, *time*, and reaction_number, for each replicate
	numChems=4; %number data colums that are actual chemicals (ie: not time or rxnNum)

	%indices for arrays (because numbers are hard to debug, and 
	%structs/classes of arrays might be slow compared to strict arrays
	n.ir=1; %cI rna
	n.ip=2; % cI protein
	n.o_r=3; %cro RNA
	n.op=4; %cro protein
	n.tm=5; %time steps
	n.rn=6; %reaction that happened at each step




	%%% initial conditions. 
	ic_none=1;
	ic_cro=2;
	ic_ci=3;
	ic_both=4;
	theIC=ic_none;
	IC(theIC,:)=zeros(1,numVars); 
	titleTxtIC{1}='All initial concentrations 0';
	%20 molecules of cro Rna
	theIC=ic_cro;
	IC(theIC,:)=zeros(1,numVars); 
	IC(theIC,n.o_r)=20; 
	titleTxtIC{theIC}='All initial concentrations 0 except cro_{rna}(0)=20';
	%20 molecules of cI Rna
	theIC=ic_ci;
	IC(theIC,:)=zeros(1,numVars); 
	IC(theIC,n.ir)=20; 
	titleTxtIC{theIC}='All initial concentrations 0 except cI_{rna}(0)=20';
	%20 molecules of botn cI and cro Rna
	theIC=ic_both;
	IC(theIC,:)=zeros(1,numVars); 
	IC(theIC,n.ir)=20;
	IC(theIC,n.o_r)=20; 
	titleTxtIC{theIC}='Both rna initially=20';
	% end of ICs
	numICs=size(IC,1);
	numRuns=runsPerIC * numICs;


	A_array=-1*ones(maxStep+1,numVars,numRuns); %add one for t=0 initial condition
	%REASON FOR the -1 inital array:
	%preallocating Acells speeds things up by 50%(!), but it gives erroneous 0 values
	%for runs that terminate early
	%a workaround is to populate with negative vals, since this is
	%a non-negative simulation. 
	%   this way the negative data can be filtered out later




	%%% model parameters
	P=struct ( ...
		'mu_ci',50, ... % (molecules/(cell*sec)) (synthesis constant for cI RNA)
		'x_ci_r',1.2, ... % (1/s)  (degradation constant for cI RNA)
		'w_ci',50, ... % (1/s) (synthesis constant for cI proten)
		'x_ci_p',1.2, ... % (1/s) (degradation constant for  cI protein)
		'k_cro',10, ... % (molecules/cell) (Cro protein concentration that occupies 50% of binding sites on cI dna)
		...	% note: this is not an error: k_cro is used in the equation for cI
		...
		... % cro equation
		'mu_cro',50, ... % (molecules/(cell*sec)) (synthesis constant for cro RNA)
		'x_cro_r',0.8, ... % (1/s)  (degradation constant for cro RNA)
		'w_cro',50, ... % (1/s) (synthesis constant for cro protein)
		'x_cro_p',0.8, ... % (1/s) (degradation constant for cro protein)
		'k_ci',10 ... %  (molecules/cell) (cI protein concentration that occupies 50% of binding sites on cro dna)
				... % this is not an error: k_cro is used in the equation for cI
	);

	tic();
	T=struct ( 'numRuns',numRuns,'maxStep',maxStep, ...
				'numVars',numVars, 'numICs',numICs, 'numChems',numChems);

	theParms=struct ( 'T',T , 'P',P,'IC',IC, 'n',n );

	if doSimulation
		A_array=runSimulation(A_array,theParms);
	end
end

%separate runs with different ICs
IC_idx=reshape(1:numRuns,numRuns/numICs,numICs)';
%the IC_idx holds the indexes A_array(:,:,IC_idx(theIC)) 
%for each IC condition

%% do aggregate stats
disp('last-iteration stats for each IC');
%TODO: this asssumes that no run ended early due to a0==0
%if a run does end early, there will be a bunch of terminating
% -1s that will get averaged in....

% endValStats=struct( eMean,[:,:] , ...
% 					eMedian,[:,:], ...
% 					eStd,[:,:], ...
% 					eMin,[:,:], ...
% 					eMax,[:,:]);
transgressors=zeros(numICs,runsPerIC); %array to hold the indices for
			%replicates that reached the unexpected equilibrium
titleTrans=cell(numICs); %text used later for figure titles

endStats=struct ('mean',zeros(numICs,numChems), ...	
				'median',zeros(numICs,numChems), ...
				'std',zeros(numICs,numChems), ...
				'min',zeros(numICs,numChems), ...
				'max',zeros(numICs,numChems)) ;
			
for theIC=1:numICs
	%endPointData=squeeze((A_array(end,1:numChems,IC_idx(theIC,:))))';
	%squeeze removes all sigleton dimensions, so  when
	%there is only one IC, endPointData is one dimensional which breaks
	%by code. using reshape instead...
	curRuns=IC_idx(theIC,:); %run numbers for current IC
	endPointData=A_array(end,1:numChems,curRuns);
	endPointData=reshape(endPointData,size(endPointData,2),size(endPointData,3));
	endPointData=endPointData';

	endStats.mean(theIC,:)=mean(endPointData);
	endStats.median(theIC,:)=median(endPointData);
	endStats.std(theIC,:)=std(endPointData);
	endStats.min(theIC,:)=min(endPointData);
	endStats.max(theIC,:)=max(endPointData);
 	
 	fprintf(fp,'%s\n',titleTxtIC{theIC});
	fprintf(fp,'------------------------------------------------------------\n');
	%% check for transgressor runs where the endpoint is at
	% the non-expected equilibrium
	% FOr each set of replicates for a given IC
	% check the median vals, and assume that if
	% median(cro_protein) > median(ic_protein), the expected equilibrium
	% is cro and vice versa (this could fail on edge-case runs)
	% Then, given we know what the expected equililbium is, find any
	% replicates where endpoint_cI_pro > endpoint_cro_pro (or vice versa)
	% These runs are the ones with unexpected equilibria: count them and 
	% flag them for plottting.

	if endStats.median(theIC,n.op)>endStats.median(theIC,n.ip) %op (cro) equil
		transgressors(theIC,:)=(endPointData(:,n.ip)>endPointData(:,n.op))';
		fprintf(fp,'%d/%d runs reached the unexpected cI equilibrium\n', ...
			sum(transgressors(theIC,:)),runsPerIC);
	else %ip (cI) equilibrium
		transgressors(theIC,:)=(endPointData(:,n.ip)<endPointData(:,n.op))';
		fprintf(fp,'%d/%d runs reached the unexpected cro equilibrium\n', ...
			sum(transgressors(theIC,:)),runsPerIC);
	end
	
	%%recalculate descriptive stats with transgressors removed
	expectedData= endPointData(~(transgressors(theIC,:)),:);
	endStats.mean(theIC,:)=mean(expectedData);
	endStats.median(theIC,:)=median(expectedData);
	endStats.std(theIC,:)=std(expectedData);
	endStats.min(theIC,:)=min(expectedData);
	endStats.max(theIC,:)=max(expectedData);
 	fprintf(fp,'mean:\t%s \n',num2str(endStats.mean(theIC,:)))
	fprintf(fp,'median :\t%s \n',num2str(endStats.median(theIC,:)))
	fprintf(fp,'std deviation:\t%s \n',num2str(endStats.std(theIC,:)))
	fprintf(fp,'*********************************************************\n');
% 	disp('min');
%	disp(endStats.min(theIC,:));
% 	disp('max');
% 	disp(endStats.max(theIC,:));
end

%% plot results
for theIC=1:numICs
    %protein vs protein
    %TODO: deal with non-positive data in log-log plot

    
    % title stings used in all plots
    titleTxt2=titleTxtIC{theIC};
    titleTxt3={['\mu=',num2str(P.mu_ci),'  \omega=', num2str(P.w_ci), ...
                '  \chi_{cI}=',num2str(P.x_ci_r), ...
                '  \chi_{cro}=',  num2str(P.x_cro_r), '  k=',num2str(P.k_ci)]};
            
            
    if doPhasePlots
    figure()
        for theRun=IC_idx(theIC,:)
            l_data=A_array(:,:,theRun); 
            l_data(0==l_data)=0.1;
            loglog(l_data(:, n.ip), l_data(:, n.op));
            hold on;
            titleTxt1={['cro_{pro} vs cI_{pro} for ',num2str(size(IC_idx,2)), ...
				' stochastic trials']; ...
				'0.1 added to zero values to allow log-scale plot'}	;
            title([titleTxt1;titleTxt2;titleTxt3]);
            ylabel ('molecules of cro protein');
			xlabel ('molecules of cI protein');

        end
	end




	titleTxt1={['Gilespie simulation of lysis/lysogeny gene model showing ', ...
		num2str(size(IC_idx,2)),' stochastic trials']};
	%this ia outside if doXYPlots, becuase logy plotting using this 
	%title also
	
	if doXYPlots
		figure()
		for theRun=IC_idx(theIC,:)
			plot(cumsum(A_array(:,n.tm,theRun)),A_array(:,1:numChems,theRun));
			%the time vector A_array(:,n.tm,theRun) contains the time
			%between reactions. So to plot conectrations over time
			%we use the cumsum() of the reaction times
			hold on;
		end
		legTxt={'cI mRNA','cI protein','cro mRNA','cro protein'};
		title([titleTxt1;titleTxt2;titleTxt3]);
		xlabel('time (s)');
		ylabel('number of molecules');
		legend(legTxt,'Location','East');
		hold off;
	end

	
    %semilogy
    if doLogPlots
		figure()
		for theRun=IC_idx(theIC,:)
			l_data=A_array(:,:,theRun); 
			l_data(0==l_data(:,1:numChems))=0.1;
			semilogy(cumsum(l_data(:,n.tm)),l_data(:,1:numChems));
				%the time vector A_array(:,n.tm,theRun) contains the time
				%between reactions. So to plot conectrations over time
				%we use the cumsum() of the reaction times
			xlabel('time (s)');
			ylabel('number of molecules');
			legend(legTxt,'Location','East');
			titleTxt1_5={	'0.1 added to zero values to allow log-scale plot'}	;
			title([titleTxt1;titleTxt1_5;titleTxt2;titleTxt3]);

			hold on;
		end
	end
    
end %loop over ICc

%do combined phase plot
if doPhasePlots
	figure()
	for theRun=1:numRuns
		l_data=A_array(:,:,theRun); 
		l_data(0==l_data)=0.1;
		loglog(l_data(:, n.ip), l_data(:, n.op));
		hold on;
	end
	titleTxt1={['cro_{pro} vs cI_{pro} for ',num2str(size(IC_idx,2)),' stochastic trials']; ...
			'0.1 added to zero values to allow log-scale plot'}	;
	titleTxt2={'showing all ICs'};
	titleTxt3={['\mu=',num2str(P.mu_ci),'  \omega=', num2str(P.w_ci), ...
				'  \chi_{cI}=',num2str(P.x_ci_r), ...
				'  \chi_{cro}=',  num2str(P.x_cro_r), '  k=',num2str(P.k_ci)]};
	title([titleTxt1;titleTxt2;titleTxt3]);
	ylabel ('molecules of cro protein');
	xlabel ('molecules of cI protein');
end

%% Part 4 degradation and switch from lysogeny 
% deterministic sim showed time-dependent switch taking 50s
% for xCI_pro=6.2, try 6.5 here for a faster response


%set lyosogenic ICs
P.x_ci_p=6.5;
%clear theParms.IC titleTxtIC
%clear-ing a member of a struct does not work

clear IC titleTxtIC
theIC=1;
IC(theIC,:)=zeros(1,numVars); %initialize all ICs to 0
IC(theIC,n.ir)=20; %update cI_{rna}(0) 
theParms.IC=IC;
clear IC
theParms.T.numICs=size(theParms.IC,1);
theParms.T.numRuns=runsPerIC * theParms.T.numICs;
titleTxtIC{theIC}='All initial concentrations 0 except cI_{rna}(0)=20';
%theParms.T.maxStep=100000; %simulation length
theParms.T.maxStep=10000; %simulation length (TODO: testing)
L_array=-1*ones(maxStep+1,theParms.T.numVars,theParms.T.numRuns); %add one for t=0 initial condition
L_array=runSimulation(L_array,theParms);
%% plot lysogenic decat simulation
