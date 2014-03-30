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


	%% initializations

if doSimulation
	%%% model variables
	maxStep=50000; % maximum steps
	curStep=0;  % current step

	% 5 vars x r replicate array x n steps
	numVars=6; %four chemicals, *time*, and reaction_number, for each replicate
	numChems=4; %number data colums that are actual chemicals (ie: not time or rxnNum)

	%indices for arrays (because numbers are hard to debug, and 
	%structs/classes of arrays might be slow compared to strict arrays
	n_ir=1; %cI rna
	n_ip=2; % cI protein
	n_or=3; %cro RNA
	n_op=4; %cro protein
	n_tm=5; %time steps
	n_rn=6; %reaction that happened at each step




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
	IC(theIC,n_or)=20; 
	titleTxtIC{theIC}='All initial concentrations 0 except cro_{rna}(0)=20';
	%20 molecules of cI Rna
	theIC=ic_ci;
	IC(theIC,:)=zeros(1,numVars); 
	IC(theIC,n_ir)=20; 
	titleTxtIC{theIC}='All initial concentrations 0 except cI_{rna}(0)=20';
	%20 molecules of botn cI and cro Rna
	theIC=ic_both;
	IC(theIC,:)=zeros(1,numVars); 
	IC(theIC,n_ir)=20;
	IC(theIC,n_or)=20; 
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



	%%% Eight reaction rates as a function of time and state
	% (two reactions (create and destroy) for each variable)
	numRxns=2*(numChems); %(the last variables, time rxnNum, dont count in reactions );
	V=zeros(1,numRxns); 

	%indexes for reactions rates
	v_ir_p=1; %cI rna positive (creation) rate
	v_ir_n=2; %cI rna negative (decay) rate
	v_ip_p=3; % cI protein positive (creation) rate
	v_ip_n=4; % cI protein negative (decay) rate
	v_or_p=5; %cro RNA positive (creation) rate
	v_or_n=6; %cro RNA negative (decay) rate
	v_op_p=7; %cro protein positive (creation) rate
	v_op_n=8; %cro protein negative (decay) rate


	%%% population update vectors for these eqs
	N=zeros(numRxns,numChems); % A(n+1,1:numChems)=A(n,1:numChems)+N(reaction_num_chosen,:);
	N(v_ir_p,n_ir)=1;
	N(v_ir_n,n_ir)=-1;
	N(v_ip_p,n_ip)=1;
	N(v_ip_n,n_ip)=-1;
	N(v_or_p,n_or)=1;
	N(v_or_n,n_or)=-1;
	N(v_op_p,n_op)=1;
	N(v_op_n,n_op)=-1;

	tic();
	%% loop across replicates  

	fprintf('starting %d runs of %d iterations each\n',numRuns,maxStep);
	for theRun=1:numRuns


		%rotate through ICs as we go through runs
		theIC=ceil((theRun/numRuns)*numICs);
		A_array(1,:,theRun)=IC(theIC,:);



		%A=Acells{theRun}; %A( time, numVars ) %would be nice if this were a pointer/reference
		A=A_array(:,:,theRun); %A( time, numVars ) %would be nice if this were a pointer/reference
											%but it's not: i'll need to update A_array(:,:,theRun)
											%when the simulation loop is done.

		a0=0; % *zero propensity*: probability that _some_ reaction happens

		V=zeros(1,numRxns); 

		curStep=1; %loop counter for simulation loop
		%% simulation loop
		notDone=true;
		while notDone

			%% current concentrations
			a_ir=A(curStep,n_ir); % current amount of $$cI_{rna}$$
			a_ip=A(curStep,n_ip); % current amount of $$cI_{pro}$$
			a_or=A(curStep,n_or); % current amount of $$cro_{rna}$$
			a_op=A(curStep,n_op); % current amount of $$cro_{pro}$$



			%% calculate V, the vector of reaction rates
			% V has dimensions (1,numRxns)
			% $cI_{rna}$
			% hill function $=\frac{a_op^2}{k_{cro}^2 + a_op^2}$
			V(v_ir_p)= P.mu_ci .* (1 - hill2(a_op,P.k_cro)); % cI_rna generator
			V(v_ir_n) = P.x_ci_r .* a_ir ; % cI_rna decay

			% $cI_{pro}$
			V(v_ip_p)= P.w_ci .* a_ir; % cI_pro generator
			V(v_ip_n)= P.x_ci_p .* a_ip; % cI_pro decay

			% $cro_{rna}$
			% hill function $=\frac{a_ip^2}{k_{cI}^2 + a_ip^2}$
			V(v_or_p)= P.mu_cro .* (1 - hill2(a_ip,P.k_ci)); % cro_rna generator
			V(v_or_n)= P.x_cro_r .* a_or ; % cro_rna decay

			% $cro_{rna}$
			V(v_op_p)= P.w_cro .* a_or; % cro_pro generator
			V(v_op_n)= P.x_cro_p .* a_op; % cro_pro decay

			% calculate a0= $\sum V_{j}
			a0=sum(V);

			if a0 %at least one reaction has non-zero probability...
				%%% calculate tNext
				u=rand();
				t_next=-log(u)./a0; 
				A(curStep+1,n_tm)=t_next; %keeps track of time steps in datavector A

				%%% calculate next reaction
				r=rand();
				V_sum=cumsum(V); 


				theRxnNum=find(V_sum > r.*a0,1,'first');
					% this _find()_ idiom is a hack from <cite>ams332 class text p. 169</cite>
					% TODO: this might fail if rounding error makes r.*a0 > V_sum
					% Safer to assign to a temp var, test if empty(), replace with
					% zero, and then assigning to theRxnNum.
				A(curStep+1,n_rn)=theRxnNum; % keeps track of which reactions happened


				%% update state of A
				A_chems_next=A(curStep,1:numChems)+N(theRxnNum,:);
				%set any negative concentrations equal to 0
				A_chems_next(A_chems_next<0)=0;
				A(curStep+1,1:numChems)=A_chems_next;

				curStep=curStep+1;
				%display something every thousand steps so we know it's not hung
				if ~mod(curStep,5000)
					fprintf('%d.',curStep);
					%flush the output paging buffer
					if is_octave
						fflush(stdout) ;
					end
				end



			else %  a0==0 no more possible transistions
				notDone=false; 
			end

			%%% check loop exit conditions
			notDone= (notDone && (curStep <= maxStep) && a0 > 0 ); % loop condition
		end % main simulation loop
	fprintf('\nfinished run %d/%d\n',theRun,numRuns);

	%update Acells
	%Acells{theRun}=A; 
	A_array(:,:,theRun)=A; 
	end % loop replicates loop


	toc()
end %if do simulation

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
	endPointData=squeeze((A_array(end,1:numChems,IC_idx(theIC,:))))';
	endStats.mean(theIC,:)=mean(endPointData);
	endStats.median(theIC,:)=median(endPointData);
	endStats.std(theIC,:)=std(endPointData);
	endStats.min(theIC,:)=min(endPointData);
	endStats.max(theIC,:)=max(endPointData);
 	
 	fprintf('%s\n',titleTxtIC{theIC});
	fprintf('------------------------------------------------------------\n');
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

	if endStats.median(theIC,n_op)>endStats.median(theIC,n_ip) %op (cro) equil
		transgressors(theIC,:)=(endPointData(:,n_ip)>endPointData(:,n_op))';
		fprintf('%d/%d runs reached the unexpected cI equilibrium\n', ...
			sum(transgressors(theIC,:)),runsPerIC);
	else %ip (cI) equilibrium
		transgressors(theIC,:)=(endPointData(:,n_ip)<endPointData(:,n_op))';
		fprintf('%d/%d runs reached the unexpected cro equilibrium\n', ...
			sum(transgressors(theIC,:)),runsPerIC);
	end
	
	%%recalculate descriptive stats with transgressors removed
	expectedData= endPointData(~(transgressors(theIC,:)),:);
	endStats.mean(theIC,:)=mean(expectedData);
	endStats.median(theIC,:)=median(expectedData);
	endStats.std(theIC,:)=std(expectedData);
	endStats.min(theIC,:)=min(expectedData);
	endStats.max(theIC,:)=max(expectedData);
 	fprintf('mean:\t%s \n',num2str(endStats.mean(theIC,:)))
	fprintf('median :\t%s \n',num2str(endStats.median(theIC,:)))
	fprintf('std deviation:\t%s \n',num2str(endStats.std(theIC,:)))
	fprintf('*********************************************************\n');
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
            loglog(l_data(:, n_ip), l_data(:, n_op));
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
			plot(cumsum(A_array(:,n_tm,theRun)),A_array(:,1:numChems,theRun));
			%the time vector A_array(:,n_tm,theRun) contains the time
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
			semilogy(cumsum(l_data(:,n_tm)),l_data(:,1:numChems));
				%the time vector A_array(:,n_tm,theRun) contains the time
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
		loglog(l_data(:, n_ip), l_data(:, n_op));
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
