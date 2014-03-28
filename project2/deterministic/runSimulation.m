%% function [timeVector,X,dX,t_ode_cells,y_ode_cells]=runSimulation(theParms);
% the parameter theParms.N.doOde45 chooses algorithm: [0,1,2]=[forward euler only,both,odeOnly]
function [timeVector,X,dX,t_ode_cells,y_ode_cells]=runSimulation(theParms)
	P=theParms.P;
	IC=theParms.IC;
	vs=theParms.vs;
	N=theParms.N;

	numRuns=size(IC,2); %IC: each row is a different var, containing ICs for different runs

	%define simulation parameters
	numVars=size(IC,1); %number of independent variables (cImRNA cIProtein, croRNA, croPRot)

	%in case we want to change ode45() parameters:
	%odeOptions=odeset('NonNegative',[1,1,1,1]);
%	if N.noNegative
%		nNegArray=ones(1,numVars);
%	else
%		nNegArray=zeros(1,numVars);
%	
%	end
%	N.odeOptions.NonNegative=nNegArray;
%	N.odeOptions.RelTol=0.000001;
%	N.odeOptions.AbsTol=0.00000001;

	%% Define the ODE system
    %use anonymous function for second degree Hill equation (models dimers)
    %%
    % 
    % $$hill(x)=x^2/(k^2 + x^2)$$
    % 
	
    hill2=@(k,x) x.^2./(k.^2 + x.^2);

	% using anonymous function handles: 
	% <cite>http://www.mathworks.com/help/matlab/ref/function_handle.html</cite>
	% so we can pass the expression to different approximation algorithms if
	% needed....
		% NOTE: variable values (eg mu, k) are stored as constants
		% in the funtion handle when it is declared, and persist when
		% the handle is passed to another function
        %%
        % 
        % $$d[cI_{RNA}/dt=\mu_{cI} (1 - hill([cro_{prot}])) - \chi_{RNA} * [cI_{RNA}]$$
        % 
	
    dIr_dt=@(ciR,ciP,croR,croP) P.muCI.*(1- hill2(P.kCro,croP)) - P.xCIRna .*ciR;% cI mRNA rate of change 
			% *YES* it should be kCro and croP here, for the cI eq.
	dIp_dt=@(ciR,ciP,croR,croP) P.wCI.*ciR - P.xCIPro.*ciP; % cI protein rate of change
	dRr_dt=@(ciR,ciP,croR,croP) P.muCro.*(1- hill2(P.kCI,ciP)) - P.xCroRna .*croR; % cro mRNA rate of change
			% *YES* it should be kCI and ciP here, for the cro eq.
	dRp_dt=@(ciR,ciP,croR,croP) P.wCro.*croR - P.xCroPro.*croP; % cro protein rate of change

	%make an cell-array of the function handles, to pass to the solver
	dX_dt{vs.ir,1}=dIr_dt;
	dX_dt{vs.ip,1}=dIp_dt;
	dX_dt{vs.o_r,1}=dRr_dt;
	dX_dt{vs.op,1}=dRp_dt;

    %% loop through the matrix of intial conditions
	for theRun=1:numRuns
		%select the initial condition vector for this run from the IC matrix:
		X0=IC(:,theRun);
		
		if (2==N.doOde45) %skip forward euler
			timeVector=0;
			X=0;
			dX=0;
		else
		%% call doForwardEuler() to simulate a run
			[timeVector(:,theRun),X(:,:,theRun),dX(:,:,theRun)]=doForwardEuler(dX_dt,X0,N.startTime,N.endTime,N.timeStep,N.noNegative);
		end

		% OPTIONAL.... for double checking results, or testing use ode45():
		% MatLab ODE45 (a runge-kutta method)
		% the idiom y(vs.ir) means the value of the previous step's cIRna.
		% y() is an array of the previous step's values.
		
        %% Call ode45() as a secondary algorithm if requested
        if (N.doOde45>0) 
			dy_dt= @(t,y) [ P.muCI *(1 - hill2(P.kCro,y(vs.op)) ) - P.xCIRna * y(vs.ir); ...
				P.wCI * y(vs.ir) - P.xCIPro * y(vs.ip); ...
				P.muCro *(1 - hill2(P.kCI,y(vs.ip)) ) - P.xCroRna * y(vs.o_r); ...
				P.wCro * y(vs.o_r) - P.xCroPro * y(vs.op) ];
			if (is_octave) %running in octave
				[t_ode,y_ode]=ode45(dy_dt,[N.startTime,N.endTime],X0,N.odeOptions);
				t_ode_cells{theRun}=t_ode(:);
				y_ode_cells{theRun}=y_ode(:,:);

			else %matlab specific
				odeSol=ode45(dy_dt,[N.startTime,N.endTime],X0,N.odeOptions);
				%deval lets you choose selected timepoints 
				t_ode_cells{theRun}=N.startTime:0.5:N.endTime;
				y_ode_cells{theRun}=deval(odeSol,t_ode_cells{theRun}).';
			end %octave check
		else % no ode45 requested...  just create dummy variables to return
			t_ode_cells={};
			y_ode_cells={};
		end %if ode45
	end %loop over initial conditions
end %function runSimulation
