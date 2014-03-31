%% Simulation of a non-linear system of ODEs
% modeling deterministic phage lysis-lysogeny
%
% The model uses two genes cI and cro, whose protein products
% inhibit each other's rna transcription:
%   croProtein -| cIRna -> cIProtein   
%   cIProtein -| croRna -> croProtein
%
% Inhibition of RNA transcription is modeled as % dimeric binding 
% of inhibitive protein transcription factors to DNA 
% using Hill-equation kinetics:
%   dCIRna/dt= (1 - croProt^2/(k1 + croProt^2) ) - degradation
%   dCroRna/dt= (1 - CIProt^2/(k2 + CIProt^2) ) - degradation
%        where k1,k2 are constants.
%
% Production of Protein from RNA uses linear kinetics: 
%   dCIProt/dt= w1 CIRna - degradation; w1 is a constant
%   dCroProt/dt= w2 CroRna - degradation ; w2 is a constant
%
% Degradation kinetics for RNA and Protein for both genes are linear: 
%   eg: dCIRna/dt = production - X * CIRna, where X is a degradation constant
%
% A forward Euler algorithm is used to numerically approximate the system.
%
%<acknowlegments>
% All code in this file, unless otherwise noted in comments, is 
% written by the following contrbutors from AMS332:
% * A. Hope Jasentuliyana, SBU ID 100043659
%</acknowlegments>
%
% git tag: TODO: tag hand-in version
%
%<copyright>
%****************************************************************************
% Copyright (c) 2014,  A. Hope Jasentuliyana.  All rights reserved.
% This file is part of homework for Stony Brook University (SBU)  course
% AMS332, Spring 2014.
%
% Latest version of source should be available at: 
% http://github.com/ahopejasen/ams332.git
% 
% Resources used include class notes, handouts, and lectures given in
% SBU course AMS 332, Spring 2014 by:
% * Prof David Green
% * Prof Giancarlo La Camera
%
% This file is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This file is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this file.  If not, see <http://www.gnu.org/licenses/>.
%***************************************************************************/
%</copyright>
%% initializaton

% add path to shared functions from project1 
% WARN: assumes without checking that these functions will be in SIMDIR
% <cite>http://www.mathworks.com/matlabcentral/newsreader/view_thread/132505</cite>
SHAREDIR='../../shared/';
STARTPATH=addpath(SHAREDIR);
% TODO: exit code/cleanup should restore STARTPATH

% by definition everything is executed. To turn off specicic parts
% use doMesh=0; doPart1=0; doPart2=0; prior to executing script
if ~exist('doMesh','var')
	doMesh=1;
	doMesh;
end

if ~exist('doPart1','var')
	doPart1=1;
end

if ~exist('doPart2','var')
	doPart2=1;
end

if ~exist('doPart3','var')
	doPart3=1;
end

%% Initialize variables
P=struct();
IC=struct();
vs=struct();
N=struct();
theParms=struct('P',P, 'IC',IC','vs',vs,'N',N);

%% for notational convenience I'm describing the system using vector notation
% dX_dt=A(X); X(0)=X0; where X is a 4-vector, and A(X) is a vector function.
% the following constants allow me to use meaningful indices for the vectors:
n_ir=1; %cI rna
n_ip=2; % cI protein
n_or=3; %cro RNA
n_op=4; %cro protein
numVars=4;

%same as above, as a struct, to pass to functions
theParms.vs.ir=1;
theParms.vs.ip=2; 
theParms.vs.o_r=3;%avoiding 'or' keyword
theParms.vs.op=4;

%% define constant parameters of the simulated system
theParms.P=struct ( ...
	'muCI',50, ... % (molecules/(cell*sec)) (synthesis constant for cI RNA)
	'xCIRna',1.2, ... % (1/s)  (degradation constant for cI RNA)
	'wCI',50, ... % (1/s) (synthesis constant for cI proten)
	'xCIPro',1.2, ... % (1/s) (degradation constant for  cI protein)
	'kCro',10, ... % (molecules/cell) (Cro protein concentration that occupies 50% of binding sites on cI dna)
	...	% note: this is not an error: kCro is used in the equation for cI
	...
	... % cro equation
	'muCro',50, ... % (molecules/(cell*sec)) (synthesis constant for cro RNA)
	'xCroRna',0.8, ... % (1/s)  (degradation constant for cro RNA)
	'wCro',50, ... % (1/s) (synthesis constant for cro protein)
	'xCroPro',0.8, ... % (1/s) (degradation constant for cro protein)
	'kCI',10 ... %  (molecules/cell) (cI protein concentration that occupies 50% of binding sites on cro dna)
			... % this is not an error: kCro is used in the equation for cI
);



%% matrix of initial conditions
%rows are molecule types, cols are runs with different 
%starting conditions.
% initial conditions
%put all initial condintions into a struct
%rows are molecule types, cols are runs with different 
%starting conditions.
clear IC
%conditions for assignment
IC(n_ir,:)=[0,0,0,0]; % (molecules/cell) starting cI mRNA concentration
IC(n_ip,:)=[0,50,0,0]; % (molecules/cell) starting cI protein concentration
IC(n_or,:)=[0,0,20,0]; % (molecules/cell) starting cro mRNA concentration
IC(n_op,:)=[0,0,0,20]; % (molecules/cell) starting cro protein concentration

%further testing of dynamics of high concentrations
%IC(n_ir,:)=[20000,0,20000,5000]; % (molecules/cell) starting cI mRNA concentration
%IC(n_ip,:)=[0,2000,0,0]; % (molecules/cell) starting cI protein concentration
%IC(n_or,:)=[0,0,300,300]; % (molecules/cell) starting cro mRNA concentration
%IC(n_op,:)=[0,0,0,0]; % (molecules/cell) starting cro protein concentration
theParms.IC=IC;
clear IC;

%% numerical approximation parameters
odeOptions=odeset('NonNegative',[1,1,1,1]); %tells ode45 to force non-neg for all four vars
odeOptions=odeset(odeOptions,'RelTol',0.001); %default
odeOptions=odeset(odeOptions,'AbsTol',0.000001); %default
%odeOptions=odeset(odeOptions,'RelTol',0.000001);
%odeOptions=odeset(odeOptions,'AbsTol',0.00000001);
theParms.N=struct ( ...
	'startTime',0, ...  % (s)
	'endTime',25, ...  % (s) when absolute change over 5s < 1, rel change < 0.1%
	'timeStep',0.01, ...  % (s) delta_t for approximation TODO: testin TODO: testing
	'noNegative',1, ...  %constant passed to doForwardEuler that forces values to be non-negative
	'doOde45',1, ... % [0 , 1 , 2] = [euler only, both , ode45 only]
	'odeOptions',odeOptions ... %
);
clear odeOptions;

%clear result arrays 
%(the solver functions they are passed to will re-dimension them to the proper size)
clear('timeVector','X','dX','t_ode','y_ode','t_ode_cells','y_ode_cells','odeSol');


%% Project Part 1
%run the approximation algorithms, once for each set of IC's

if doPart1
	[timeVector,X,~,t_ode_cells,y_ode_cells]=runSimulation(theParms);

	%% plot results
	doLog=1; %request semilogy plot instead of secondary axis plot
	plotSimulation(timeVector,X,t_ode_cells,y_ode_cells,theParms,doLog); 
end
clear doLog;
%% Project Part 2
% $$ Plot cI_{prot} vs cro_{prot} for ICs 0-20 $$


% $$ Plot cI_{prot} vs cro_{prot} for ICs 0-2000 $$
if doPart2
	icStart=0;
	icEnd=2000;
	icStep=500;
	numICs=floor((icEnd-icStart)/icStep+1); %floor b/c if we have a partial step, it doesn't count

	%% create sequences of ICs
		%see: <cite>http://www.mathworks.com/matlabcentral/answers/16270-create-vector-of-repeating-elements-sort-of</cite>
	clear theParms.IC IC
	IC(n_ir,:)=repmat([icStart:icStep:icEnd],1,numICs); % 1 2 3 4 ... n 1 2 3 4 ... n 1 2 3 4 ... n (n times)
	IC(n_ip,:)=zeros(1,numICs^2);
	tmpV=repmat([icStart:icStep:icEnd],numICs,1);
	IC(n_or,:)=tmpV(:)' ; % 1 1 1 1 2 2 2 2 .... n n n n 
	IC(n_op,:)=zeros(1,numICs^2);
	theParms.IC=IC;
	clear IC;

	%% run simulation
	theParms.N.doOde45=1; % don't double check results
	theParms.N.timeStep=0.01;
	theParms.N.endTime=25;
	[timeVector,X,dX,t_ode_cells,y_ode_cells]=runSimulation(theParms);
	%doOde45: [0,1,2]==[only forward euler,both,only ode45]
	%% $$cro_{p} vs cI_{p} plot$$
	numRuns=size(timeVector,2) ; % timeVector is timeSteps x runs
	X1=X+1; %for semi-log or log-log plots
	figure();
	for theRun=1:numRuns
		loglog(X1(n_ip,:,theRun).',X1(n_op,:,theRun).');
		hold on;
		%plot(X(n_ip,:,theRun).',X(n_op,:,theRun).');
	end %runs
	hold off;

	%% plot all runs of different ICs on same plot
	plotAll(timeVector,X,t_ode_cells,y_ode_cells, theParms)
end %if doPart2

%% mesh plot of equilibrium vals (doMeshPlot recalculates data)
%% Shows that the boundary between the stable equlibria lies approximately
%% on the straight line $$ cro_{rna}(t=0) = 0.015 cI_{rna}(t=0) $$
%% The boundary appears to decrease in slope as cI_{rna} increases
%% showing that at high concentrations, the ratio of initial cI:cro
%% leading to the cI equilibrium increases.
%% It may continue to decelarate, or else approach a straight line.
if doMesh
	cISteps=[0:2:100 110:10:200]; %initial cI_{RNA}
	croSteps=[0:0.2:5 6:10];%initial cro{RNA}
	%cISteps=[0:2000:20000]; % $$initial cI_{RNA}$$
	%croSteps=0.015*cISteps;% $$initial cro{RNA}$$
	theParms.N.doOde45=2;
%	if is_octave %use forward euler for now to avoid bug
%		theParms.N.doOde45=0; % only use euler
%	else
%		theParms.N.doOde45=2; % only use ode45
%	end
	theParms.N.endTime=20;
	theParms.N.timeStep=0.01;
	[ir0,or0,croAnsEul,cIAnsEul,croAnsOde,cIAnsOde ]= ...
		doMeshPlot(theParms,croSteps,cISteps);
end %if doMesh

%% Part 3

if doPart3
	%% set up stable cI equilibrium
	% this was a real state from a simulation run
	theParms.IC(n_ir,1)=4.1660e+01;
	theParms.IC(n_ip,1)=1.7358e+03;
	theParms.IC(n_or,1)=2.0742e-03;
	theParms.IC(n_op,1)=1.2964e-01;
	theParms.IC=theParms.IC(:,1); %delete any other colums
	%run simulation for 25 more seconds and ensure stability
%	clear timeVector X dX t_ode_cells y_ode_cells
%	[timeVector , X, dX,t_ode_cells, y_ode_cells]=runSimulation(theParms);
%	%% verify straight line plot
%	plotSimulation(timeVector, X, t_ode_cells , y_ode_cells , theParms , 1);
%
%	%% verify relative and absolute change is negligible
%	relDiff(X,1,theParms.N.timeStep); % difference over 1s, simulation timestep 0.01
%	% results: abs change < 0.01/s rel change < 0.001% /s
%
	% now change simulation parameters to cause cI decay
	% Results: switch occurs near $\chi_{cI,Pro}=6.01175$
	% $\chi_{cI,Pro}=6.2 switches to cro after 50s $
	% $\chi_{cI,Pro}=6.05 switches to cro after 110s $
	% $\chi_{cI,Pro}=6.02 switches to cro after 230s $
	% $\chi_{cI,Pro}=6.015 switches to cro after 370s $
	% $\chi_{cI,Pro}=6.013 switches to cro after 560s $
	% $\chi_{cI,Pro}=6.012 switches to cro after 970s $
	% $\chi_{cI,Pro}=6.01175 switches to cro after 1370s $
	% $\chi_{cI,Pro}=6.0115 gives stable mix 260:7 cI:cro to 2000 secs$
	% $\chi_{cI,Pro}=6.011 gives stable mix 262:6 cI:cro to 2000 secs$
	% $\chi_{cI,Pro}=6.01 gives stable mix 263:7 cI:cro to 1000 secs$
	% $\chi_{cI,Pro}=6.0 gives stable mix 300:7 cI:cro to 300 secs$
	% xCIProVect=[6.011:0.001:6.014];

	%skip ode45 for now, too many crosses on plot
	theParms.N.doOde45=0;
	xCIProVect=[6.01175];
	%theParms.N.endTime=25;
	theParms.N.endTime=2000; %try a long sim
	clear timeVector X dX t_ode_cells y_ode_cells
	for theRun=[1:size(xCIProVect,2)]
		theParms.P.xCIPro=xCIProVect(theRun);
		%note it's ok to overwite timeVector since all runs have the same
		%but t_ode_cells can change
		[timeVector , X(:,:,theRun), dX,t_ode_cells{theRun}, ...
				y_ode_cells{theRun}]=runSimulation(theParms);
		plotSimulation(timeVector, X(:,:,theRun).', t_ode_cells{theRun} , y_ode_cells{theRun} , theParms , 1);
	end
end %if part3
