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
function deterministic
%%initializaton

% add path to shared functions from project1 
% WARN: assumes without checking that these functions will be in SIMDIR
% <cite>http://www.mathworks.com/matlabcentral/newsreader/view_thread/132505</cite>
SHAREDIR='../shared/';
STARTPATH=addpath(SHAREDIR);
% TODO: exit code/cleanup should restore STARTPATH


% for notational convenience I describing the system using vector notation
% dX_dt=A(X); X(0)=X0; where X is a 4-vector, and A(X) is a vector function.
% the following constants allow me to use meaningful indices for the vectors:
n_ir=1; %cI rna
n_ip=2; % cI protein
n_or=3; %cro RNA
n_op=4; %cro protein

%same thing as a struct 
vs=struct (  ... %variable names
	'ir',1, ...
	'ip',2, ...
	'o_r',3, ... %avoid 'or' keyword
	'op',4 ...
);
% define constant parameters of the system
% cI equation
P=struct ( ...
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




% initial conditions
cIRna=[0,0,0,0]; % (molecules/cell) starting cI mRNA concentration
cIPro=[0,50,0,0]; % (molecules/cell) starting cI protein concentration
croRna=[0,0,20,0]; % (molecules/cell) starting cro mRNA concentration
croPro=[0,0,0,20]; % (molecules/cell) starting cro protein concentration
%put all initial condintions into a struct
%rows are molecule types, cols are runs with different 
%starting conditions.
%IC=struct ( ...
%	'ir',cIRna, ... %cI RNA
%	'ip',cIPro, ... %cI Protein
%	'or',croRna, ... %cro RNA
%	'op',croPro ... %cro Protein
%);
%%put all initial condintions into a matrix
%rows are molecule types, cols are runs with different 
%starting conditions.
IC(n_ir,:)=cIRna;
IC(n_ip,:)=cIPro;
IC(n_or,:)=croRna;
IC(n_op,:)=croPro;

% numerical approximation parameters
N=struct ( ...
	'startTime',0, ...  % (s)
	'endTime',20, ...  % (s)
	'timeStep',0.01, ...  % (s) delta_t for approximation
	'noNegative',1 ...  %constant passed to doForwardEuler that forces values to be non-negative
);

%clear result arrays 
%(the solver functions they are passed to will re-dimension them to the proper size)
clear('X','timeVector','t_ode','y_ode','t_ode_cells','y_ode_cells','odeSol');


doOde45=1; %use ode45() solver as a check
%% Project Part 1
%run the approximation algorithms, once for each set of IC's
[timeVector,X,t_ode_cells,y_ode_cells]=runSimulation(P,IC,vs,N,doOde45);

%%plot results
plotSimulation(timeVector,X,t_ode_cells,y_ode_cells,P,IC,vs,N,doOde45);
end %main function

function [timeVector,X,t_ode_cells,y_ode_cells]=runSimulation(P,IC,vs,N,doOde45);
	numRuns=size(IC,2); %IC: each row is a diff var, containing ICs for different runs
	%for theRun=1:numRuns

	%define simulation parameters
	numVars=size(IC,1); %number of independent variables (cImRNA cIProtein, croRNA, croPRot)
	%make sure initial condition vectors are same column size:
	%assert(size(cIRna,2)==numRuns && size(cIPro,2)== numRuns && ...
	%	 size(croRna,2)==numRuns && size(croPro,2)  == numRuns); 

	%in case we want to change ode45() parameters:
	odeOptions=odeset('NonNegative',[1,1,1,1]);
	odeOptions.RelTol=0.001;
	odeOptions.AbsTol=0.001;

	%use anonymous function for second degree Hill equation (models dimers)
	hill2=@(k,x) x.^2./(k.^2 + x.^2);

	% using anonymous function handles: 
	% <cite>http://www.mathworks.com/help/matlab/ref/function_handle.html</cite>
	% so we can pass the expression to different approximation algorithms if
	% needed....
		% NOTE: variable values (eg mu, k) are stored as constants
		% in the funtion handle when it is declared, and persist when
		% the handle is passed to another function
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

	for theRun=1:numRuns
		%select the initial condition vector for this run from the IC matrix:
		X0=IC(:,theRun);
		
		%call doForwardEuler() to simulate a run
		[timeVector(:,theRun),X(:,:,theRun)]=doForwardEuler(dX_dt,X0,N.startTime,N.endTime,N.timeStep,N.noNegative);

		% OPTIONAL.... for double checking results, use ode45():
		% same thing but compare with MatLab ODE45 (a runge-kutta method)
		% the idiom y(vs.ir) means the value of the previous step's cIRna.
		% y() is an array of the previous step's values.
		if (doOde45) 
			dy_dt= @(t,y) [ P.muCI *(1 - hill2(P.kCro,y(vs.op)) ) - P.xCIRna * y(vs.ir); ...
				P.wCI * y(vs.ir) - P.xCIPro * y(vs.ip); ...
				P.muCro *(1 - hill2(P.kCI,y(vs.ip)) ) - P.xCroRna * y(vs.o_r); ...
				P.wCro * y(vs.o_r) - P.xCroPro * y(vs.op) ];
			if (is_octave) %running in octave
				[t_ode,y_ode]=ode45(dy_dt,[N.startTime,N.endTime],X0,odeOptions);
				t_ode_cells{theRun}=t_ode(:);
				y_ode_cells{theRun}=y_ode(:,:);

			else %matlab specific
				odeSol=ode45(dy_dt,[N.startTime,N.endTime],X0,odeOptions);
				t_ode_cells{theRun}=N.startTime:0.5:N.endTime;
				y_ode_cells{theRun}=deval(odeSol,t_ode_cells{theRun});
			end %octave check
		else % no ode45 requested...  just create dummy variables to return
			t_ode_cells={};
			y_ode_cells={};
		end %if ode45
	end %loop over initial conditions
end %function runSimulation


%%plotting function
function plotSimulation(timeVector,X,t_ode_cells,y_ode_cells,P,IC,vs,N,doOde45);
numVars=size(IC,1);
numRuns=size(IC,2);
for theRun=1:numRuns
	figure();
	hold on;
	plot(timeVector(:,theRun).',X(:,:,theRun).'); %plot() wants column vectors
	if (doOde45)
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
			num2str(N.timeStep),'s']; ...
		'crosses show MatLab ode45() solver';};
	titleTxt3={['ICs: cI_{r}=',num2str(IC(vs.ir,theRun)),'  cI_{p}=',num2str(IC(vs.ip,theRun)), ...
			'  cro_{r}=',num2str(IC(vs.o_r,theRun)),  '  cro_{p}=',num2str(IC(vs.op,theRun))]; ...
		['\mu=',num2str(P.muCI),'  \omega=', num2str(P.wCI),'  \chi_{cI}=',num2str(P.xCIRna), ...
			'  \chi_{cro}=',  num2str(P.xCroRna), '  k=',num2str(P.kCI)]};
	
	hT=title([titleTxt1;titleTxt2;titleTxt3]);
	fixTitle(hT);
    xlabel('time (s)');
	ylabel('concentration (molecules/cell)')
    hold off;

	%do plotyy to plot variable with largest values on a seperate axis
	figure();
    thisX=X(:,:,theRun); %seperating this out allows the nice max(thisX(:)) trick below
	%see <cite>http://stackoverflow.com/questions/2635120/how-can-i-find-the-maximum-or-minimum-of-a-multi-dimensional-matrix-in-matlab</cite>
	[xLargest,xPos]=max(thisX(:)); %xPos uses 1-D (collum-major) indexing
	[maxVarRow,maxTimeCol]=ind2sub(size(thisX),xPos); %convert to 2-d indexing
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
end %looping over runs to make plots
end %function plot


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
