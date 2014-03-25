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

% define constant parameters of the system
% cI equation
muCI=50; % (molecules/(cell*sec)) (synthesis constant for cI RNA)
xCIRna=1.2; % (1/s)  (degradation constant for cI RNA)
wCI=50; % (1/s) (synthesis constant for cI proten)
xCIPro=1.2; % (1/s) (degradation constant for  cI protein)
kCro=10; % (molecules/cell) (Cro protein concentration that occupies 50% of binding sites on cI dna)
	% note: this is not an error: kCro is used in the equation for cI

% cro equation
muCro=50; % (molecules/(cell*sec)) (synthesis constant for cro RNA)
xCroRna=0.8; % (1/s)  (degradation constant for cro RNA)
wCro=50; % (1/s) (synthesis constant for cro protein)
xCroPro=0.8; % (1/s) (degradation constant for cro protein)
kCI=10; %  (molecules/cell) (cI protein concentration that occupies 50% of binding sites on cro dna)
	% this is not an error: kCro is used in the equation for cI


%use anonymous function for second degree Hill equation (models dimers)
hill2=@(k,x) x.^2./(k.^2 + x.^2);


% initial conditions
cIRna=[0,0,0,0]; % (molecules/cell) starting cI mRNA concentration
cIPro=[0,50,0,0]; % (molecules/cell) starting cI protein concentration
croRna=[0,0,20,0]; % (molecules/cell) starting cro mRNA concentration
croPro=[0,0,0,20]; % (molecules/cell) starting cro protein concentration


numRuns=size(cIRna,2); 
numVars=4; %number of independent variables (cImRNA cIProtein, croRNA, croPRot)
%make sure initial condition vectors are same column size:
assert(size(cIRna,2)==numRuns && size(cIPro,2)== numRuns && ...
	 size(croRna,2)==numRuns && size(croPro,2)  == numRuns); 

% numerical approximation parameters
startTime=0; % (s)
endTime=20; % (s)
timeStep=0.001; % (s) delta_t for approximation
%timeStep=0.1; % FOR TESTING ONLY
%constant passed to doForwardEuler that forces values to be non-negative
noNegative=1;

%in case we want to change ode45() parameters:
odeOptions=odeset('NonNegative',[1,1,1,1]);
odeOptions.RelTol=0.001;
odeOptions.AbsTol=0.001;
%odeOptions.NonNegative=1;

% for notational convenience I describing the system using vector notation
% dX_dt=A(X); X(0)=X0; where X is a 4-vector, and A(X) is a vector function.
% the following constants allow me to use meaningful indices for the vectors:
n_ir=1; %cI rna
n_ip=2; % cI protein
n_or=3; %cro RNA
n_op=4; %cro protein

%put all initial condintions into a matrix
%rows are molecule types, cols are runs with different 
%starting conditions.
IC(n_ir,:)=cIRna;
IC(n_ip,:)=cIPro;
IC(n_or,:)=croRna;
IC(n_op,:)=croPro;

%clear result arrays 
%(the solver functions they are passed to will re-dimension them to the proper size)
clear('X','timeVector','t_ode','y_ode','t_ode_cells','y_ode_cells','odeSol');



%% Project Part 1
%run the approximation algorithms, once for each set of IC's
for theRun=1:numRuns

	% using anonymous function handles: 
	% <cite>http://www.mathworks.com/help/matlab/ref/function_handle.html</cite>
	% so we can pass the expression to different approximation algorithms if
	% needed....
		% NOTE: variable values (eg mu, k) are stored as constants
		% in the funtion handle when it is declared, and persist when
		% the handle is passed to another function
	dIr_dt=@(ciR,ciP,croR,croP) muCI.*(1- hill2(kCro,croP)) - xCIRna .*ciR;% cI mRNA rate of change 
			% *YES* it should be kCro and croP here, for the cI eq.
	dIp_dt=@(ciR,ciP,croR,croP) wCI.*ciR - xCIPro.*ciP; % cI protein rate of change
	dRr_dt=@(ciR,ciP,croR,croP) muCro.*(1- hill2(kCI,ciP)) - xCroRna .*croR; % cro mRNA rate of change
			% *YES* it should be kCI and ciP here, for the cro eq.
	dRp_dt=@(ciR,ciP,croR,croP) wCro.*croR - xCroPro.*croP; % cro protein rate of change

	%make an cell-array of the function handles, to pass to the solver
	dX_dt{n_ir,1}=dIr_dt;
	dX_dt{n_ip,1}=dIp_dt;
	dX_dt{n_or,1}=dRr_dt;
	dX_dt{n_op,1}=dRp_dt;

	%select the initial condition vector for this run from the IC matrix:
	X0=IC(:,theRun);
	
	%call doForwardEuler() to simulate a run
	[timeVector(:,theRun),X(:,:,theRun)]=doForwardEuler(dX_dt,X0,startTime,endTime,timeStep,noNegative);

	% same thing but compare with MatLab ODE45 (a runge-kutta method)
	dy_dt= @(t,y) [ muCI *(1 - hill2(kCro,y(n_op)) ) - xCIRna * y(n_ir); ...
			wCI * y(n_ir) - xCIPro * y(n_ip); ...
			muCro *(1 - hill2(kCI,y(n_ip)) ) - xCroRna * y(n_or); ...
			wCro * y(n_or) - xCroPro * y(n_op) ];
	if (is_octave) %running in octave
		%[t_ode(:,theRun),y_ode(:,:,theRun)]=ode45(dy_dt,[startTime,endTime],X0,odeOptions);
		%doesn't word because number of time steps subject to change across runs
		%need to use cell-arrays

		[t_ode,y_ode]=ode45(dy_dt,[startTime,endTime],X0,odeOptions);
		t_ode_cells{theRun}=t_ode(:);
		y_ode_cells{theRun}=y_ode(:,:);

	else %matlab specific
		odeSol=ode45(dy_dt,[startTime,endTime],X0,odeOptions);
		t_ode_cells{theRun}=startTime:0.5:endTime;
		y_ode_cells{theRun}=deval(odeSol,t_ode_cells{theRun});
	end %octave check
end %loop over initial conditions


%%plot results
for theRun=1:numRuns
	figure();
	hold on;
	plot(timeVector(:,theRun).',X(:,:,theRun).'); %plot() wants column vectors
	%plot ode45 as crosses
	%plot(t_ode(:,theRun),y_ode(:,:,theRun).','+')
	plot(t_ode_cells{theRun},y_ode_cells{theRun}(:,:),'+')
	%	t_ode(:,theRun),y_ode(:,n_ip,theRun).','g+', ...
	%	t_ode(:,theRun),y_ode(:,n_or,theRun).','r+', ...
	%	t_ode(:,theRun),y_ode(:,n_op,theRun).','y+');
	legTxt={'cI mRNA','cI protein','cro mRNA','cro protein'};
	legend(legTxt,'Location','East');
	% trick to do multiline title *with* variable values:
	% <cite>http://mechatronics.me.wisc.edu/labresources/MatlabTipsNTricks.htm</cite>
	% <cite>http://www.mathworks.com/matlabcentral/answers/93295</cite>
	titleTxt1={'Numerical solution to lysis gene model'};
	titleTxt2={	['Solid lines show forward Euler with time-step: ', ...
			num2str(timeStep),'s']; ...
		'crosses show MatLab ode45() solver';};
	titleTxt3={['ICs: cI_{r}=',num2str(cIRna(theRun)),'  cI_{p}=',num2str(cIPro(theRun)), ...
			'  cro_{r}=',num2str(croRna(theRun)),  '  cro_{p}=',num2str(croPro(theRun))]; ...
		['\mu=',num2str(muCI),'  \omega=', num2str(wCI),'  \chi_{cI}=',num2str(xCIRna), ...
			'  \chi_{cro}=',  num2str(xCroRna), '  k=',num2str(kCI)]};
	
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
end %main function



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
