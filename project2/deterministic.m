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

%%initializaton

% add path to ode simulation functions from project1 
% WARN: assumes without checking that these functions will be in SIMDIR
% <cite>http://www.mathworks.com/matlabcentral/newsreader/view_thread/132505</cite>
SIMDIR='../project1/autoregODE/';
STARTPATH=addpath(SIMDIR);
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
xCroPro=0.8; % (1/s) (degradation constant for cro protein)
wCro=50; % (1/s) (synthesis constant for cro protein)
kCI=10; %  (molecules/cell) (cI protein concentration that occupies 50% of binding sites on cro dna)
	% this is not an error: kCro is used in the equation for cI


%use anonymous function for second degree Hill equation (models dimers)
hill2=@(k,x) x.^2./(k.^2 + x.^2);


% initial conditions
cIRna=[0,50,0]; % (molecules/cell) starting cI mRNA concentration
cIPro=[0,0,0]; % (molecules/cell) starting cI protein concentration
croRna=[0,0,20]; % (molecules/cell) starting cro mRNA concentration
croPro=[0,0,0]; % (molecules/cell) starting cro protein concentration


numRuns=size(cIRna,2); 
numVars=4; %number of independent variables (cImRNA cIProtein, croRNA, croPRot)
%make sure initial condition vectors are same column size:
assert(size(cIRna,2)==size(cIProt,2)==size(croRna,2)==size(croProt,2)==numRuns); 

% numerical approximation parameters
startTime=0; % (s)
endTime=20; % (s)
timeStep=0.01; % (s) delta_t for approximation

% for notational convenience I describing the system using vector notation
% dX_dt=A(X); X(0)=X0; where X is a 4-vector, and A(X) is a vector function.
% the following constants allow me to use meaningful indices for the vectors:
_ir=1; %cI rna
_ip=2; % cI protein
_or=3; %cro RNA
_op=4; %cro protein

%put all initial condintions into a matrix
%rows are molecule types, cols are runs with different 
%starting conditions.
IC(_ir)=cIRna;
IC(_ip)=cIPro
IC(_or)=croRna
IC(_op)=croPro;

%clear result arrays 
%(the solver functions they are passed to will re-dimension them to the proper size)
clear('X','timeVector','t_ode','y_ode');


%in case we want to change ode45() parameters:
odeOptions=odeset();
%odeOptions.RelTol=0.001;
%odeOptions.AbsTol=0.001;
%odeOptions.NonNegative=1;

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
	dIr_dt=@(ciR,croP) muCI.*(1- hill2(kCro,croP)) - xCIRna .*ciR;% cI mRNA rate of change 
			% *YES* it should be kCro and croP here, for the cI eq.
	dIp_dt=@(ciR,ciP) wCI.*ciR - xCIPro.*ciP; % cI protein rate of change
	dRr_dt=@(croR,croC) muCro.*(1- hill2(kCI,ciP)) - xCroRna .*croR; % cro mRNA rate of change
			% *YES* it should be kCI and ciP here, for the cro eq.
	dRp_dt=@(croR,croP) wCro.*croR - xCroPro.*croP; % cro protein rate of change

	%make an cell-array of the function handles, to pass to the solver
	dX_dt={_ir,1}=dIr_dt;
	dX_dt={_ip,1}=dIp_dt;
	dX_dt={_or,1}=dRr_dt;
	dX_dt={_op,1}=dRp_dt;

	%select the initial condition vector for this run from the IC matrix:
	X0=IC(theRun);
	
	%call doForwardEuler() to simulate a run
	[timeVector(:,theRun),X(:,:,theRun)]=doForwardEuler(dX_dt,X0,startTime,endTime,timeStep);

	% same thing but compare with MatLab ODE45 (a runge-kutta method)
	%dy_dt(_ir)=@(t,y) mu *(1 - hill2(kCro,y(_op)) ) - xCIRna * y(_ir);
	%dy_dt(_ip)=omega * y(i_rna) - chi_p * y(i_pro)]; 
	%[t_ode(:,theRun),y_ode(:,:,theRun)]=ode45(dy_dt,[startTime:0.5:endTime],X0,odeOptions);
end %run


%%plot results
for theRun=1:numRuns
	figure();
	hold on;
	plot(timeVector(:,theRun).',X(:,:,theRun.')); %plot() wants column vectors
	%plot ode45 as crosses
	%plot(t_ode(:,theRun),y_ode(:,i_rna,theRun),'b+',t_ode(:,theRun),y_ode(:,i_pro,theRun),'g+');
	legend({'mRNA','protein'},'Location','East');
	% trick to do multiline title *with* variable values:
	% <cite>http://mechatronics.me.wisc.edu/labresources/MatlabTipsNTricks.htm</cite>
	% <cite>http://www.mathworks.com/matlabcentral/answers/93295</cite>
	title({'Numerical solution to auto-regulatory gene model';['Solid lines show forward Euler with time-step: ',num2str(timeStep),' s'];'crosses show MatLab ode45() solver';['mu=',num2str(mu),'(mM/s) omega=',num2str(omega),'(1/s) chi\_r=',num2str(chi_r),'(1/s) chi\_p=',num2str(chi_p),'(1/s) k=',num2str(k),'(mM)']})
	xlabel('time (s)');
	ylabel('concentration (mM)')
end

%%project part 3: phase plane
%% phase plane plot
