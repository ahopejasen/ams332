%% Simulation of a non-linear system of ODEs
% modeling an autoregulatory gene whose dimeric protein product
% is a necesary transcription factor.
%
% Model uses the hill equation for co-operative binding of 
% mRNA transcription rates.
%
% A forward Euler algorithm is used to numerically approximate the system.
%
% todo: grouping plots for the legend is still not working
% in MatLab, although it works in Octave.
%<copyright>
%****************************************************************************
% Copyright (c) 2014,  A. Hope Jasentuliyana.  All rights reserved.
% This file is part of homework for Stony Brook University (SBU)  course
% AMS332, Spring 2014.
%
% Latest version of source should be available at: 
% https://github.com/ahopejasen/ams332.git
% 
% All code in this file, unless otherwise noted in comments, is 
% written by the following contrbutors from AMS332:
% * A. Hope Jasentuliyana, SBU ID 100043659
%
% Resources used include class notes, handouts, and lectures given in
% SBU course AMS 332, Spring 2014 by:
% * Prof David Green
% * Prof Giancarlo La Camera
% as well as MatLab documentation at: http://www.mathworks.com/help/matlab
%
% The code was developed in whole or in part using GNU Octave
% and testing for MatLab compatibility
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
%todo: some warnings in MatLab:
%Warning: Plot empty. 
%> In legend at 286
  %In autoReg at 213 
%Warning: Plot empty. 
%> In legend at 286
  %In autoReg at 210 

%%initializaton

% equation parameters
mu=1; % (mM/s) (synthesis constant for mRNA)
omega=1; % (1/s) (synthesis constant for protein: mM_protein/(s*mM_mRNA))
chi_r=1; % (1/s)  (degradation constant for RNA)
chi_p=1; % (1/s) (degradation constant for protein)
k=0.33; % (mM) (ligand concentration that occupies 50% of binding sites)

% initial conditions
r0=[0.5,0.0,0.0,0.2,0.5]; % (mM) starting mRNA concentration
p0=[0.5,0.2,0.5,0.0,0.0]; % (mM) starting protein concentration
numRuns=size(r0,2); 
numVars=2; %number of independent variables (mRNA & Protein)
assert(size(p0,2)==numRuns); %make sure initial condition vectors are same size

% numerical approximation parameters
startTime=0; % (s)
endTime=20; % (s)
timeStep=0.01; % (s) delta_t for approximation

% result array index constants 
i_rna=1;
i_pro=2;

%clear result arrays (the solver functions will dimension them)
clear('X','timeVector','t_ode','y_ode');

%in case we want to change ode45() parameters:
odeOptions=odeset();
%odeOptions.RelTol=0.001;
%odeOptions.AbsTol=0.001;
%odeOptions.NonNegative=1;

%% Project Part 1 & 2
%run the approximation algorithms, once for each set of IC's
for theRun=1:numRuns
	X0=[r0(theRun);p0(theRun)]; % initial conditions

	% try using anonymous function handles: 
	% <cite>http://www.mathworks.com/help/matlab/ref/function_handle.html</cite>
	% so we can pass the expression to different approximation algorithms if
	% needed....
		% NOTE: variable values (eg mu, k) are stored as constants
		% in the funtion handle when it is declared, and persist when
		% the handle is passed to another function
	dr_dt=@(r,p) mu.*p.^2./(k.^2 + p.^2) - chi_r.*r; % mRNA rate of change
	dp_dt=@(r,p) omega.*r - chi_p.*p; % protein rate of change

	dX_dt={dr_dt; dp_dt}; % this is a cell-array of function handles

	[timeVector(:,theRun),X(:,:,theRun)]=doForwardEuler(dX_dt,X0,startTime,endTime,timeStep);

	% compare with MatLab ODE45 (a runge-kutta method)
	dy_dt=@(t,y) [mu *y(i_pro)^2/(k^2 + y(i_pro)^2) - chi_r * y(i_rna);omega * y(i_rna) - chi_p * y(i_pro)]; 
	[t_ode(:,theRun),y_ode(:,:,theRun)]=ode45(dy_dt,[startTime:0.5:endTime],X0,odeOptions);
end %run


%%plot results
for theRun=1:numRuns
	figure();
	hold on;
	plot(timeVector(:,theRun).',X(:,:,theRun.')); %plot() wants column vectors
	%plot ode45 as crosses
	plot(t_ode(:,theRun),y_ode(:,i_rna,theRun),'b+',t_ode(:,theRun),y_ode(:,i_pro,theRun),'g+');
	legend({'mRNA','protein'},'Location','East');
	% trick to do multiline title *with* variable values:
	% <cite>http://mechatronics.me.wisc.edu/labresources/MatlabTipsNTricks.htm</cite>
	title({'Numerical solution to auto-regulatory gene model';['Solid lines show forward Euler with time-step: ',num2str(timeStep),' s'];'crosses show MatLab ode45() solver';['mu=',num2str(mu),'(mM/s) omega=',num2str(omega),'(1/s) chi\_r=',num2str(chi_r),'(1/s) chi\_p=',num2str(chi_p),'(1/s) k=',num2str(k),'(mM)']})
	xlabel('time (s)');
	ylabel('concentration (mM)')
end

%%project part 3: phase plane

clear('X','timeVector','icVector','numRuns');
%initial concentrations in mM
icVector=[0.0:0.1:0.2 0.2:0.2:1.4];%shows unstable equil better
%icVector=[0.0:0.5:1.4];%for quick testing
%icVector=[0.0:0.2:1.4]; %requested for assignment
numRuns=size(icVector,2)^2;
theRun=1; %index variable for X
figure();
hold on;
for r0=icVector %initial RNA concentrations in mM
	for p0=icVector %initial protein concentrations in mM
		X0=[r0;p0];
		[timeVector(:,theRun),X(:,:,theRun)]=doForwardEuler(dX_dt,[r0;p0],startTime,endTime,timeStep);
		%plot results
		plot(timeVector(:,theRun).',X(:,:,theRun).');		
		theRun=theRun+1;
	end
end
title({'Concentration vs time for autoregulatory model'; ...
	['for ',num2str(numRuns),' combinations of initial conditions']});
xlabel('time (s)');
ylabel('concentration (mM)');


%% phase plane plot

%phase portrait with color showing slope direction
%	red: protein & rna increasing (NorthEast)
%	yellow: protein decreasing, rna increasing (SE)
%	magenta: protein increasing, rna decreasing (NW)
%	yellow: protein decreasing, rna decreasing (SW)
%
%	Note: for legend tweaking, see:
%	<cite>http://www.mathworks.com/help/matlab/creating_plots/controlling-legends.html</cite>
figure();
hold on;

%graphics handles to group lines according to direction
%fixme: This is using ideas from the MatLab documentation for
%grouping legends, but although it works in Octave, it generates
%a warning and no legend in MatLab itself.
hNEGroup=hggroup;
hNWGroup=hggroup;
hSEGroup=hggroup;
hSWGroup=hggroup;
%clear vars used in loop
clear('XRNA','XProt','dRNA','dProtein','dRNAUp','dProteinUp');
clear('d_rup_pup','d_rup_pdn','d_rdn_pup','d_rdn_pdn');
clear('hHandle');
for theRun=1:numRuns
	XRNA=squeeze(X(i_rna,:,theRun)); %vector of RNA concentrations over time
	XProt=squeeze(X(i_pro,:,theRun));  %vector of protein concentration ove time
	dRNA=(dr_dt(X(i_rna,:,theRun),(X(i_pro,:,theRun)))); %vector of dRNA/dt
	dProtein=(dp_dt(X(i_rna,:,theRun),(X(i_pro,:,theRun)))); %vector of dProtein/dt
	dRNAUp=(dRNA>0); % logical array of positive RNA slopes
	dProteinUp=(dProtein>0); %logical array of positive Protein slopes
	d_rup_pup=dRNAUp & dProteinUp; %logical array of slopes to NE
	d_rup_pdn=dRNAUp & ~ dProteinUp; % SE slopes
	d_rdn_pup= ~ dRNAUp & dProteinUp; % NW slopes
	d_rdn_pdn= ~ dRNAUp & ~ dProteinUp; % SW slopes

	hHandle=plot(XProt(d_rup_pup).',XRNA(d_rup_pup).' ,'-r'); %phase plot of points with NE slope
	if (hHandle) %group it
		set(hHandle,'Parent',hNEGroup);
	end
	hHandle=plot( XProt(d_rup_pdn).',XRNA(d_rup_pdn).','-y'); %SE slopes
	if (hHandle) %group it
		set(hHandle,'Parent',hSEGroup);
	end
	hHandle=plot( XProt(d_rdn_pup).',XRNA(d_rdn_pup).','-m'); %NW
	if (hHandle) %group it
		set(hHandle,'Parent',hNWGroup);
	end
	hHandle=plot( XProt(d_rdn_pdn).',XRNA(d_rdn_pdn).','-b'); %SW
	if (hHandle) %group it
		set(hHandle,'Parent',hSWGroup);
	end
end


%set legends
legend('RNA up, Prot up (NE)', ...
	'RNA down, Prot up (SE)', ...
	'RNA up, Prot down (NW)', ...
	'RNA down, Prot down (SW)');
legend('Location','SouthOutside');

%plot nullclines
clear('nullcline','nullClineRange')

nullcline{i_rna}=@(p) mu .* p .^ 2 ./ (chi_r .* (k .^ 2 + p .^ 2)) ;
%corresponds to y=f(x)

nullcline{i_pro}=@(r) omega .* r ./  chi_p ;
%corresponds to x=f(y): need to plot flipped

nullClineRange=[0:0.001:1.4];
plot(nullcline{i_pro}(nullClineRange),nullClineRange,'-k'); 
plot(nullClineRange,nullcline{i_rna}(nullClineRange), '-k');
xlabel('protein (mM)');
ylabel('mRNA (mM)');
title({'Phase portrait of autoregulatory simulation'; ...
	'Trajectory color shows direction (see legend)'; ...
	'Nullclines in Black'});
%todo flesh-out chart titles
hold off;

