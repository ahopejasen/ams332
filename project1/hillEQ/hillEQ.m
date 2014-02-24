%% Script to explore hill equation
%<acknowlegments>
% All code in this file, unless otherwise noted in comments, is 
% written by the following contrbutors from AMS332:
% * A. Hope Jasentuliyana, SBU ID 100043659
%</acknowlegments>
%
% git tag: handin1.0 was handed in for grade credit.
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

%OLD iterative code, now vectorized
	%lResult=zeros(numSteps,numCoef); %loop result
	%for i=[1:numCoef] %iterate over hill exponents
		%coef=h(i); %current hill coefficient
		%for j = [1:numSteps]
			%conc=substrateRange(j); %current concentration
			%lResult(j,i)=vMax*conc^coef/(k^coef + conc^coef); %hill eq. Units: mM/s
		%end %concentrations
	%end %coeficients 
%END OLD CODE

	
function hillEQ()
	%constants
	kDef=20.0; %(mM) default value
	vMaxDef=5.0; %(mM/s) default
	hDef=2; % Hill coefficient dfault	
	substrateStep=0.1; %step size for subtrate range
	theTitle={'Reaction rate vs concentration for hill equation'; ...
		'using different hill coefficients.'; ...
		['step size= ',num2str(substrateStep), ...
			'(mM); k=',num2str(kDef), ...
			'(mM); vMax=',num2str(vMaxDef),'(mM/s)'] };


	h=[1 2 10];%exponent of Hill eqn
	%we repeat solving eq and plotting it three times
	%so lets use a local function for all that...
	%Try passing it a struct containing all parameters needed
	stHill=struct( 'h',h, ... %exponent of Hill eq (# of subunits)
			'vMax',vMaxDef, ... % units mM/s
			'k',kDef, ... % units mM (concentration where rate=1/2 vMax)
			'sR',[0:substrateStep:100].', ... % units mM (col vector)
			'cTitle',theTitle, ... %plot title
			'legChar','h',... %character to use for legend
			'itVar', h ); %whichever parameter gets iterated over
					%goes here again so that hillPlot(stHill) will
					%know what to use for the legend


	%result array, one column per hill coefficient
	clear('vResult','hPlot');
	[hPlot,vResult]=hillPlot(stHill); %local function caclulates and  plots data
	title(theTitle); %todo: haven't figured out how to update 
			%a struct-array's cell array
			%so it's easier to set the title outside of hillPlot()
			%but I'm leaving stHill as a cell array to keep
			%compatibility



	%%part 2 iterate over k

	stHill(1).h=hDef; %will this clear the former array?
	stHill(1).k=[10.0 20.0 40.0]; % units mM (concentration where rate=1/2 vMax)
	stHill(1).legChar='k'; %character to use for legend
	stHill(1).itVar=stHill(1).k; %whichever parameter gets iterated over
	%stHill.cTitle=theTitle; ... %plot title
	clear('vResult','hPlot');
	[hPlot,vResult]=hillPlot(stHill); %local function caclulates and  plots data
	theTitle={'Reaction rate vs concentration for hill equation'; ...
		'using different k coefficients.'; ...
		['step size= ',num2str(substrateStep), ...
			'(mM); h=',num2str(hDef), ...
			'; vMax=',num2str(vMaxDef),'(mM/s)'] };
	title(theTitle);



	%%part 3 iterate over vMax

	stHill(1).k=20.0; % units mM (concentration where rate=1/2 vMax)
	stHill(1).vMax=[2.0 5.0 10.0]; %mM/s
	stHill(1).legChar='vMax'; %character to use for legend
	stHill(1).itVar=stHill(1).vMax; %whichever parameter gets iterated over
	%stHill.cTitle=theTitle; ... %plot title
	clear('vResult','hPlot');
	[hPlot,vResult]=hillPlot(stHill); %local function caclulates and  plots data
	theTitle={'Reaction rate vs concentration for hill equation'; ...
		'using different vMax coefficients.'; ...
		['Substrate concentration step size: ',num2str(substrateStep),' mM']; ...
		['h=',num2str(hDef),':k=',num2str(kDef),'(mM)'] };
	theTitle={'Reaction rate vs concentration for hill equation'; ...
		'using different vMax coefficients.'; ...
		['step size= ',num2str(substrateStep), ...
			'(mM); h=',num2str(hDef), ...
			'; k=',num2str(kDef),'(mM)'] };
	title(theTitle);
end %function hillEQ()

function [hHandle,vResults]=hillPlot(stP)
	myStp=stP(1); %because stP contains a cellarray and thus is a struct array
	vResults=hill(myStp.h,myStp.vMax,myStp.k,myStp.sR);
	hHandle=figure();
	plot(myStp.sR , vResults );

	%do titles outside of function for now
	%this cell array stuff is non-intuitive

	%todo: figure out how to redimension a structarray
	%and update the title cell array

	%fixme:also need to put in the multiline title fix
	%so the title doesn't disapear
	%can't figure a way to do this directly:
	%for c=1:size(stP(:))
		%theTitle{c}=stP(c).cTitle;
	%end;
	%title(theTitle); %adds title to current plot

	%make legend
	numCoef=size(myStp.itVar,2);
	for i=[1:numCoef]
		legTxt(i)={strcat(myStp.legChar,'=',num2str(myStp.itVar(i)))}; %need curly brackets
	end
	legend(legTxt,'location','eastOutside');
	xlabel('Concentration (mM)');
	ylabel('Reaction rate (mM/s)');
end %function hillPlot

	
