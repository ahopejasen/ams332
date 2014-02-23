%% Script to explore hill equation
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
	substrateStep=0.1; %step size for subtrate range
	theTitle={'Reaction rate vs concentration for hill equation';'using different hill coefficients.';['Substrate concentration step size: ',num2str(substrateStep),' mM']};
	h=[1 2 10];%exponent of Hill eq (# of subunits) 

	%we repeat solving eq and plotting it three times
	%so lets use a local function for all that...
	%Try passing it a struct containing all parameters needed
	stHill=struct( 'h',h, ... %exponent of Hill eq (# of subunits)
			'vMax',5.0, ... % units mM/s
			'k',20.0, ... % units mM (concentration where rate=1/2 vMax)
			'sR',[0:substrateStep:100].', ... % units mM (col vector)
			'cTitle',theTitle, ... %plot title
			'legChar','h',... %character to use for legend
			'itVar', h ); %whichever parameter gets iterated over
					%goes here again so that hillPlot(stHill) will
					%know what to use for the legend

	%numCoef=size(stHill.h,2); %number of hill coefficients
		%^^DOesn't work because member 'cTitle' is a cell array,
		%which makes stHill a struct array
		%so we have to index into it. :-/
		%the non cell-array elements get broadcast, so it 
		%doesn't matter what index to use
	numCoef=size(stHill(1).h,2); %number of hill coefficients
	numSteps=size(stHill(1).sR,1); %height of array

	%result array, one column per hill coefficient
	clear('vResult','hPlot');
	[hPlot,vResult]=hillPlot(stHill); %local function caclulates and  plots data
	title(theTitle); %todo: haven't figured out how to update 
			%a struct-array's cell array
			%so it's easier to set the title outside of hillPlot()
			%but I'm leaving stHill as a cell array to keep
			%compatibility



	%%part 2 iterate over k

	stHill(1).h=2; %will this clear the former array?
	stHill(1).k=[10.0 20.0 40.0]; % units mM (concentration where rate=1/2 vMax)
	stHill(1).legChar='k'; %character to use for legend
	stHill(1).itVar=stHill(1).k; %whichever parameter gets iterated over
	%stHill.cTitle=theTitle; ... %plot title
	clear('vResult','hPlot');
	[hPlot,vResult]=hillPlot(stHill); %local function caclulates and  plots data
	theTitle={'Reaction rate vs concentration for hill equation';'using different k coefficients.';['Substrate concentration step size: ',num2str(substrateStep),' mM']};
	title(theTitle);



	%%part 3 iterate over vMax

	stHill(1).k=20.0; % units mM (concentration where rate=1/2 vMax)
	stHill(1).vMax=[2.0 5.0 10.0]; %mM/s
	stHill(1).legChar='vMax'; %character to use for legend
	stHill(1).itVar=stHill(1).vMax; %whichever parameter gets iterated over
	%stHill.cTitle=theTitle; ... %plot title
	clear('vResult','hPlot');
	[hPlot,vResult]=hillPlot(stHill); %local function caclulates and  plots data
	theTitle={'Reaction rate vs concentration for hill equation';'using different vMax coefficients.';['Substrate concentration step size: ',num2str(substrateStep),' mM']};
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

	
