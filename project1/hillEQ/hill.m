
%****************************************************************************
% Copyright (c) 2014,  A. Hope Jasentuliyana.  All rights reserved.
% This file is part of homework for Stony Brook University (SBU)  course
% AMS332, Spring 2014.
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
% This file free software: you can redistribute it and/or modify
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

%plotting Hill eq

%constants
h=[1 2 10]; %exponent of Hill eq (# of subunits)
numCoef=size(h,2); %number of hill coefficients
vMax=5.0; % units mM/s
k=20.0; % units mM
substrateStep=0.1; %step size for subtrate range

%use column vecor for substrate concentrations for easy plotting
% [:].' gives a col vector
substrateRange=[0:substrateStep:100].'; % units mM 
numSteps=size(substrateRange,1); %height of array

%result array, one column per hill coefficient
result=zeros(numSteps,numCoef);


for i=[1:numCoef] %iterate over hill exponents
	for j = [1:numSteps]
		conc=substrateRange(j); %current concentration
		result(j,i)=vMax*conc^i/(k^i + conc^i); %hill eq. Units: mM/s
	end %concentrations
		
end %hill 

plot(substrateRange,result);
title(['Reaction rate vs concentration for hill equation';'using different hill coefficients.';'Substrate concentration step size: ',num2str(substrateStep),' mM']);

%make an array for the legend
legTxt={''}; %curly brackets for cell array.
%cell arrays: http://www.mathworks.com/help/matlab/matlab_prog/create-a-cell-array.html

for i=[1:numCoef] 
	legTxt(i)={strcat('h=',num2str(h(i)))}; %need curly brackets 
end
legend(legTxt,'location','east');
xlabel('Concentration (mM)');
ylabel('Reaction rate (mM/s)');
