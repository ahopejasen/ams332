% forward Euler algotihm
% function [timeVect,X] doForwardEuler(dX_dtFunct,X0,startTime,endTime,timeStep);
% Returns:
%	X:  mxn matrix, m where corresponds to the number
%		of dependent variables in the system
%		and n=(startTime-endTime)/timeStep
%	timeVect: 1xn [start:timeStep:end] 
%
% Parameters:
%	dX_dtFunct is a mxvector of function handles corresponding to 
%		rates of change for each dependent var.
%	X0 is mx1 initial condition vector
%
% inputs are not checked for sanity, use at your own risk.

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

function [timeVect,X] =doForwardEuler(dX_dtFunct,X0,startTime,endTime,timeStep)


	%sanity checking 
	if ( (size(dX_dtFunct,1)<1) || (size(dX_dtFunct,2)~=1) || (size(dX_dtFunct) ~= size(X0)) )
		error('input vectors to doForwardEuler have invalid size');
	end

	if ( (endTime< startTime) || (endTime < (startTime+timeStep)) )
		error('invalid times passed to doForwardEuler: start %d, end %d, step %d',startTime,endTime,timeStep);
	end 

	%intializations (ALL_CAPS vars should be consts, but MatLab doesn't do const w/o classes)
	NUM_STEPS= floor((endTime-startTime)/timeStep); % WARNING ROUNDING ERROR HERE?
	NUM_VARS= size(X0,1); %number of dependent vars
	timeVect=[startTime:timeStep:endTime] ; % dependent variable


	X=zeros(NUM_VARS,NUM_STEPS);
	X(:,1)=X0;
	curSlope=zeros(NUM_VARS,1);
	%iterate       
	for (i=1:NUM_STEPS-1) %loop to find i+1th value based on ith val
		for j=[1:NUM_VARS]
			curSlope(j)=dX_dtFunct{j}(X(1,i),X(2,i));
			%curSlope(j)=dX_dtFunct{j}(X(:,i).');  %cant get this to work
		end
		X(:,i+1)=X(:,i) + curSlope*timeStep;	 %forward Euler algorithm using vector addition
	end

%debuggin

end %doForwardEuler()