%constants use by forwardEulerTest.m
%<acknowlegments>
% All code in this file, unless otherwise noted in comments, is 
% written by the following contrbutors from AMS332:
% * A. Hope Jasentuliyana, SBU ID 100043659
%</acknowlegments>
%
% git tag: handin1.0 was handed in for grade credit.
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
%There is no good way to use constants in MATLAB
%this discussion:
%<cite>http://stackoverflow.com/questions/1773850/constants-in-matlab</cite>
%advocates using a classdef, although it says performance issues
%may or may not make this unsuitable.
%MatLAB docs <cite>http://www.mathworks.com/help/matlab/matlab_oop/properties-with-constant-values.html</cite>
% say "MATLAB evaluates the expressions when loading the class (when you first reference a constant property from that class)." Hopefully that's
% good enough in a global scope, s.t. there is no overhead for using
%classdef constants in a loop.
%
%this discussion shows the old way (using function calls):
%<cite>http://blogs.mathworks.com/loren/2006/09/13/constants/</cite>

classdef myConstants
	properties (Constant)
		NUM_STEPS=100; %steps to iterate
		TIME_STEP=0.01; 
		NUM_COLS=2; %number of dependent vars
		T_ZERO=0; %start time
		X_ZERO=[1,1]; % 1xNUM_COLS initial condition vector
		R=3; % x(t)=exp(R*t);
	end
end %myConstants
