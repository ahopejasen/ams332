%% vectorized hill function
% function result = hill(h,vMax,k,substrateRange)
% inputs: 
%		h (1xm)
%		vMax (scalar)
%		k (scalar)
%		substrateRange (n,1) 
% outputs:
%		result(n,m)

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
%h=[1 2 10]; %exponent of Hill eq (# of subunits)
%vMax=5.0; % units mM/s
%k=20.0; % units mM

%result array, one column per hill coefficient
%result=zeros(numSteps,numCoef);
function result = hill(h,vMax,k,substrateRange)
	subExp=substrateRange.^h; 
	result=vMax.*subExp ./(k.^h + subExp);%hill eq. Units: mM/s
end %hill()
