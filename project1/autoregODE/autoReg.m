% Simulation of a non-linear system of ODEs
% modeling an autoregulatory gene whose dimeric protein product
% is a necesary transcription factor.
%
% Model uses the hill equation for co-operative binding of 
% mRNA transcription rates.
%
% A forward Euler algorithm is used to numerically approximate the system.

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

%%global constants

% equation parameters
mu=1; % (mM/s) (synthesis constant for mRNA)
omega=1; % (1/s) (synthesis constant for protein: mM_protein/(s*mM_mRNA))
chi_r=1; % (1/s)  (degradation constant for RNA)
chi_p=1; % (1/s) (degradation constant for protein)
k=0.33; % (mM) (ligand concentration that occupies 50% of binding sites)

% initial conditions
r0=0.5; % (mM) starting mRNA concentration
p0=0.5; % (mM) starting protein concentration

% numerical approximation parameters
startTime=0; % (s)
endTime=20; % (s)
timeStep=0.01; % (s) delta_t for approximation

% simulation variables
X0=[r0;p0]; % initial conditions

% try using anonymous functions: 
% http://www.mathworks.com/help/matlab/ref/function_handle.html
% so we can pass them to different approximation algorithms if
% needed....
dr_dt=@(r,p) mu*p^2/(k^2 + p^2) - chi_r*r; % mRNA rate of change
dp_dt=@(r,p) omega*r - chi_p*p; % protein rate of change

dX_dt={dr_dt; dp_dt}; % this is a cell-array of function handles

[timeVector,X]=doForwardEuler(dX_dt,X0,startTime,endTime,timeStep);


%% off by one error from doForwardEuler
plot(timeVector(1,1:end-1).',X.'); %plot() wants column vectors
legend({'mRNA','protein'});



