%test of forward Euler algotihm

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

%globals (should be consts, but MatLab doesn't do const w/o classes)
NUM_STEPS= 100; %steps to iterate
TIME_STEP= 0.01; 
NUM_COLS= 2; %number of dependent vars
T_ZERO= 0; %start time
X_ZERO= [1,1]; % 1xNUM_COLS initial condition vector
R= 3; % x(t)=exp(R*t);

%initializations
x=zeros(NUM_STEPS,NUM_COLS); % first col euler, second is analytic
t=zeros(NUM_STEPS,1);
t(1)=T_ZERO;
for (j=1:NUM_COLS)
       x(1,j)=X_ZERO(j);
end
       
%iterate       
for (i=1:NUM_STEPS-1) %loop to find i+1th value based on ith val
	t(i+1)=t(i)+TIME_STEP;
	curTime=t(i+1);
	prevVal=x(i);
	dx_dt=R*prevVal; % dx/dt=3x
	curVal=prevVal+dx_dt*TIME_STEP; %forward Euler algorithm
	x(i+1,1)=curVal;
	x(i+1,2)=exp(R * curTime); %analytic
end;

%compute sumsqaures error
err_cols=NUM_COLS-1; 
errSqMtx=zeros(NUM_STEPS,err_cols); % squared error
ssError=zeros(1,err_cols); % (x(i,j+1) - x(i,j))^2
relError=zeros(NUM_STEPS,err_cols); % abs(x(i,j+1)-x(i,j))/x(i,j)
for (i=1:NUM_STEPS)
    for (j=1:err_cols)
        curDiff=abs(x(i,j+1)-x(i,j));
        relError(i,j)=curDiff/abs(x(i,j));
        curErrSquare=curDiff^2;
        ssError(j)=ssError(j)+curErrSquare;
        errSqMtx(i,j)=curErrSquare;
    end % errcols
end

ssError

%plot
semilogy(t,x,t,relError,'r');
legend('euler','analytic','relative error');
title({'e^{3t} vs t computed using forward euler and analytically.';'relative error is in red'});
