%test of forward Euler algotihm

%****************************************************************************
% Copyright (c) 2014,  A. Hope Jasentuliyana.  All rights reserved.
% This file is part of homework for Stony Brook University (SBU)  course
% AMS332, Spring 2014.
% 
% All code in this file, unless otherwise noted in comments, is 
% written by A. Hope Jasentuliyana, SBU ID 100043659.
% Resources used include class notes, handouts, and lectures given in
% SBU course AMS 332, Spring 2014.
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




C=myConstants();
x=zeros(C.NUM_STEPS,C.NUM_COLS); % first col euler, second is analytic
t=zeros(C.NUM_STEPS,1);
t(1)=C.T_ZERO;
for (i=1:C.NUM_COLS)
       x(i)=C.X_ZERO(i);
end
       
       
for (i=1:C.NUM_STEPS-1) %loop to find i+1th value based on ith val
	t(i+1)=t(i)+C.TIME_STEP;
	curTime=t(i+1);
	prevVal=x(i);
	dx_dt=C.R*prevVal; % dx/dt=3x
	curVal=prevVal+dx_dt*C.TIME_STEP; %forward Euler algorithm
	x(i+1,1)=curVal;
	x(i+1,2)=exp(C.R * curTime); %analytic
end;

%compute sumsqaures error
ERR_COLS=C.NUM_COLS-1; 
errMtx=zeros(C.NUM_STEPS,ERR_COLS);
ssError=zeros(1,ERR_COLS); % (x(i,j+1) - x(i,j))^2
relError=zeros(C.NUM_STEPS,ERR_COLS); % (x(i,j+1)-x(i,j))/x(i,j)
for (i=1:(C.NUM_STEPS-1))
    for (j=1:ERR_COLS)
        curDiff=x(i,j+1)-x(i,j);
        relError(i,j)=curDiff/x(i,j);
        curError=curDiff^2;
        ssError(j)=ssError(j)+curError;
        errMtx(i,j)=curError;
    end % errcols
end

ssError

semilogy(t,x);
title('exponential function approximation')

plot(t,relError);
title('relative error as a function of time')
