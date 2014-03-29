%% hill2 second degree hill eq
function [Y]=hill2(X,k)
	Xsq=X.^2;
	Y=Xsq./(k^2+Xsq);
end
