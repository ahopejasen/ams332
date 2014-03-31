function [ endStats, transgressors ] = plotSimulation( A_array, theParms, titleTxtIC, outFileHandle )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%unpack parameter using v2struct() from:
%<cite>http://blogs.mathworks.com/pick/2011/08/05/converting-variables-to-structures-and-vice-versa/</cite>
v2struct(theParms); %creates T,P,n,IC and any others
v2struct(T); %creates T.* as scalars
v2struct(n); %indexes into A_array(:,n_*,:);
v2struct(doFlags);

if (~exist('outFileHandle','var') || outFileHandle < 1) %invalid output file
   outFileHandle=1;%stdout 
end



%% do aggregate stats

%separate runs with different ICs
IC_idx=reshape(1:numRuns,numRuns/numICs,numICs)';
%the IC_idx holds the indexes A_array(:,:,IC_idx(theIC)) 
%for each IC condition


disp('last-iteration stats for each IC');
%TODO: this asssumes that no run ended early due to a0==0
%if a run does end early, there will be a bunch of terminating
% -1s that will get averaged in....

% endValStats=struct( eMean,[:,:] , ...
% 					eMedian,[:,:], ...
% 					eStd,[:,:], ...
% 					eMin,[:,:], ...
% 					eMax,[:,:]);
transgressors=zeros(numICs,runsPerIC); %array to hold the indices for
			%replicates that reached the unexpected equilibrium
titleTrans=cell(numICs); %text used later for figure titles

endStats=struct ('mean',zeros(numICs,numChems), ...	
				'median',zeros(numICs,numChems), ...
				'std',zeros(numICs,numChems), ...
				'min',zeros(numICs,numChems), ...
				'max',zeros(numICs,numChems)) ;
			
for theIC=1:numICs
	%endPointData=squeeze((A_array(end,1:numChems,IC_idx(theIC,:))))';
	%squeeze removes all sigleton dimensions, so  when
	%there is only one IC, endPointData is one dimensional which breaks
	%by code. using reshape instead...
	curRuns=IC_idx(theIC,:); %run numbers for current IC
	endPointData=A_array(end,1:numChems,curRuns);
	endPointData=reshape(endPointData,size(endPointData,2),size(endPointData,3));
	endPointData=endPointData';

	endStats.mean(theIC,:)=mean(endPointData);
	endStats.median(theIC,:)=median(endPointData);
	endStats.std(theIC,:)=std(endPointData);
	endStats.min(theIC,:)=min(endPointData);
	endStats.max(theIC,:)=max(endPointData);
 	
 	fprintf(outFileHandle,'%s\n',titleTxtIC{theIC});
	fprintf(outFileHandle,'------------------------------------------------------------\n');
	%% check for transgressor runs where the endpoint is at
	% the non-expected equilibrium
	% FOr each set of replicates for a given IC
	% check the median vals, and assume that if
	% median(cro_protein) > median(ic_protein), the expected equilibrium
	% is cro and vice versa (this could fail on edge-case runs)
	% Then, given we know what the expected equililbium is, find any
	% replicates where endpoint_cI_pro > endpoint_cro_pro (or vice versa)
	% These runs are the ones with unexpected equilibria: count them and 
	% flag them for plottting.

	if endStats.median(theIC,n.op)>endStats.median(theIC,n.ip) %op (cro) equil
		transgressors(theIC,:)=(endPointData(:,n.ip)>endPointData(:,n.op))';
		fprintf(outFileHandle,'%d/%d runs reached the unexpected cI equilibrium\n', ...
			sum(transgressors(theIC,:)),runsPerIC);
	else %ip (cI) equilibrium
		transgressors(theIC,:)=(endPointData(:,n.ip)<endPointData(:,n.op))';
		fprintf(outFileHandle,'%d/%d runs reached the unexpected cro equilibrium\n', ...
			sum(transgressors(theIC,:)),runsPerIC);
	end
	
	%%recalculate descriptive stats with transgressors removed
	expectedData= endPointData(~(transgressors(theIC,:)),:);
	endStats.mean(theIC,:)=mean(expectedData);
	endStats.median(theIC,:)=median(expectedData);
	endStats.std(theIC,:)=std(expectedData);
	endStats.min(theIC,:)=min(expectedData);
	endStats.max(theIC,:)=max(expectedData);
 	fprintf(outFileHandle,'mean:\t%s \n',num2str(endStats.mean(theIC,:)))
	fprintf(outFileHandle,'median :\t%s \n',num2str(endStats.median(theIC,:)))
	fprintf(outFileHandle,'std deviation:\t%s \n',num2str(endStats.std(theIC,:)))
	fprintf(outFileHandle,'*********************************************************\n');
% 	disp('min');
%	disp(endStats.min(theIC,:));
% 	disp('max');
% 	disp(endStats.max(theIC,:));
end

%% plot results
for theIC=1:numICs
    %protein vs protein
    %TODO: deal with non-positive data in log-log plot

    
    % title stings used in all plots
    titleTxt2=titleTxtIC{theIC};
    titleTxt3={['\mu=',num2str(P.mu_ci),'  \omega=', num2str(P.w_ci), ...
                '  \chi_{cI,}=',num2str(P.x_ci_p), ...
                '  \chi_{cro,p}=',  num2str(P.x_cro_p), '  k=',num2str(P.k_ci)]};
            
            
    if doPhasePlots
    figure()
        for theRun=IC_idx(theIC,:)
            l_data=A_array(:,:,theRun); 
            l_data(0==l_data)=0.1;
            loglog(l_data(:, n.ip), l_data(:, n.op));
            hold on;
            titleTxt1={['cro_{pro} vs cI_{pro} for ',num2str(size(IC_idx,2)), ...
				' stochastic trials']; ...
				'0.1 added to zero values to allow log-scale plot'}	;
            title([titleTxt1;titleTxt2;titleTxt3]);
            ylabel ('molecules of cro protein');
			xlabel ('molecules of cI protein');

        end
	end




	titleTxt1={['Gilespie simulation of lysis/lysogeny gene model showing ', ...
		num2str(size(IC_idx,2)),' stochastic trials']};
	%this ia outside if doXYPlots, becuase logy plotting using this 
	%title also
	
	if doXYPlots
		figure()
		for theRun=IC_idx(theIC,:)
			plot(cumsum(A_array(:,n.tm,theRun)),A_array(:,1:numChems,theRun));
			%the time vector A_array(:,n.tm,theRun) contains the time
			%between reactions. So to plot conectrations over time
			%we use the cumsum() of the reaction times
			hold on;
		end
		legTxt={'cI mRNA','cI protein','cro mRNA','cro protein'};
		title([titleTxt1;titleTxt2;titleTxt3]);
		xlabel('time (s)');
		ylabel('number of molecules');
		legend(legTxt,'Location','East');
		hold off;
	end

	
    %semilogy
    if doLogPlots
		figure()
		for theRun=IC_idx(theIC,:)
			l_data=A_array(:,:,theRun); 
			l_data(0==l_data(:,1:numChems))=0.1;
			semilogy(cumsum(l_data(:,n.tm)),l_data(:,1:numChems));
				%the time vector A_array(:,n.tm,theRun) contains the time
				%between reactions. So to plot conectrations over time
				%we use the cumsum() of the reaction times
			xlabel('time (s)');
			ylabel('number of molecules');
			legend(legTxt,'Location','East');
			titleTxt1_5={	'0.1 added to zero values to allow log-scale plot'}	;
			title([titleTxt1;titleTxt1_5;titleTxt2;titleTxt3]);

			hold on;
		end
	end
    
end %loop over ICc

%do combined phase plot
if doPhasePlots
	figure()
	for theRun=1:numRuns
		l_data=A_array(:,:,theRun); 
		l_data(0==l_data)=0.1;
		loglog(l_data(:, n.ip), l_data(:, n.op));
		hold on;
	end
	titleTxt1={['cro_{pro} vs cI_{pro} for ',num2str(size(IC_idx,2)),' stochastic trials']; ...
			'0.1 added to zero values to allow log-scale plot'}	;
	titleTxt2={'showing all ICs'};
	titleTxt3={['\mu=',num2str(P.mu_ci),'  \omega=', num2str(P.w_ci), ...
				'  \chi_{cI}=',num2str(P.x_ci_r), ...
				'  \chi_{cro}=',  num2str(P.x_cro_r), '  k=',num2str(P.k_ci)]};
	title([titleTxt1;titleTxt2;titleTxt3]);
	ylabel ('molecules of cro protein');
	xlabel ('molecules of cI protein');
end %if

end %function

