function [ir0,or0,croAnsEul,cIAnsEul,croAnsOde,cIAnsOde ]=doMeshPlot(theParms,croSteps,cISteps)

vs=theParms.vs;


%% create x-y grid of ICs to allow surface plotting of equilibrium values vs ICs

[ir0,or0]=meshgrid(cISteps,croSteps); %each is a vector of ICs 
croAnsOde=[];
cIAnsOde=[];
cIAnsEul=[];
croAnsEul=[];
for croIdx=[1:size(croSteps,2)]
	for cIIdx=[1:size(cISteps,2)]
		theParms.IC=[ir0(croIdx,cIIdx);0;or0(croIdx,cIIdx);0];


		%% run the simulation for all ICs in theParms
        [tvect,Z,dZ,t_ode_cells,y_ode_cells]=runSimulation(theParms);
        if theParms.N.doOde45>0 %ode45 
            %y_ode_cells{} is contains  time x vars array
            croAnsOde(croIdx,cIIdx)=y_ode_cells{1}(end,vs.op);
            cIAnsOde(croIdx,cIIdx)=y_ode_cells{1}(end,vs.ip);
		end
        
		
		if  theParms.N.doOde45<2 %forward Euler data 
			%Z is vars x time
            croAnsEul(croIdx,cIIdx)=Z(vs.op,end); %selects last val of $cro_{prot}$
            cIAnsEul(croIdx,cIIdx)=Z(vs.ip,end); %selects last val of $cI_{prot}$
        end
	end
end
%% meshplot the end (equilibrium) values of protein vs cro and cI initial conditions
titleTxt1={'Boundary between cro (red) and cI (teal) equilibria'};
titleTxtEul={['Forward Euler: timeStep=',num2str(theParms.N.timeStep),'s']};
titleTxtOde={['ode45() relTol=',num2str(theParms.N.odeOptions.RelTol), ...
				'absTol=', num2str(theParms.N.odeOptions.AbsTol)]};
xlabelTxt='initial molecules cI_{RNA}';
ylabelTxt='initial molecules cro_{RNA}';
zlabelTxt='equibibrium molecules of protein';
%% plot forward euler data if requested
if theParms.N.doOde45<2
	figure
	mesh(ir0,or0,croAnsEul)
	hold on
	mesh(ir0,or0,cIAnsEul)
	hold off
	title([titleTxt1;titleTxtEul]);
	xlabel(xlabelTxt);
	ylabel(ylabelTxt);
	zlabel(zlabelTxt);
end

%% plot od45 data if requested
if theParms.N.doOde45>0
	figure
	mesh(ir0,or0,croAnsOde)
	hold on
	mesh(ir0,or0,cIAnsOde)
	hold off
	title([titleTxt1;titleTxtOde]);
	xlabel(xlabelTxt);
	ylabel(ylabelTxt);
	zlabel(zlabelTxt);
end %if
end %doMeshPlot
