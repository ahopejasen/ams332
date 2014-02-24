%show divergence between forward Euler and ode45
%near critical point between zero and non-zero stability
%note, behavior of the two simulations coverge when
%doForwardEuler is given a timeStep<=0.005

% result array index constants 
i_rna=1;
i_pro=2;

% equation parameters
mu=1; % (mM/s) (synthesis constant for mRNA)
omega=1; % (1/s) (synthesis constant for protein: mM_protein/(s*mM_mRNA))
chi_r=1; % (1/s)  (degradation constant for RNA)
chi_p=1; % (1/s) (degradation constant for protein)
k=0.33; % (mM) (ligand concentration that occupies 50% of binding sites)

%set up equations
dr_dt=@(r,p) mu.*p.^2./(k.^2 + p.^2) - chi_r.*r; % mRNA rate of change
dp_dt=@(r,p) omega.*r - chi_p.*p; % protein rate of change
dX_dt={dr_dt; dp_dt}; % this is a cell-array of function handles
%setup for ode45
dy_dt=@(t,y) [mu *y(i_pro)^2/(k^2 + y(i_pro)^2) - chi_r * y(i_rna);omega * y(i_rna) - chi_p * y(i_pro)]; 

% numerical approximation parameters
startTime=0; % (s)
endTime=20; % (s)
timeStep=0.01; % (s) delta_t for approximation

%for p0=[0.21544:0.00000025:0.215441]
for p0=[0.5]
	X0=[0;p0];
	[timeVector,Z]=doForwardEuler(dX_dt,X0,0,endTime,timeStep);
	[t_ode,y_ode]=ode45(dy_dt,[0 endTime],X0);
	hold on;
	plot(timeVector ,Z);
	plot(t_ode,y_ode,'r-')
end
%do smaller timestep ode45 sim
clear t_ode y_ode;
[t_ode,y_ode]=ode45(dy_dt,[0 endTime],X0);
plot(t_ode,y_ode,'m-')
 ylabel ('concentration(mM)');
 xlabel('time (s)')
 title({'forward Euler (blue/green) and ode45 (red) simulations';'Protein(t=0)=[0.21544:0.00000025:0.215441]';'Euler timestep=0.01s, ode45 using defaults'});
 title({'forward Euler (blue/green) and ode45 simulations';'mRNA(t=0)=0';'Protein(t=0)=[0.21544:0.00000025:0.215441]';'Euler timestep=0.01s, ode45 using defaults'});
 hH=gcf
 %saveas(hH,'critical.png');

