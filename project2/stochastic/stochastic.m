%% stochastic simulation of lysis/lysogeny using gilespie algorithm


%% initializations

%%% model variables
% Four concentrations and time 
% I'm keeping all time steps in memory to allow doing
% stats (mean, variance) across replicates and across time for stable states

%% simulation parameters
numRuns=1; %replicates... running simultaneously 
maxStep=50000; % maximum steps
curStep=0;  % current step

% 5 vars x r replicate array x n steps
numVars=5; %four chemicals and time, for each replicate
A=zeros(numVars,numRuns,1); %initial conditions
						%note: A(i,:,j) is a *ROW vector*
a0=zeros(numRuns,1); % *zero propensity* probability that _some_ reaction happens

runsInPlay=ones(numRuns,1); %filter for the cols still meeting loop conditions

%indices for arrays 
n_ir=1; %cI rna
n_ip=2; % cI protein
n_or=3; %cro RNA
n_op=4; %cro protein
n_tm=5; %time steps


%%% model parameters
mu_ci=1.2;  % $\mu_{cI}$ cI RNA production
mu_cro=1.2; % $\mu_{cro}$ cro RNA production
w_ci=1.2; % $\omega_{cI}$ cI protein production
w_cro=1.2; % $\omega_{cro}$ cro protein production
x_ci_r=1.2; % $\chi_{cI,rna}$ ci rna degradation
x_ci_p=1.2; % $\chi_{cI,pro}$ ci pro degradation
x_cro_r=0.8;% $\chi_{cro,rna}$ cro rna degradation
x_cro_p=0.8; % $\chi_{cro,pro}$ cro pro degradation
k_ci=10;
k_cro=10;



%%% Eight reaction rates as a function of time and state
% variable  x replicate
% two reactions (create and destroy) for each variable
numRxns=2*(numVars-1); %(the last variable, time, doesn't count );
V=zeros(2*numRxns,numRuns); 

%indexes for reactions rates
v_ir_p=1; %cI rna positive (creation) rate
v_ir_n=2; %cI rna negative (decay) rate
v_ip_p=3; % cI protein positive (creation) rate
v_ip_n=4; % cI protein negative (decay) rate
v_or_p=5; %cro RNA positive (creation) rate
v_or_n=6; %cro RNA negative (decay) rate
v_op_p=7; %cro protein positive (creation) rate
v_op_n=8; %cro protein negative (decay) rate


%%% population update vectors for these eqs
N=zeros(numRxns,numVars);
N(v_ir_p,n_ir)=1;
N(v_ir_n,n_ir)=-1;
N(v_ip_p,n_ip)=1;
N(v_ip_n,n_ip)=-1;
N(v_or_p,n_or)=1;
N(v_or_n,n_or)=-1;
N(v_op_p,n_op)=1;
N(v_op_n,n_op)=-1;

%% main loop
notDone=true;
while notDone

	%% current concentrations across replicates (row vectors)
	a_ir=A(n_ir,:,curStep); % current amount of $$cI_{rna}$$
	a_ip=A(n_ip,:,curStep); % current amount of $$cI_{pro}$$
	a_or=A(n_or,:,curStep); % current amount of $$cro_{rna}$$
	a_op=A(n_op,:,curStep); % current amount of $$cro_{pro}$$



	%% calculate V, the vector of reaction rates
	% V has dimensions (numRxns,numRuns)
	% $cI_{rna}$
	% hill function $=\frac{a_op^2}{k_{cro}^2 + a_op^2}$
	V(v_ir_p)= mu_ci .* (1 - hill2(a_op,k_cro)); % cI_rna generator
	V(v_ir_n) = x_ci_r .* a_ir ; % cI_rna decay

	% $cI_{pro}$
	V(v_ip_p)= w_ci .* a_ir; % cI_pro generator
	V(v_ip_n)= x_ci_p .* a_ip; % cI_pro decay

	% $cro_{rna}$
	% hill function $=\frac{a_ip^2}{k_{cI}^2 + a_ip^2}$
	V(v_or_p)= mu_cro .* (1 - hill(a_ip,k_ci)); % cro_rna generator
	V(v_or_n)= x_cro_r .* a_or ; % cro_rna decay

	% $cro_{rna}$
	V(v_op_p)= w_cro .* a_or; % cro_pro generator
	V(v_op_n)= x_cro_p .* a_op; % cro_pro decay

	% calculate a0= $\sum V_{j}
	% V() is (numVars-1)xnumRuns
	% we want a0 to be numRuns
	a0=sum(V,1); % *row vector* of dim numRuns
	a0=a0.*runsInPlay ; % restores previous zero columns
	runsInPlay=a0>0;

	if any(a0(:)) %at least one reaction has non-zero probability...
		% calculate tNext
		u=rand(numRuns,1);
		t_next=-log(u)./rand(numRuns,1); 
		t_next(isequal(t_next,inf))=0; %handles division by zero!
			%<cite>www.mathworks.com/matlabcentral/newsreadr/view_thread/32253</cite>

		%%% calculate next reaction
		r=rand(1,numRuns);
		V_sum=cumsum(V,2*numVars,numRuns); 


		for idx=[1:numRuns]
			theRxnNum(idx)=find(V_sum(:,idx) > r.*a0(idx),1,'first');
				% a hack from ams332 class text p. 169
				% TODO: this might fail if rounding error makes r.*a0 > V_sum
				% safer to assign to a temp var, test if empty(), replace with
				% zero, and then assigning to theRxnNum.
				
		end
		
		%% update state of A
		curStep=curStep+1;



	else %  a0==0 no more possible transistions
		notDone=false; 
	end
	notDone= (notDone && (curStep <= maxStep) && any(a0(:)) > 0 ) % loop condition
end % main loop
