%% runs main simulation loop
%A_array is (timesteps x variables x runs)
%it contains theParms.T.numVars variables (colums), of which
%theParms.T.numChems==(numVars-2) are actual chemical concentrations,
%and the last two are the time interval for a reaction to happen
%and the type of reaction, respectively.
%A_array's colums are indexed by by the structure theParms.n for
%easy mnemonic, (vs using numbers)
%	n.ir=1; %cI rna
%   n.ip=2; % cI protein
%   n.o_r=3; %cro RNA
%   n.op=4; %cro protein
%   n.tm=5; %time steps
%   n.rn=6; %reaction that happened at each step
function [A_array]=runSimulation(A_array,theParms)


IC=theParms.IC; %initial condition matirix: runs x vars
n=theParms.n; % indicies for variables: eg A_array(end,n.ir,1)== 
				% last timestep of ir (cI_{rna}) variable, first run
T=theParms.T;
P=theParms.P;


	%%% Eight reaction rates as a function of time and state
	% (two reactions (create and destroy) for each variable)
	numRxns=2*(T.numChems); %(the last variables, time rxnNum, dont count in reactions );
	V=zeros(1,numRxns); 

	%indexes for reactions rates
	v.ir_p=1; %cI rna positive (creation) rate
	v.ir_n=2; %cI rna negative (decay) rate
	v.ip_p=3; % cI protein positive (creation) rate
	v.ip_n=4; % cI protein negative (decay) rate
	v.or_p=5; %cro RNA positive (creation) rate
	v.or_n=6; %cro RNA negative (decay) rate
	v.op_p=7; %cro protein positive (creation) rate
	v.op_n=8; %cro protein negative (decay) rate


	%%% population update vectors for these eqs
	N=zeros(numRxns,T.numChems); % A(n+1,1:T.numChems)=A(n,1:T.numChems)+N(reaction_num_chosen,:);
	N(v.ir_p,n.ir)=1;
		%A=Acells{theRun}; %A( time, numVars ) %would be nice if this were a pointer/reference
		%A=Acells{theRun}; %A( time, numVars ) %would be nice if this were a pointer/reference
	N(v.ir_n,n.ir)=-1;
	N(v.ip_p,n.ip)=1;
	N(v.ip_n,n.ip)=-1;
	N(v.or_p,n.o_r)=1;
	N(v.or_n,n.o_r)=-1;
	N(v.op_p,n.op)=1;
	N(v.op_n,n.op)=-1;

	%% loop across replicates  

	fprintf('starting %d runs of %d iterations each\n',T.numRuns,T.maxStep);
	%flush the output paging buffer
	if is_octave
		fflush(stdout) ;
	end
	for theRun=1:T.numRuns


		%rotate through ICs as we go through runs
		theIC=ceil((theRun/T.numRuns)*T.numICs);
		A_array(1,:,theRun)=IC(theIC,:);



		A=A_array(:,:,theRun); %A( time, numVars ) %would be nice if this were a pointer/reference
						%but it's not: i'll need to update A_array(:,:,theRun)
						%when the simulation loop is done.

		a0=0; % *zero propensity*: probability that _some_ reaction happens

		V=zeros(1,numRxns); 

		curStep=1; %loop counter for simulation loop
		%% simulation loop
		notDone=true;
		while notDone

			%% current concentrations
			a_ir=A(curStep,n.ir); % current amount of $$cI_{rna}$$
			a_ip=A(curStep,n.ip); % current amount of $$cI_{pro}$$
			a_or=A(curStep,n.o_r); % current amount of $$cro_{rna}$$
			a_op=A(curStep,n.op); % current amount of $$cro_{pro}$$



			%% calculate V, the vector of reaction rates
			% V has dimensions (1,numRxns)
			% $cI_{rna}$
			% hill function $=\frac{a_op^2}{k_{cro}^2 + a_op^2}$
			V(v.ir_p)= P.mu_ci .* (1 - hill2(a_op,P.k_cro)); % cI_rna generator
			V(v.ir_n) = P.x_ci_r .* a_ir ; % cI_rna decay

			% $cI_{pro}$
			V(v.ip_p)= P.w_ci .* a_ir; % cI_pro generator
			V(v.ip_n)= P.x_ci_p .* a_ip; % cI_pro decay

			% $cro_{rna}$
			% hill function $=\frac{a_ip^2}{k_{cI}^2 + a_ip^2}$
			V(v.or_p)= P.mu_cro .* (1 - hill2(a_ip,P.k_ci)); % cro_rna generator
			V(v.or_n)= P.x_cro_r .* a_or ; % cro_rna decay

			% $cro_{rna}$
			V(v.op_p)= P.w_cro .* a_or; % cro_pro generator
			V(v.op_n)= P.x_cro_p .* a_op; % cro_pro decay

			% calculate a0= $\sum V_{j}
			a0=sum(V);

			if a0 %at least one reaction has non-zero probability...
				%%% calculate tNext
				u=rand();
				t_next=-log(u)./a0; 
				A(curStep+1,n.tm)=t_next; %keeps track of time steps in datavector A

				%%% calculate next reaction
				r=rand();
				V_sum=cumsum(V); 


				theRxnNum=find(V_sum > r.*a0,1,'first');
					% this _find()_ idiom is a hack from <cite>ams332 class text p. 169</cite>
					% TODO: this might fail if rounding error makes r.*a0 > V_sum
					% Safer to assign to a temp var, test if empty(), replace with
					% zero, and then assigning to theRxnNum.
				A(curStep+1,n.rn)=theRxnNum; % keeps track of which reactions happened


				%% update state of A
				A_chems_next=A(curStep,1:T.numChems)+N(theRxnNum,:);
				%set any negative concentrations equal to 0
				A_chems_next(A_chems_next<0)=0;
				A(curStep+1,1:T.numChems)=A_chems_next;

				curStep=curStep+1;
				%display something every thousand steps so we know it's not hung
				%if ~mod(curStep,5000)
				if ~mod(curStep,T.maxStep/5)
					fprintf('%d.',curStep);
					%flush the output paging buffer
					if is_octave
						fflush(stdout) ;
					end
				end



			else %  a0==0 no more possible transistions
				notDone=false; 
			end

			%%% check loop exit conditions
			notDone= (notDone && (curStep <= T.maxStep) && a0 > 0 ); % loop condition
		end % main simulation loop
	fprintf('\nfinished run %d/%d\n',theRun,T.numRuns);

	%update Acells
	%Acells{theRun}=A; 
	A_array(:,:,theRun)=A; 
	end % loop replicates loop


	toc()

end %run simulation
