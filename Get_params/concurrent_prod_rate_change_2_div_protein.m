% Find parameters for the threshold level of the other division protein
% at constriction. The program returns the loss function value (see paper 
% for how to calculate the loss function) for input parameter

% clear

function tot_rd = concurrent_prod_rate_change_2_div_protein(params, mean_exp, i)

i

cells=[];
% Input parameters
v_inp=1.4;
tau= 72.5; % generation time in minutes
gr1= 0.0087; % growth rate 1 in min^-1
gr2 = 0.0126; % growth rate 2 in min^-1
ngen=24; % number of generation
mtcd = 18.3; % time from constriction to division
del_tcd = 1.5; % change in constriction timing before and after overproduction
sbc = 1.8+0.69; % size added from birth to constriction
sbc_1= 1.8+params(1); % size added from birth to constriction for the other division protein
s0 = 1; % production rate of FtsZ
s0_1=s0;% prod rate of other division prot
k = 0.029; % min^-1
k_1=k;
sorig = s0;
sorig_1=s0_1;
mthresh = sbc*sorig/gr1;
mthresh_1 = sbc_1*sorig_1/gr1;
kp_cells= 5*tau; % number of cells averaged over in the experiment
ch_cd= 12*tau; % the bbd changes to dbbd um at this time in the exp.
is0 = 0.1; % initial overprodution
ds0 = 0.85; % new production rate of FtsZ 
time_dly = 11; % time delay before FtsZ overproduction
s0 = s0+is0; sorig = s0;

%-----Noise parameters
cvl1=0.2*gr1; cvl2=0.2*gr2; % std dev in growth rate
sigtcd = 0.2*mtcd; 
cvt2= 0.2; cvt2_1=0.2; % CV of production rate of FtsZ
cvtthr = 0.2*mthresh; cvtthr_1=0.2*mthresh_1;

%-----Outputs
td_pop=NaN(1,6000); % Td timing of the population
lb_pop=NaN(1,6000); % Length at birth of the population
ld_pop=NaN(1,6000); % Length at division of the population
rate_pop=NaN(1,6000); % rates of the population
cm_pop=NaN(1,6000); % conc of the matur prot.
c_pop=NaN(1,6000); % conc of the tot prot
lc_pop=NaN(1,6000); % Length at constriction of the population
tc_pop=NaN(1,6000); % Time at constriction of the population
cm_tc_pop=NaN(1,6000); % Conc of mat prot. at constriction of the population
c_tc_pop=NaN(1,6000); % Conc at constriction of the population
rep_pop = NaN(1,6000); % other div process limit
tot_pop = NaN(1,6000); % tot cells
times=NaN(1,6000); % timing of the simulation

lmin=gr1/(5);

for no_in=1:75
x= Cell_const_div(v_inp);
x.rate= gr1;

x.nthresh = mthresh + randn()*cvtthr;
if x.nthresh < x.n
    x.nthresh = mthresh + randn()*cvtthr;
end

x.nthresh_1 = mthresh_1 + randn()*cvtthr_1;
if x.nthresh_1 < x.n_1
    x.nthresh_1 = mthresh_1 + randn()*cvtthr_1;
end

%-----For keeping track of events

x.vNextConst = 0;

s_in=s0+randn()*cvt2*s0;
while s_in<0
    s_in=s0+randn()*cvt2*s0;
end
x.s= s_in;

s_in_1=s0_1+randn()*cvt2_1*s0_1;
while s_in_1<0
    s_in_1=s0_1+randn()*cvt2_1*s0_1;
end
x.s_1= s_in_1;

cells=[cells x];
end

%-----Advance in time for %gens generations
gens = ngen*tau; % total time in mins for which the
tStep = 0.5; %in units of mins
simTime = 0;
post_ch= 0; % counter tells if ch has been changed
cnt=1;
tot = 0;
rep = 0;

while simTime < gens
    %-----Step
    simTime = simTime + tStep;

     if simTime>ch_cd+time_dly && post_ch==0
        s0=s0+ds0;
        for ct=1:length(cells)
            x=cells(ct);
            x.s = x.s/sorig*s0;
        end
        mtcd=mtcd+del_tcd;
        post_ch= 1;
    end
    
    for ct=1:length(cells)
    x=cells(ct);
    x.t = x.t + tStep;
    x.v = x.v*exp(tStep*x.rate); % grow exponentially
    x.n = x.n+x.s*x.v*tStep;
    x.nm= x.nm+k*(x.n-x.nm)*tStep;
    x.n_1 = x.n_1+x.s_1*x.v*tStep;
    x.nm_1= x.nm_1+k_1*(x.n_1-x.nm_1)*tStep;
    x.tNextDivs = x.tNextDivs-tStep;
    end
    %-----Perform events
    
    cellsNew=[];

    for ct=1:length(cells)
    x=cells(ct);
    
    %-----Constriction
    %Upon accumulating enough volume
    if x.vNextConst==0 && x.n-x.nthresh>=0 && x.n_1-x.nthresh_1>=0 

        if x.n_1-x.nthresh_1<x.s_1*x.v*tStep
            rep = rep+1;
            x.rep_cell = 1;
            x.ftsz_cell= 0;
        elseif x.n-x.nthresh<x.s*x.v*tStep
            x.rep_cell = 0;
            x.ftsz_cell= 1;
%         else
%             "Both are not limiting or both are equally limiting"
%             break
        end
        x.tot_cell = 1;
        tot = tot+1;

        x.vNextConst = 1;

        x.vOfConst = [x.vOfConst [x.v; x.t-x.tLastDiv; x.nm/x.v; x.n/x.v]];

        x.tNextDivs= mtcd+randn()*sigtcd;
        while x.tNextDivs<0
            x.tNextDivs= mtcd+randn()*sigtcd;
        end

        x.rate = max(gr2 + randn()*cvl2,lmin);

    end

    %-----Division
    %Upon accumulating enough volume
    if x.tNextDivs<=0

        x.td_cell = x.t-x.tLastDiv;
        x.lb_cell = x.vb;
        x.ld_cell = x.v;
        x.cm_cell = x.nm/x.v;
        x.c_cell = x.n/x.v;
        x.rate_cell= x.rate;
        x.lc_cell = x.vOfConst(1,1);
        x.tc_cell = x.vOfConst(2,1);
        x.cm_tc_cell = x.vOfConst(3,1);
        x.c_tc_cell = x.vOfConst(4,1);
        
        %-----Update cell
        
        x.tLastDiv = x.t;
        x.vd = x.v;
        rc = 0.5+randn()*0.03;
        x.v = x.vd*rc;
        x.vb = x.v;
        x.n = x.n*0.5;
        x.nm= x.nm*0.5;
        x.n_1 = x.n_1*0.5;
        x.nm_1= x.nm_1*0.5;
        x.rec_data=1;
        x.rate = max(gr1 + randn()*cvl1,lmin);
        x.tNextDivs = nan;
        x.vOfConst(:,1)=[];

        x.nthresh = mthresh + randn()*cvtthr;
        if x.nthresh < x.n
            x.nthresh = mthresh + randn()*cvtthr;
        end
        s_in=s0+randn()*cvt2*s0;
        while s_in<0
            s_in=s0+randn()*cvt2*s0;
        end
        x.s= s_in;

        x.nthresh_1 = mthresh_1 + randn()*cvtthr_1;
        if x.nthresh_1 < x.n_1
            x.nthresh_1 = mthresh_1 + randn()*cvtthr_1;
        end
        s_in_1=s0_1+randn()*cvt2_1*s0_1;
        while s_in_1<0
            s_in_1=s0_1+randn()*cvt2_1*s0_1;
        end
        x.s_1= s_in_1;
        x.vNextConst = 0;
        
        %----------other daughter
        if simTime<=kp_cells
        y=copy(x);
        y.td_cell= NaN;
        y.lb_cell = NaN;
        y.ld_cell = NaN;
        y.rate_cell = NaN;
        y.cm_cell = NaN;
        y.c_cell = NaN;
        y.lc_cell = NaN;
        y.tc_cell = NaN;
        y.cm_tc_cell = NaN;
        y.c_tc_cell = NaN;
        y.rep_cell = 0;
        y.ftsz_cell= 0;
        y.tot_cell = 0;

        y.rec_data=0;
        y.v = y.vd*(1-rc);
        y.vb = y.v;
        y.rate= max(gr1 + randn()*cvl1,lmin);

        y.nthresh = mthresh + randn()*cvtthr;
        if y.nthresh < y.n
            y.nthresh = mthresh + randn()*cvtthr;
        end
        s_in=s0+randn()*cvt2*s0;
        while s_in<0
            s_in=s0+randn()*cvt2*s0;
        end
        y.s= s_in;
        
        y.nthresh_1 = mthresh_1 + randn()*cvtthr_1;
        if y.nthresh_1 < y.n_1
            y.nthresh_1 = mthresh_1 + randn()*cvtthr_1;
        end
        s_in_1=s0_1+randn()*cvt2_1*s0_1;
        while s_in_1<0
            s_in_1=s0_1+randn()*cvt2_1*s0_1;
        end
        y.s_1= s_in_1;

        cellsNew =[cellsNew y];
        end
    end

    
    
    end

    cells = [cells cellsNew];
    
    
    ind_in = find([cells.rec_data]==1);
    td_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).td_cell];
    lb_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).lb_cell];
    ld_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).ld_cell];
    cm_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).cm_cell];
    c_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).c_cell];
    lc_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).lc_cell];
    tc_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).tc_cell];
    cm_tc_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).cm_tc_cell];
    c_tc_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).c_tc_cell];
    rep_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).rep_cell];
    tot_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).tot_cell];
    rate_pop(cnt:cnt+length(ind_in)-1)= [cells(ind_in).rate_cell];
    times(cnt:cnt+length(ind_in)-1)= (simTime-ch_cd)*ones(1,length(ind_in));
    cnt=cnt+length(ind_in);
    for cnt_ind=1:length(ind_in)
        cells(ind_in(cnt_ind)).rec_data=0;
    end
    
end

mean_lc_bf_sim = mean(lc_pop(times>-400 & times<0));
mean_lc_af_sim = mean(lc_pop(times> 400 & times<800));

tot_rd = sqrt((mean_lc_bf_sim/mean_exp(1)-1)^2 + (mean_lc_af_sim/mean_exp(2)-1)^2);

end
