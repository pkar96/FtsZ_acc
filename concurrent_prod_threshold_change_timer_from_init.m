% Concurrent processes model is simulated over a population of cells. 
% The concurrent process is at the constriction. Constriction is controlled
% by a timer from initiation and accumulation of an division protein from 
% birth. Division happens after a constant time from constriction. Growth 
% rate is different before and after constriction. At certain time,
% the threshold to which the division protein accumulates changes.

clear

cells=[];
% Input parameters
v_inp=1.4; % initial size
tau= 72.5; % generation time in minutes
gr1= 0.0087; % growth rate 1 in min^-1
gr2 = 0.0126; % growth rate 2 in min^-1
ngen=24; % number of generation
bii = 0.6; % length per origin between initiations (um)
mtcd = 18.3; % time from constriction to division (min) 
mtic = 75.6; % time between DNA replication initiation and constriction (min)
s0 = 1; % production rate of FtsZ  
k = 0.029; % maturation rate of the reporter protein min^-1 
mthresh = 286.2069; % Threshold at constriction for division protein
mthresh_orig = mthresh;
kp_cells= 5*tau; % Time (min) for reaching a population steady state 
ch_cd= 12*tau; % the bbd changes to dbbd um at this time in the simulations
dthresh = 0.07; % percentage change in threshold level 
time_dly = 37; % time delay before FtsZ overproduction

%-----Noise parameters
cvl1=0.2*gr1; cvl2=0.2*gr2; % std dev in growth rate
cvt = 0.1*bii; % std dev in initiation size added
sigtic = 0.2*mtic; % std dev of the time elapsed between DNA replication initiation and constriction
sigtcd = 0.2*mtcd; % std dev of the time elapsed between constriction and cell division
cvt2= 0.2; % CV of production rate of FtsZ
cvtthr = 0.2*mthresh; % std dev of the threshold

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
rep_pop = NaN(1,6000); % replication limiting
tot_pop = NaN(1,6000); % tot cells
times=NaN(1,6000); % timing of the simulation

lmin=gr1/(5);

for no_in=1:100
x= Cell_const(v_inp);
x.rate= gr1;

x.nthresh = mthresh + randn()*cvtthr;
if x.nthresh < x.n
    x.nthresh = mthresh + randn()*cvtthr;
end

%-----For keeping track of events

ii_in = bii+ randn()*cvt;
while ii_in<0
    ii_in=bii+ randn()*cvt;
end
x.vNextInit = x.v + x.oris*(ii_in);

x.tNextConst_in(1)= mtic+randn()*sigtic;
while x.tNextConst_in(1)<0
    x.tNextConst_in(1)= mtic+randn()*sigtic;
end

x.vNextConst = 0;

s_in=s0+randn()*cvt2*s0;
while s_in<0
    s_in=s0+randn()*cvt2*s0;
end
x.s= s_in;

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

    % Threshold changes after a certain time
     if simTime>ch_cd+time_dly 
        mthresh_prev = mthresh/mthresh_orig;
        mthresh=mthresh_orig*(1+dthresh*(1-exp(-(simTime-(ch_cd+time_dly))/tau))); 
        for ct=1:length(cells)
            x=cells(ct);
            x.nthresh = x.nthresh*mthresh/(mthresh_prev*mthresh_orig);
        end
        post_ch= 1;
    end
    
    for ct=1:length(cells)
    x=cells(ct);
    x.t = x.t + tStep;
    x.v = x.v*exp(tStep*x.rate); % grow exponentially
    x.n = x.n+x.s*x.v*tStep;
    x.nm= x.nm+k*(x.n-x.nm)*tStep;
    x.tNextConst_in=x.tNextConst_in-tStep;
    x.tNextDivs = x.tNextDivs-tStep;
    end
    %-----Perform events
    
    cellsNew=[];

    for ct=1:length(cells)
    x=cells(ct);
    
    %-----Constriction
    
    if x.vNextConst==0 && x.n-x.nthresh>=0 && ~isempty(x.tNextConst_in) && x.tNextConst_in(1)<=0

        if x.tNextConst_in(1)>-tStep
            rep = rep+1;
            x.rep_cell = 1;
            x.tot_cell = 1;
        else
            x.rep_cell = 0;
            x.tot_cell = 1;
        end
        tot = tot+1;

        x.vOfInits(1)=[];
        x.tOfInits(1)=[];
        x.oOfInits(1)=[];

        x.tNextConst_in(1)=[];
        x.vNextConst = 1;

        x.vOfConst = [x.vOfConst [x.v; x.t-x.tLastDiv; x.nm/x.v; x.n/x.v]];

        x.tNextDivs= mtcd+randn()*sigtcd;
        while x.tNextDivs<0
            x.tNextDivs= mtcd+randn()*sigtcd;
        end

        x.rate = max(gr2 + randn()*cvl2,lmin);

    end

    %-----Division
    
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
        x.c_tc_cell = x.nthresh;
        
        %-----Update cell
        
        x.tLastDiv = x.t;
        x.vd = x.v;
        rc = 0.5+randn()*0.03;
        x.v = x.vd*rc;
        x.vb = x.v;
        x.n = x.n*0.5;
        x.nm= x.nm*0.5;
        x.rec_data=1;
        x.rate = max(gr1 + randn()*cvl1,lmin);
        x.oris = x.oris/2;
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
        x.vNextConst = 0;

        %----------Update vNextInit and vNextDiv (see Methods)
        x.vNextInit = x.vNextInit/2 + (rc-0.5)*x.vi;
        
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
        y.vNextInit = y.vNextInit + (1-2*rc)*y.vi;
        cellsNew =[cellsNew y];
        end
    end

    %-----Initiate
    %Upon accumulating enough volume
    if x.v - x.vNextInit >= 0
        %-----Record event
        x.tLastInit = x.t;  %%% why duplicated with x.tOfInits(-1)?
        x.vi = x.v;
        x.vOfInits= [x.vOfInits x.v];
        x.tOfInits= [x.tOfInits x.t];
        x.oOfInits= [x.oOfInits 2*x.oris];
        
        %-----Update cell
        x.oris = 2*x.oris; %double number of origins
        ii_in = bii+ randn()*cvt;
        while ii_in<0
            ii_in=bii+ randn()*cvt;
        end
        x.vNextInit = x.v + x.oris*(ii_in);

        xtic=mtic+randn()*sigtic;
        while xtic<0
            xtic=mtic+randn()*sigtic;
        end
        x.tNextConst_in=[x.tNextConst_in xtic];
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

%%

% Plotting 

figure
x= times; y= tc_pop;
[bin,da,yfit, P, errbr]=binning_with_error_1(y,x,-450:30:800);
hold on
s=scatter(bin,da, 100, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p=plot(bin,da, 'LineWidth',2, 'Color', 'k');
p=xline(0, 'LineWidth',2, 'Color', 'k', 'LineStyle','--');
eb = errorbar(bin,da,errbr, 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k');
set(get(get(eb,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('T_d (min)');
xlabel('Time (min)');
box on
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
xlim([-400, 800])

figure
x= times; y= ld_pop;
[bin,da,yfit, P, errbr]=binning_with_error_1(y,x,-450:30:800);
hold on
s=scatter(bin,da, 100, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p=plot(bin,da, 'LineWidth',2, 'Color', 'k');
p=xline(0, 'LineWidth',2, 'Color', 'k', 'LineStyle','--');
eb = errorbar(bin,da,errbr, 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k');
set(get(get(eb,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('L_d (\mum)');
xlabel('Time (min)');
box on
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
xlim([-400, 800])

figure
x= times-td_pop+tc_pop; y= lc_pop;
[bin,da,yfit, P, errbr]=binning_with_error_1(y,x,-450:30:800);
hold on
s=scatter(bin,da, 100, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p=plot(bin,da, 'LineWidth',2, 'Color', 'k');
p=xline(0, 'LineWidth',2, 'Color', 'k', 'LineStyle','--');
eb = errorbar(bin,da,errbr, 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k');
set(get(get(eb,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel('L_c (\mum)');
xlabel('Time (min)');
box on
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
xlim([-400, 800])

figure
x= times-td_pop+tc_pop; y= rep_pop;
[bin1,da1,yfit, P, errbr]=binning_with_error_1(y,x,-450:30:800);
x= times-td_pop+tc_pop; y= tot_pop;
[bin2,da2,yfit, P, errbr]=binning_with_error_1(y,x,-450:30:800);
hold on
s=scatter(bin1,1-da1./da2, 100, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p=plot(bin1,1-da1./da2, 'LineWidth',2, 'Color', 'k');
p=xline(0, 'LineWidth',2, 'Color', 'k', 'LineStyle','--');
ylabel('Prob. FtsZ limiting');
xlabel('Time (min)');
box on
set(gca, 'FontSize', 30)
set(gcf, 'Position',[276,42,777,602])
xlim([-400, 800])
ylim([0,1])

