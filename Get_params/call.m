% find the best parameters corresponding to the threshold level of division
% protein at constriction and the time between DNA replication initiation
% and division
clear

% Parameters to be compared against in the experiment. Corresponds to 
% [Lc,bf   Lc,af]
mean_exp = [2.5564, 2.8623];

% Parameters cycled through- sbc is related to threshold level and fac is
% related to the time between initiation and constriction. See later for a
% conversion to protein number and time in mins
sbc = 0.65:0.01:0.77; 
fac= 0.7:0.01:0.74; 

% Create a grid
V={sbc, fac};
C = cell(1,numel(V));
[C{:}] = ndgrid(V{:});
C = cellfun(@(X) reshape(X,[],1),C,'UniformOutput',false);
C = horzcat(C{:});

dist=zeros(1, length(C));

% Run the simulations in parallel and calculate the loss function
parfor i=1:length(dist)
     [dist(i)]=concurrent_prod_rate_change_timer_from_init(C(i,:),mean_exp,i);
end

% Min residual values
[M,ind] = min(dist);
[C(ind,1), C(ind,2)]

% Plot
[X,Y] = meshgrid((log(3.6/(2*0.6))*72.5/log(2)-21.3)*fac,sbc+1.8);
Z = reshape(dist,length(sbc),length(fac));

figure
s= surf(X, Y, Z, 'FaceAlpha',0.5,'EdgeColor','none');

%%

% find the best parameters corresponding to the threshold level of the
% other division protein at constriction 

clear

% Parameters to be compared against in the experiment. Corresponds to 
% [Lc,bf   Lc,af]
mean_exp = [2.5564, 2.8623];

% Parameters cycled through- sbc is related to threshold level of the other
% division protein keeping the division protein of FtsZ fixed (determined
% by the simulation in the above section). See later for a conversion to 
% protein number
sbc = 0.34:0.01:0.44;

dist=zeros(1, length(sbc));

% Run the simulations in parallel and calculate the loss function
parfor i=1:length(dist)
     [dist(i)]=concurrent_prod_rate_change_2_div_protein(sbc(i),mean_exp,i);
end

% Min residual values
[M,ind] = min(dist);
[sbc(ind)]

figure
s= plot(sbc+1.8, dist);

