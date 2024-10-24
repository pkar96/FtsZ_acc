function [bin,da, yfit, P, err]= binning_with_error_1(y,x, tot_bin)
[N, edges]= histcounts(x,tot_bin);
bin=zeros(1,length(edges)-1);
da=zeros(1,length(bin));
err=zeros(1,length(bin));
for i=1:length(bin)
    da(i)=mean(y(edges(i)<=x & edges(i+1)>=x),"omitnan");
    err(i)=std(y(edges(i)<=x & edges(i+1)>=x))/sqrt(length(y(edges(i)<=x & edges(i+1)>=x)));
    bin(i)=(edges(i)+edges(i+1))/2;
end
rem= find(N<10);
if ~isempty(rem)
    da(rem)=[];
    err(rem)=[];
    bin(rem)=[];
end
P = polyfit(x,y,1);
yfit = P(1)*x+P(2);
end