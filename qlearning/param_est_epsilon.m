%this code is to estimate a value for epsilon given a set of payoff and
%safety setting values. first i will try to produce a histogram for each
%epsilon value at each possible safety setting. then i will try to extract 
%a certain probability (i.e. the likelihood) of attaining a certain payoff
%value pi given a greediness setting, epsilon -> p(pi|epsilon).

load likelihood_variables_epsilon.mat

bins = -50:0.05:1;
values_hist = zeros(10,length(sE),length(bins));

payoff = 0.7;%t_cost_lambda;
S = 1; %<-safety setting with S being (S-1)*10) - (S*10 ).



for n=1:10 %safety settings
    for nE=1:length(sE)
        values = reshape(t_cost(nE,range_values_indices{n,nE},:),1,[]);
        values_hist(n,nE,:) = histc(values,bins);
        values_hist(n,nE,:) = values_hist(n,nE,:)/sum(values_hist(n,nE,:));%normalise to extract probability [0 1].
        v = reshape(values_hist(n,nE,:),1,[]);
        
        if payoff<0
            index0=find(bins<=0);
            P0(n,nE)=sum(v(index0));
        end
            
        index=find(bins==payoff);
        
        if isempty(index)==1 
        
            tmp = [bins payoff];
            tmp=sort(tmp);
            index=find(tmp==payoff);
            index = index-1;
        
        end
        P(n,nE)=v(index);
        %figure(1);bar(bins,v);
    end
end

max_likelihood = max(P(S,:));
if max_likelihood == 0; ML = 0; disp('likelihood is zero for all possible epsilon at this safety setting...Posterior probs. are for negative payoffs overall'); P = P0;
else
    
index=find(P(S,:)==max(P(S,:)));

ML = P(S,index)

sE(index)

end


P_epsilon = 1/length(sE); % 33 possible values for epsilon so the prior or 'naive' probability (P(epsilon)) should be 1/33.

P_payoff = sum(P(S,:) * P_epsilon); %p(pi) = SUM[(p(pi|epsilon(n)) x p(epsilon(n)))] (p(epsilon) is the same for ALL possible values.)

for n=1:length(sE);
    PP(n)=(P(S,n)*(P_epsilon))/P_payoff; %p(epsilon|pi) = [p(pi|epsilon) x p(epsilon)]/p(pi)
end





        




