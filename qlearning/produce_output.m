%produce output
%1) correlation - prob. of taking right turn vs. average trap risk
%2) percentage of maximum payoff

%1.

correlRiskVSProbRight = zeros(35,196);
Risk = cell(35,196);
ProbRight_ = cell(35,196);

d_corr=zeros(35,196);
d = zeros(35,196);
d_eucl=zeros(35,196);

for ep=1:35
    for lam=1:196
    
    c = zeros(250,1);
    
    for n=1:250; 
        c(n) = prob_right{ep,lam,n}(4);
    end
    
    ProbSafe = AverageTrapCount(lam,ep);
    
    Risk{ep,lam} = [1 ProbSafe 1 ProbSafe 1];
    
    ProbRight_{ep,lam} = [mean(c(1:50)),mean(c(51:100)),mean(c(101:150)),mean(c(151:200)),mean(c(201:250))];
    
    var1 = Risk{ep,lam}';
    
    var2 = ProbRight_{ep,lam}';
    
    %[c p] = corr(var1,var2);
    
    d_corr(ep,lam) = pdist([var1';var2'],'corr');
    
    d(ep,lam) = acos(dot(var1,var2)/(norm(var1)*norm(var2)));
    
    d_eucl(ep,lam) = pdist([var1';var2'],'euclidean');
    
    %correlRiskVSProbRight(ep,lam) = c;

    end
end

figure;subplot(3,1,1);imagesc(nL,sE,d);
hold on;subplot(3,1,2);imagesc(nL,sE,d_eucl);
subplot(3,1,3);imagesc(nL,sE,d_corr);



%2.
t_cost_risky1 = zeros(35,196);
t_cost_risky2 = zeros(35,196);
t_cost_risky3 = zeros(35,196);
t_cost_risky4 = zeros(35,196);
t_cost_risky5 = zeros(35,196);

t_cost_per_cycle = zeros(35,196,5);

perc_max_payoff = cell(35,196);

for ep=1:35
    for lam=1:196
        
        p = AverageTrapCount(lam,ep);
                
        max_payoff = (p * 0.9) + ((1-p)*0.78);

        t_cost_per_cycle(ep,lam,1) = mean(t_cost(ep,lam,1:50))/0.9;
        t_cost_per_cycle(ep,lam,2) = mean(t_cost(ep,lam,51:100))/max_payoff;
        t_cost_per_cycle(ep,lam,3) = mean(t_cost(ep,lam,101:150))/0.9;
        t_cost_per_cycle(ep,lam,4) = mean(t_cost(ep,lam,151:200))/max_payoff;
        t_cost_per_cycle(ep,lam,5) = mean(t_cost(ep,lam,201:250))/0.9;
        
        perc_max_payoff{ep,lam} = [t_cost_risky1(ep,lam)/0.9 t_cost_risky2(ep,lam)/max_payoff t_cost_risky3(ep,lam)/0.9 t_cost_risky4(ep,lam)/max_payoff t_cost_risky5(ep,lam)/0.9];
        
        perc_max_payoff{ep,lam}(perc_max_payoff{ep,lam}<0) = 0;

    end
end

t_cost_per_cycle(t_cost_per_cycle<0) = 0;

figure;
subplot(3,2,1);imagesc(nL,sE,t_cost_per_cycle(:,:,1)*100,[0 100]);title('Episodes 1-50: Safe','FontWeight','bold');
hold on;
subplot(3,2,2);imagesc(nL,sE,t_cost_per_cycle(:,:,2)*100,[0 100]);title('Episodes 51-100: Trap set','FontWeight','bold');
subplot(3,2,3);imagesc(nL,sE,t_cost_per_cycle(:,:,3)*100,[0 100]);title('Episodes 101-150: Safe','FontWeight','bold');
subplot(3,2,4);imagesc(nL,sE,t_cost_per_cycle(:,:,4)*100,[0 100]);title('Episodes 151-200: Trap set','FontWeight','bold');
subplot(3,2,5);imagesc(nL,sE,t_cost_per_cycle(:,:,5)*100,[0 100]);title('Episodes 201-250: Safe','FontWeight','bold');

        
        
        
        