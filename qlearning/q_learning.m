
sA=[exp(-0.075*(0:1:37))*25 0.5 0];

nL = [0.2:0.01:1 1.2:0.2:24];

tEpisodes = 250;

AverageTrapCount = zeros(length(nL),length(sA));

action_taken = cell(length(sA),length(nL),tEpisodes);

prob_right = cell(length(sA),length(nL),tEpisodes);

t_cost = zeros(length(sA),length(nL),tEpisodes);

t_cost_lambda = zeros(length(sA),length(nL));

prediction_error = zeros(length(sA),length(nL),tEpisodes);

for nA=1:length(sA)
disp(nA);
tic;   
for nlambda = 1:length(nL)

lambda = nL(nlambda);

t_int = poissrnd(lambda,10000,1);

t_int(t_int == 0) = [];

t_timings = cumsum(t_int);

T = zeros(tEpisodes,1);

TrapCount = 0;
nTEpisodes = 0;

for n = 2:2:floor(tEpisodes/50)
 nTEpisodes = nTEpisodes + 1;
 lower = ((n-1) * 50)+1;
 upper = n*50;
 T(t_timings(lower<=t_timings & t_timings<=upper)) = -0.5;
 TrapCount = TrapCount + length(find(T(lower:upper)==0))/length(T(lower:upper));
 
end

AverageTrapCount(nlambda,nA) = TrapCount/nTEpisodes;

gamma = 0.8;

B = 0;
Q = ones(42,4)*B;

Probs = zeros(42,4);

c = 0;   

k = -0.005; %decay constant for exponential decay (epsilon)

N = [-1 0];
S = [1 0];
E = [0 1];
W = [0 -1];

action_set = {N W S E};

alpha = 1;
%epsilon = sE(nE);
alpha_SM = sA(nA);

%t_cost = zeros(tEpisodes,1);

for episode=1:tEpisodes
%sprintf('Now in episode %d ... \n', episode)

%reward matrix is time-dependent

state = [6 5];
Term = 'N';

c = c + 1;

P = T(c);

%env_model;%Rmap(5,6)

while Term == 'N';

%select action 

%epsilon = exp(k * episode); %decay function for epsilon

state_i = map_index(state(1),state(2));
%sprintf('[%d] state: %d \n', episode, state_i)
available_actions = find(TM(state_i,:));

Q_Max = max(Q(state_i,available_actions));

max_action = find(Q(state_i,available_actions) == Q_Max);

%softmax action selection

x = zeros(1,4); 

for nAction = 1:length(available_actions)
    
x(available_actions(nAction)) = exp(alpha_SM * Q(state_i, available_actions(nAction)));

end

x=x./sum(x);

%Probs(state_i,:) = x(:); 

action = randsample(1:length(action_set),1,true,x);

state1 = state + action_set{action};
state1_i = map_index(state1(1),state1(2));



%sprintf('[%d] taking action: %d \n', episode, action)

%update function
q_old = Q(state_i, action);
Q(state_i, action) = Q(state_i, action) + alpha * (R(state_i,action) + gamma * max(Q(state1_i,:)) - Q(state_i,action));

%chronology(episode,c) = state_i;

if [state_i action] == [29 4]
prediction_error(nA,nlambda,episode) = Q(state_i,action)-q_old;
end


if state_i == 29
    action_taken{nA,nlambda,episode}(end+1) = action_set{action}(1,2);
    prob_right{nA,nlambda,episode} = x; 
end



%t_cost(episode) = t_cost(episode) + R(state_i,action);
t_cost(nA,nlambda,episode) = t_cost(nA,nlambda,episode) + R(state_i,action);

state = state1;
mapD = map;
mapD(state(1),state(2)) = 5; mapD(5,6) = P + 1;
%figure(1);imagesc(mapD);
%figure(2);imagesc(Q);
if state == [1 6]; Term = 'Y';end

end

%Qlog(episode,:) = Q(29,:);

end

%t_cost_lambda(nlambda,nA) = mean(t_cost)';
t_cost_lambda(nlambda,nA) = mean(t_cost(nA,nlambda,:))';


end
duration_per_full_cycle(nA)=toc;

end
