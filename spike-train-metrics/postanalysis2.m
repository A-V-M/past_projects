clc;
TAB = char(9);
disp(sprintf('\n'));
disp(sprintf('\t [1] Metric progression over time'));
disp(sprintf('\t [2] Metric comparison (for each pair-wise combination of observed units) \n \t within specific time frame'));
disp(sprintf('\t [3] Spike count progression over time'));
disp(sprintf('\n'));
choice = input('Choice: ');
switch choice
    case 1
        n = 1;
        c = 1;
while n <= combinations
    hold on
    figure(ceil(n/10));   
    subplot(10,1,c);
    c = c + 1;    
    bar(cumsum(t),metrics(:,n));
    s = sprintf('%s / %s ',spiketrainData{spiketrainCombinations(n,1)},spiketrainData{spiketrainCombinations(n,2)});
    title(s);
    if c > 10;c = 1;xlabel('Time (s)','fontweight','b');end
    hold off
    n = n + 1;
end
    case 2
        
        l = input('Enter lower time limit [in seconds from origin]: ');
        u = input('Enter upper time limit [in seconds from origin]: ');
        disczero = input('Discard zero values for distance matrix? ', 's');
        
for a = 1:combinations
    %have to insert check here for values lower than temporal resolution
    lower = find(round(cumsum(t)) == l);
    if u == 9999; endtime = cumsum(t);u=round(endtime(end));end
    upper = find(round(cumsum(t)) == u);
    if l == 0;timeFrame = 1:upper;else timeFrame = (lower+1):upper;end
    i = find(metrics(timeFrame,a)>0); % zero values discarded
    if disczero == 'y';avg_metrics(a) = mean(metrics(timeFrame(i),a));else avg_metrics(a) = mean(metrics(timeFrame,a));end
end
    
for j=1:length(spiketrainData)
[rows, columns] = find(spiketrainCombinations == j);
rows = sort(rows);
n = 1;
for i = 1:length(spiketrainData)
if j == i
d_matrix1(j, i) = 0;
else
d_matrix1(j, i) = avg_metrics(rows(n));
n = n + 1;
end
end
end
figure;bar(d_matrix1,'stack');legend(spiketrainData,'Location','NorthEastOutside');

disp(sprintf('\n\t Data report for pair-wise metrics for interval %ds to %ds [%0.5fs] \n', l,u,sum(t(timeFrame))));
disp(sprintf('\t %s', spiketrainData{1:observedUnits}));
disp('Mean:');disp(sprintf('\t %0.3f',mean(d_matrix1)));
disp('Std:');disp(sprintf('\t %0.3f',std(d_matrix1)));
disp('Variance:');disp(sprintf('\t %0.3f',var(d_matrix1)));
disp('Spike rate (spikes/s):');disp(sprintf('\t %0.3f',sum(spikeCount(timeFrame,:))/sum(t(timeFrame))));
for c = 1:observedUnits
    disp(spiketrainData{c});
    disp(sprintf('\t %0.3f',d_matrix1(:,c)));
end
disp(sprintf('\n'));

    case 3
        figure;
        for a = 1:observedUnits
            subplot(observedUnits,1,a);
            bar(cumsum(t),spikeCount(1:(end-1),a));
            s = spiketrainData{a};
            title(s);
        end

    otherwise
        disp('Invalid option');
        break;
end
