%AM, CTCN/2009
%client code

%connect net-com

serverName = 'localhost';
disp(sprintf('Connecting to %s...', serverName));
succeeded = NlxConnectToServer(serverName);
if succeeded ~= 1
    disp(sprintf('FAILED connect to %s. Exiting script.', serverName));
    return;
else
    disp(sprintf('Connected to %s.', serverName));
end

%Identify this program to the server we're connected to.
succeeded = NlxSetApplicationName('andreas');
if succeeded ~= 1
    disp 'FAILED set the application name'
else
    disp 'PASSED set the application name'
end

[succeeded, cheetahObjects, cheetahTypes] = NlxGetCheetahObjectsAndTypes;
if succeeded == 0
    disp 'FAILED get cheetah objects and types'
else
    disp 'PASSED get cheetah objects and types'
end

%open up a stream for all objects
for index = 1:length(cheetahObjects)
    succeeded = NlxOpenStream(cheetahObjects(index));
    if succeeded == 0
        disp(sprintf('FAILED to open stream for %s', char(cheetahObjects(index))));
        break;
    end
end;
if succeeded == 1
    disp 'PASSED open stream for all current objects'
end

%objects should be post-spike sorting phase - so each object = neuron/unit

[succeeded, cheetahReply] = NlxSendCommand('-PostEvent "Test Event" 10 11');
if succeeded == 0
    disp 'FAILED to send command'
else
    disp 'PASSED send command'
end



inp_option = input('File or console [c]? (default is file) ','s');
if inp_option == 'c'
    objects = input('Tetrode objects: '); % custom setting for what objects are to be streamed 
observedUnits = 0;
spiketrIndex = 0;
i = 1;
for n = 1:length(objects) 
    objectcellID{n} = input('Unit cluster IDs: ');   
    [r, neurons] = size(objectcellID{n});
    observedUnits = observedUnits + neurons;
    while i <= neurons % loop to construct cell array which contains descriptions for observed neuron
        spiketrainData{spiketrIndex + i} = [char(cheetahObjects(objects(n))),'-', num2str(objectcellID{n}(i))]; 
        i = i + 1;
    end
    spiketrIndex = spiketrIndex + neurons; % set index for next object
    i = 1; %reset counter
end
else
    disp('Loading data file...');
fid = fopen('nlxclient.cfg');
objects = str2num(fgetl(fid));
observedUnits = 0;
spiketrIndex = 0;
i = 1;
for n = 1:length(objects) 
    objectcellID{n} = str2num(fgetl(fid));   
    [r, neurons] = size(objectcellID{n});
    observedUnits = observedUnits + neurons;
    while i <= neurons % loop to construct cell array which contains descriptions for observed neuron
        spiketrainData{spiketrIndex + i} = [char(cheetahObjects(objects(n))),'-', num2str(objectcellID{n}(i))]; 
        i = i + 1;
    end
    spiketrIndex = spiketrIndex + neurons; % set index for next object
    i = 1; %reset counter
end
disp(sprintf('PASSED observed unit numbers acquired [%d total]', observedUnits));
end
lambda = input('Lambda value: ');
disp(sprintf('TIME set at %s', datestr(now)));

%initialise variables
spiketrainCombinations = combnk((1:observedUnits),2); % calculate all possible pair-wise combinations of spike trains
[combinations, columns] = size(spiketrainCombinations);
NoMoreData = 0;
recs_returned = zeros(length(cheetahObjects));
streamPass = 1; % while loop counter

z_(1:combinations,1) = 0;
z_(1:combinations,2) = 0;
z_(1:combinations,3) = -1000;
spike_mm = zeros(9,2);
last_spike = zeros(9,1);
spCount_ = zeros(9,1);
t_zero = 0; % check if first block of spikes
t = 0;
t_int = 1;
k  = -1;
DataInThePipeline = 0;    
while NlxAreWeConnected == 1
    tic;
    
    spiketrainIndex = 1;
    
    c = 0;
    for objectIndex = 1:length(objects)
            objectToRetrieve = char(cheetahObjects(objects(objectIndex)));
            [succeeded, dataArray, timeStampArray, spikeChannelNumberArray, cellNumberArray, featureArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewTTData(objectToRetrieve);
            if succeeded == 0
                 disp(sprintf('FAILED to get new data for stream %s on pass %d', objectToRetrieve, pass));
                 NoMoreData = 1;
                 break;
            else
                disp(sprintf('[%s] Retrieved %d records for %s with %d dropped.', datestr(now), numRecordsReturned, objectToRetrieve, numRecordsDropped));
                for spI = 1:length(objectcellID{objectIndex})
                spiketr{spiketrainIndex} = timeStampArray(find(cellNumberArray==objectcellID{objectIndex}(spI))); 
                spikeCount(streamPass,spiketrainIndex) = length(spiketr{spiketrainIndex});
                spike_train{spiketrainIndex} = (double(spiketr{spiketrainIndex})/1000000);               
                         
                         if k > 0
                         s_ = spike_train{spiketrainIndex};
                         b = t_zero + (k * t_int);
                         s2{spiketrainIndex} = s_(find(s_<b));
                         buffer_ = s_buffer{spiketrainIndex};
                         appendum = buffer_(find(buffer_<b & buffer_>(b-t_int)));
                         st_ = [s1{spiketrainIndex} s2{spiketrainIndex} appendum];
                         st{spiketrainIndex} = sort(st_);
                        
                         s1{spiketrainIndex} = s_(find(s_>b & s_<(b+t_int)));
                         buffer_1 =  s_(find(s_>(b+t_int)));
                         buffer_2 =  buffer_(find(s_buffer{spiketrainIndex}>(b+t_int)));
                         s_buffer{spiketrainIndex} = [buffer_1 buffer_2];
                         
                if spiketrainIndex == 7 && DataInThePipeline == 1 
                    obs_n{streamPass} = st{spiketrainIndex};
                end
                if spiketrainIndex == 2 && DataInThePipeline == 1
                    obs_n2{streamPass} = st{spiketrainIndex};
                end
                         end       
                spiketrainIndex = spiketrainIndex + 1; 
                
                end
                
                
                recs_returned(objectIndex) = numRecordsReturned;
                if objectIndex == length(objects) %check if at final iteration            
                   if find(recs_returned>0);
                       DataInThePipeline = 1;
                       if streamPass == 1
                           spCount_ = spikeCount;
                       else
                           spCount_ = sum(spikeCount);
                       end
                       
                       if sum(spikeCount(streamPass,:)) > 0
                                        
                           if t_zero == 0 % compute tzero
                            for p=1:observedUnits
                             s_ = spike_train{p}; 
                             if isempty(s_) == 0
                                spike_mm(p,1) = (max(s_));
                                spike_mm(p,2) = (min(s_));                   
                             end                             
                            end
                            t_max = max(spike_mm(:,1));
                            t_min = min(spike_mm(find(spike_mm > 0)));
                            t_zero = (t_min + t_max)/2;
                            k = 0;
                            for p=1:observedUnits   
                                s_ = spike_train{p};
                                s1{p} = s_(find(s_>t_zero & s_<(t_zero + t_int)));
                                s_buffer{p} = s_(find(s_>(t_zero + t_int)));                       
                            end
                           end
                                              
                       end
                       
                       for i=1:combinations 
                         n1 = spiketrainCombinations(i,1);
                         n2 = spiketrainCombinations(i,2);
                         
                         if k > 0
                         x = st{n1}; %x = (double(x)/1000000);  %sorting may be unnecessary with proper data
                         y = st{n2}; %y = (double(y)/1000000); 
                         %D = mxnorm2(x,y,lambda);
                         %metrics(streamPass,i) = sqrt(D);
          
                         
                                                
                         z = z_(i,:);
                         mxnorm2b(st{n1},st{n2},lambda,b,z);
                         z_(i,:) = z;
                         
                         metrics_d(streamPass,i) = sqrt(z_(i,1));                         
                         spDist(streamPass,i) = (spCount_(n1) - spCount_(n2))^2;
                         
                            
                         end
                                                           
                       end
                       
                       if k > 0;
                       s_metrics = cmdscale(metrics_d(streamPass,:));
                       s_dims1 = size(s_metrics);
                       if s_dims1(2) >= 2
                       figure(1);subplot(2,1,1);plot(s_metrics(:,1),s_metrics(:,2),'o');text(s_metrics(:,1),s_metrics(:,2),spiketrainData);
                       else
                           figure(1);subplot(2,1,1);plot(s_metrics, 1:s_dims1(1)', 'o');text(s_metrics,1:s_dims1(1)',spiketrainData);
                       end
                       
                       s_count = cmdscale(spDist(streamPass,:));
                       s_dims2 = size(s_count);
                       if s_dims2(2) >= 2
                       figure(1);subplot(2,1,2);plot(s_count(:,1),s_count(:,2),'o');text(s_count(:,1),s_count(:,2),spiketrainData);
                       else
                           figure(1);subplot(2,1,2);plot(s_count, 1:s_dims2(1)', 'o');text(s_count,1:s_dims2(1)',spiketrainData);
                       end
                       end
                   else
                       if DataInThePipeline == 0
                           disp('Waiting for data...');
                           NoMoreData = 0;
                           streamPass = 1;
                       else
                       NoMoreData = 1;% this is to end data streaming when then there is no data in the stream
                       end
                   end
                end
            end
    end
       if NoMoreData == 1
           break;
       end
       
       streamPass = streamPass + 1;
       if DataInThePipeline == 1; proc_time = toc; pause(t_int-proc_time);t(streamPass) = streamPass * t_int; end
       if k>=0; k = k + 1;end
       
end %<-- end of while loop

for index = 1:length(cheetahObjects)
    succeeded = NlxCloseStream(cheetahObjects(index));
    if succeeded == 0
        disp(sprintf('FAILED to close stream for %s', char(cheetahObjects(index))));
        break;
    end
end;

if succeeded == 1
    disp 'PASSED close stream for all current objects'
end


%Disconnects from the server and shuts down NetCom
succeeded = NlxDisconnectFromServer();
if succeeded ~= 1
    disp 'FAILED disconnect from server'
else
    disp 'PASSED disconnect from server'
end

%following code is for some post-analysis of metrics as a function of time. still in developing stages
%

%see postanalysis.m
