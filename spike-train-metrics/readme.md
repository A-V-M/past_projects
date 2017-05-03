# spike-train-metrics
online computation of van Rossum metric from incoming spike trains


At a fundamental level it is a widely accepted belief that neurons transmit information in discrete all-or-nothing events 
elicited through a rapid influx of positively charged ions which is then followed by a refractory period. 
These events can be sequenced in a series over a relatively long period of time – these temporal series of neural events 
are termed as spike trains. These spike trains can be analysed on the basis of their temporal 
arrangement whereby different mathematical methods can be implemented for the extraction of information that might be encoded. 

an interface was written within the MATLAB environment to analyse spike train data. 
Raw neurophysiological data was recorded from rat prefrontral cortex, via a standard tetrode apparatus, and acquired from 
within a computer interface via CHEETAH® Neuralynx for the purposes of this project. 

Spike sorting was applied and the extracted spike train data was thereby used for analysis. The purpose of the software 
in question was to provide for a real-time visualisation of the variability and hence the dissimilarity between pairs of
recorded units. Metrics were computed via the van Rossum metric with software written by Dr. Thomas Wennekers. 
