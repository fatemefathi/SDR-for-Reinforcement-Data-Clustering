 clc;
 close all;
 clear all;

%load data and label and numcluster



nn=20;
accarrours=zeros(nn,1);
NMIours=zeros(nn,1);
accarrours2=zeros(nn,1);
NMIours2=zeros(nn,1);
for numiter=1:nn 
   numiter
clu2ours=RLClustering(X',c);
result_our = ClusteringMeasure(label, clu2ours);
Accours1=result_our(1);
Accours2=result_our(2);
accarrours(numiter)=Accours1;
NMIours(numiter)=Accours2;
end
meanaccours=mean(accarrours);
stdaccours=std(accarrours);
meanNMI=mean(NMIours);
stdNMI=std(NMIours);
fprintf(' MeanAccuracyours: %.4f + stdAccuracyours: %.4f \n',meanaccours,stdaccours);
fprintf(' meanNMI: %.4f + stdNMI: %.4f \n',meanNMI,stdNMI);


