function d=finddimension(data)
affinity = CalculateAffinity(data);
% compute the degree matrix
for i=1:size(affinity,1)
    D(i,i) = sum(affinity(i,:));
end
%%Define initial Pij based on covarince between data
for i=1:size(data,1) 
    for j=1:size(data,1) 
      p=corrcoef(data(i,:),data(j,:));
      Pij(i,j)=p(1,2);
    end
end
for i=1:size(data,1) 
    for j=1:size(data,1) 
      Pij(i,j)=abs(Pij(i,j))/abs(sum(Pij(i,:)));
    end
end

%%%Compute Derivation of laplacian with

 for i=1:size(data,1) 
     for j=1:size(data,1)    
        Lij(i,j)=(-1/(D(i,i)*D(j,j)))*((affinity(i,j)/(D(i,i).^2))*sum((data(i,:)-data(j,:))*(data(i,:)-data(j,:))'*affinity(i,j),1)+(affinity(i,j)/(D(j,j).^2))*sum((data(i,:)-data(j,:))*(data(i,:)-data(j,:))'*affinity(i,j),2)-affinity(i,j)*(data(i,:)-data(j,:))*(data(i,:)-data(j,:))');
    end   
 end

Lij1=Lij;
Lij=Lij.*Pij;
%Update Pij befor compute eignvector
af=zeros(size(data,1),size(data,1)); 
landa=.2;beta=0.1;miu=.002;%h is coeffient of entropy
epsilon=.5;
while (epsilon<0.01)
    for i=1:N
        for j=i:N
             Pij1(i,j)=Pij(i,j)-miu*(-(beta-landa)*lambertw((-Lij(i,j)*exp(-beta/(beta-landa)))/(beta-landa))/Lij(i,j));           
             affinity(i,j)=Lij1(i,j).* Pij1(i,j);             
        end
    end     
       epsilon=Pij1(i,j)-Pij(i,j);
       Pij(i,j)=Pij1(i,j);
       af=affinity;
end
% compute the unnormalized graph laplacian matrix
L = D - af; 
[eigVectors,eigValues] = eigs(L,10,'smallestabs');

maxDim=10;
for i = 2:maxDim
    if (eigValues(i-1,i-1) > 0 && eigValues(i,i) > 0)
        gap(i) = abs(eigValues(i-1,i-1) - eigValues(i,i));
    else
        gap(i) = 0;
    end
end
[~,Dim] = max(gap(1:maxDim));
d = min(Dim, maxDim);
