function [clu2]=RLClustering(reduceddata,numclu)

P=reduceddata;
lb=min(P);
ub=max(P);
ncluster=numclu;
np=size(P,1);%number of points
dp=size(P,2);%dimension of points
center=zeros(numclu,dp);
newcenter=zeros(numclu,dp);
%%%%%%%%%%%%%%%%%%%%%%%%%%Q-learning%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InitialEpsilon = 0.9;
EpsilonDamp = 0.995;
alpha = 0.1; % Learning Rate
alphaDamp=0.995;
gamma = 0.9; % Discount Rate
RC=ncluster;% Number Of States
nA =3 ; % Number Of Actions(add/remove/no change data to cluster)
Q = zeros(RC, nA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:ncluster
    center(i,:)=unifrnd(lb,ub);
end
Fit=zeros(50000,1);
iter=1;
Epsilon=InitialEpsilon;
oldreward=ones(ncluster,1);
cellArray = cell(ncluster,dp);
newcenter=zeros(ncluster,dp);
while true
    IndxD=[];
    dis=zeros(np,1);
    clu2=zeros(np,1);   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%find center%%%%%%%%%%%%%%%%%%%%%%%%%%
if iter==1
    for j=1:np 
        disPcenter2=zeros(ncluster,1);
        for i=1:ncluster  
             disPcenter2(i)=norm((center(i,:)-P(j,:)));%.*W(i,:));
        end
        [dis2(i),clu2(j)]=min(disPcenter2);
    end
     for i=1:ncluster 
             ui=find(clu2==i);
             up=P(ui,:);   
             cellArray{i}=up;
     end
end
iter=iter+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Within cluster scatter 
     for l=1:ncluster  
         up=cellArray{l};
       for j = 1 : size(up,2)
           distoceterW = norm((up(:,j)- center(l,j)));
           SSW(l)=sum(distoceterW);
       end
     end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % Between cluster scatter             
        for i=1:ncluster           
            for j=1:ncluster
                if (i==j)
                    continue;
                else
                 SSB(i)=(norm(center(i,:)- center(j,:)));
                end
            end
                [SSBl(i),Indxssb]=max(SSB);
        end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find reward of each cluster                      
        for i=1:ncluster
            cellsize=size(cellArray{1},1);
           if (cellsize==0)
               cellsize=1;
           end
             density(i)=SSW(i)./cellsize; 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% density other cluster
       for l=1:ncluster 
        for j=1:ncluster
            if (l==j)
                densityanothercenter2(l,j)=-1000;
            else
             densityanothercenter2(l,j)=density(j);
            end
        end
      end
        
        %%%%%%%%%%%%%%%%%%%%%step 1-find farthest data in each cluster
    fardistocenterIndx=[]; 
    for l=1:ncluster  
       up1=[];
       fardistocenter1=[];
       up1=cellArray{l};
       if ((size(up1,1)==0)||(size(up1,1)==1))
           fardistocenter1(l)=0;IndxFar=0;
           fardistocenterIndx(l)=0;           
       else
       for k=1:size(up1,1)                                
              fardistocenter1(k)= norm((up1(k,:)- center(l,:)));%.*W(l,:));              
       end 
       
        [~,IndxFar]=max(fardistocenter1);
        fardistocenterIndx(l)=IndxFar; 
       
       end
    end
     %%%%%%%%%%%%%%%%%%%%%1 find distance of far data in each cluster to another  and dense other cluster  
      disWorstanothercenter2=zeros(ncluster,ncluster);
      densityanothercenter2=zeros(ncluster,ncluster); 
     for l=1:ncluster  
         up1=cellArray{l};
        for j=1:ncluster
            if (l==j)||(fardistocenterIndx(l)==0)
                disWorstanothercenter2(l,j)=1000;
            else
                %size(up1)
               % fardistocenterIndx(l)
             disWorstanothercenter2(l,j)=norm(up1(fardistocenterIndx(l),:)- center(j,:));  
            end
        end 
      end

  
 randlist=[];
     for l=1:ncluster  
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2-find nearest center except i to this worst data
         [disWorst,cluW]=min(disWorstanothercenter2(l,:));  
         [disWorsts,cluWs]=sort(disWorstanothercenter2(l,:),2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3-find cluster with max density
         [densityD,IndxD]=max(densityanothercenter2(l,:));    
         %IndxD = IndxD(randi(numel(IndxD)));
         [densityDs,IndxDs]=sort(densityanothercenter2(l,:),2,'descend');
     
     
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    % 4- find cluster that  satisfiy    2&3
       
         [FIndxD,IndxDF1]=find(IndxD==cluW);
     
        if isempty(IndxDF1)           
            [SIndxD,IndxDF2]=find(IndxDs(2)==cluWs(2));
            if isempty(IndxDF2)
                for kk=1:ncluster
                if kk==l
                continue;
                else
                randlist=[randlist,kk];
                end 
                end
                nextcenter=randlist(randi(size(randlist,2))); %randlist(aa);
               % disp('hi3');
            else
              nextcenter=cluWs(2);
              %disp('hi2');
            end
            
        else
              nextcenter=cluW;
              %disp('hi1');
        end
                         
% i
% nextcenter
         naction = EpsilonGreedyActionSelection(l, Q, Epsilon);
         if naction==1
             newreward=SSW(l);
         elseif naction==2
             newreward=-SSW(l);
         else 
             newreward=0;%+1.5*norm(center(i,:)- newcenter(i,:));
         end
    
        
     if naction==1 %%%%%%%%%%%%%epsilon greedy must used here if naction=1 (add data for next cluster)
        up=cellArray{l};
        newcenter(l,:)=mean(up);
        up2=cellArray{nextcenter};
       if ((size(up,1)==0)||(size(up,1)==1)|| (fardistocenterIndx(l)==0))
           fardistocenterIndx(l)=0;  
           adddata=[];
       else
       % fardistocenterIndx(l)
        %size(up)
        adddata=up(fardistocenterIndx(l),:);
       end
        up3=[up2;adddata];
      %  cellArray{IndxD}=[];
        cellArray{nextcenter}=up3;
        newcenter(nextcenter,:)=mean(up3);
     
        elseif naction==2   %%%%%%%%%%%%% if naction=2 (remove data)
        up=cellArray{l};      
        %fardistocenterIndx(l)
        %size(up)
       if ((size(up,1)==0)||(size(up,1)==1)||(fardistocenterIndx(l)==0))
           fardistocenterIndx(l)=0;  
           adddata=[];
       else
        adddata=up(fardistocenterIndx(l),:);
        up(fardistocenterIndx(l),:)=[];
       end
        
   %     cellArray{i}=[];
        cellArray{l}=up; 
        newcenter(l,:)=mean(up);
        up2=cellArray{nextcenter};
        up3=[up2;adddata];
         %  cellArray{IndxD}=[];
       cellArray{nextcenter}=up3;   
       newcenter(nextcenter,:)=mean(up3);
    else 
        up=cellArray{l};      %%%%%%%%%%%%% if naction=3 (no action)           
        newcenter(l,:)=mean(up);
     end
 % nextcenter=[];
 %max(Q(nextcenter,:))
 %size(clu2)
%%%%%%%%%%%i is number of currnet center and IndxD is number of next center
        currentcenter=l;   
        Q(currentcenter, naction) = Q(currentcenter, naction) + alpha*(newreward + gamma*max(Q(nextcenter,:)) - Q(currentcenter,naction));                  
     end

if (iter>1)
    for j=1:np 
        disPcenter2=zeros(ncluster,1);
        for r=1:ncluster 
                disPcenter2(r)=norm((newcenter(r,:)-P(j,:)));
        end
        [dis2(r),clu2(j)]=min(disPcenter2);
    end
  %center 
 for r=1:ncluster 
         ui=find(clu2==r);
         up=P(ui,:);   
         cellArray{r}=up;
 end
end
 
    alpha = alpha * alphaDamp;
    Epsilon = Epsilon * EpsilonDamp;
    if isequal(newcenter,center)
        break
    else
        center=newcenter;
    end
if iter==1000
    break
end
end    

end%(while)