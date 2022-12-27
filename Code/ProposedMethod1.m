%function [W, y, S, evs] = PropsedMethod1(X, c, d, k, r)
 clc;
 close all;
 clear all;
%Input DataSet
%Determine X and label and number of class

% X: dim*num data matrix, each column is a data point
% c: number of clusters
% k: number of neighbors to determine the initial graph, and the parameter r if r<=0
% r: paremeter, which could be set to a large enough value. If r<0, then it is determined by algorithm with k
% A: num*num learned symmetric similarity matrix
% evs: eigenvalues of learned graph Laplacian in the iterations


if nargin < 4
    r = -1;
end;
if nargin < 3
    k = 15;
end;


d = c-1;
r = -1;
[R N]=size(X);
%%%%%%%%%compute initial Sij and Ls and P 
for i=1:N
    for j=1:N
%         size(X(:,1))
distX(i,j) = computeprob(X(:,i),X(:,j));
    end
end
[distX1, idx] = sort(distX,2);
S = zeros(N);
rr = zeros(N,1);
for i = 1:N
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end;

 

%%%%%%%%%%%%%compute F
if r <= 0
    r = mean(rr);
end;
lambda = r;

S0 = (S+S')/2;
D0 = diag(sum(S0));
L0 = D0 - S0;
[F, temp, evs]=eig1(L0, c, 0);
%%%%%%%%%%%%compute  W
M = (X*L0*X');
W = eig1(M, d, 0, 0);

H = eye(N)-1/N*ones(N);
St = X*H*X';
invSt = inv(St);
M = (X*L0*X');
W = eig1(M, d, 0, 0);

%landa=.2;
beta=0.9;
aa=normalize(X,'norm',1);
landaij=rand(N,N)/N;
Iter=30;
for iter = 1:Iter
       iter
    for i=1:N
      for j=1:N
          Probnew=[];s=S(i,j);landa=landaij(i,j);
        for kk=1:R
           Probnew(kk)= exp((landa-beta-(aa(kk,i)-aa(kk,j)).^2*s)/beta);
        end      
            Probnew(:)=Probnew(:)/sum(Pnew);
            P=Probnew';
            final=sum(((X(:,i)-X(:,j))).^2.*P);
            difdist(i,j)=final;
       end
    end
    %size(W'*difdist)
    distx2 =difdist; %L2_distance_1(W'*difdist,W'*difdist);
    S = zeros(N);
    for i=1:N
        idxa0=1:N;
        dxi = distx2(i,idxa0);
        ad = -(dxi)/(2*r);
        S(i,idxa0) = EProjSimplex_new(ad);
     end;
         
    S = (S+S')/2;
    D = diag(sum(S));
    L = D-S;    
    F_old = F;
    [F, temp, ev]=eig1(L, c, 0);
    evs(:,iter+1) = ev;

    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > 0.000000001
        lambda = 2*lambda;
    elseif fn2 < 0.00000000001
        lambda = lambda/2;  F = F_old;
    else
        break;
    end;

end;



% call RLClustering function
