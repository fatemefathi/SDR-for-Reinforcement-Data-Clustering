function d=computeprob(a,b)

M=size(a,1);
sum0=0;
aa=normalize(a,'norm',1);bb=normalize(a,'norm',1);
for k=1:M
    sum0=sum0+(abs(aa(k)-bb(k)));
end

for k=1:M
  if (sum0==0)
      P(k)=1/M;
  else
      P(k)=(1/(M-1))*(1-((abs(aa(k)-bb(k)))/sum0));
  end 
end
P=P';
if (size(a,1) == 1)
  a = [a; zeros(1,size(a,2))]; 
  b = [b; zeros(1,size(b,2))]; 
end
a=a.*P;b=b.*P;
aa=sum(a.*a); bb=sum(b.*b); ab=a'*b; 
d = repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab;
d = real(d);
d = max(d,0);

% % force 0 on the diagonal? 
% if (df==1)
%   d = d.*(1-eye(size(d)));
% end
 
