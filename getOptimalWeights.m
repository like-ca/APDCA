function [p, Ws]= getOptimalWeights(Hs, alpha)
K= length(Hs);
[sortedVal, idx]= sort(Hs,'ascend');
p=K; found=1; gamma=0;
while p>0 && found
    gamma= (sum(sortedVal(1:p))+2*alpha)/p;
    if (gamma - sortedVal(p))>0
        found=0;
    else
        p= p-1;
    end
end
newWs= zeros(K,1);
for ii=1:p
    newWs(ii)= (gamma - sortedVal(ii)) / (2*alpha);
end
Ws= zeros(K,1);
Ws(idx)= newWs;
end