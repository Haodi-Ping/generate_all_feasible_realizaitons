function [MOCC,radio] = merge(subgraph,G,ConnectivityM)
[a,b]=size(subgraph);
[c,d]=size(ConnectivityM);
subnode=zeros(a,c);
for i=1:a
    id = find(subgraph(i,:));
    subnode(i,G.Edges.EndNodes(id,1)')=1;
    subnode(i,G.Edges.EndNodes(id,2)')=1;
end

num = a;
lastnum = 0;
redundant = zeros(1,a)+1;
radio = zeros(1,a);
for i=1:num
    r=find(subnode(i,:));
    [lll,nnn]=size(r);
    edge=sum(sum(ConnectivityM(r,r)))/2;
    radio(i)=(edge-2*nnn+3)*2/nnn/(nnn-1);
end
while(lastnum~=num)
    ans=0;
    ansi=0;
    ansj=0;
    lastnum = num;
    for i = 1:a
       for j = i+1:a
           overnode = sum(subnode(i,:).*subnode(j,:));
           if(overnode>2)
                test=subnode(i,:)+subnode(j,:);
                r=find(test);
                [lll,nnn]=size(r);
                edge=sum(sum(ConnectivityM(r,r)))/2;
                MergeRadio=(edge-2*nnn+3)*2/nnn/(nnn-1);
               if(MergeRadio>radio(i)&&MergeRadio>(radio(j))&&MergeRadio>ans)
                   ans= MergeRadio;
                   ansi=i;
                   ansj=j;
               end
           end
           
       end
    end
    if(ans>0)
        radio(ansi)=ans;
        radio(ansj)=0;
        subnode(ansi,:)=subnode(ansi,:)+subnode(ansj,:);
        idl=find(subnode(ansi,:));
        subnode(ansi,idl)=1;
        subnode(ansj,:)=0;
        num=num-1;
    end
    
    
end
id=find(radio);
MOCC=subnode(id',:);
radio = radio(id);




end

