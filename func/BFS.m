function [ subResult ] = BFS( Cmartix )
%search connected graph by BFS
%disp('BFS Begin')
%SizeInput = size(Cmartix)
subResult = [];
[n,~]=size(Cmartix);
graphnum=1;
tonum =0;
total = zeros(1,n);
while (tonum<n)
    Mark=zeros(1,n);
    for i=1:n
        if total(i)<1
            rootID=i;
            break;
        end
    end
    lastMark = Mark;
    Mark = Cmartix(rootID,:);
    nodenum = sum(Mark);
    while nodenum > 0
        Nnode=find((Mark-lastMark)>0);
        lastMark =Mark;
        for i = Nnode
            Mark=Mark+Cmartix(i,:);
        end
        Mark(find(Mark>0))=1;
        nodenum = sum(Mark-lastMark);
    end
    subResult(graphnum,:)=Mark;
    tonum = sum(sum(subResult));
    total = total + Mark;
    graphnum=graphnum+1;
end

end

