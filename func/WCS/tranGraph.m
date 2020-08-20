function nodeSub = tranGraph(G, ConnectivityM,n,subgraph)
[a,b]=size(subgraph);
nodesub = zeros(a,n);
for i=1:a
   ID = find(subgraph(i,:)); 
   for j=ID
       nodesub(i,G.Edges.EndNodes(j,1))=1;
       nodesub(i,G.Edges.EndNodes(j,2))=1;       
   end
end
nodeSub=nodesub;
end

