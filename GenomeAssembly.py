'''
Input: A sequence path of k-mers Pattern_1, … ,Pattern_n such that the last k - 1 symbols of Pattern_i are equal to the first k-1 symbols of Pattern_i+1 for 1 ≤ i ≤ n-1.
Output: A string Text of length k+n-1 such that the i-th k-mer in Text is equal to Pattern_i (for 1 ≤ i ≤ n).
'''     
def PathToGenome(path: list[str]) -> str:
    k=len(path[0])

    for i in range(len(path)):
        if i == 0:
            # take the whole string in the forst k-mer as a start
            genome=path[0]
        else:
            # For the remaining k-mers, add the last letter to the sequence one by one
            genome=genome+path[i][k-1]

    return genome

'''
Input: A collection of k-mers Patterns.
Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).
'''
def DeBruijn(k_mers: list[str]) -> dict[str: list[str]]:
    k=len(k_mers[0])
    mer_dict={}
    
    # the prefix of each k-mer is a node in the DeBruijn graph (a key in the dictionary)
    for str in k_mers :
        prefix=str[:k-1]
        # if the key has not been created, create this prefix as a key.
        if prefix not in mer_dict.keys():
            mer_dict[prefix]=[]
        # append the suffix of this k-mer to the value list of the prefix
        mer_dict[prefix].append(str[1-k:])
            
    return mer_dict

'''
Input: The adjacency list of an Eulerian directed graph.
Output: An Eulerian cycle in this graph.   
'''
def EulerianCycle(G: dict[str, list[str]]) -> list[str]:
    # pick a random start to find cycles in the graph
    key_list=list(G.keys())
    next_node=key_list[0]
    path=[]

    # when there's no more node to explore, stop the cycle searching process
    while next_node != None:
        # each time find a cycle starting from next_node, update the visited path and the remaining graph
        path, G, next_node = Findcycle(path,G,next_node)

    # append the first node to the end to form an actual circle
    path.append(path[0])

    return path

'''
Based on the key, the current path and the graph, find a cycle starting from the key
append the new cycle to the visited path, return the updated cycle, the remaining graph and the next node to be explored if there is any.
'''
def Findcycle(path: list[str], G: dict[str, list[str]], key: str):
    # if the path is empty, simply append the new visited nodes to the path.
    # however if the path is not empty, swap the original path so that the path will start and end from THE key.
    # And then append the newly visited nodes after THE key.
    if len(path)!=0:
        keyindex=path.index(key)
        new_path=path[keyindex:]
        new_path.extend(path[:keyindex])
        path=new_path

    # from the initial node keep visiting its child until it find itself again
    next_node=key
    while G[next_node][0] != key:
        # Whenever a node is visited, append it to the path and remove the edge from the graph
        path.append(next_node)
        next_node=G[next_node].pop(0)
    path.append(next_node)
    G[next_node].pop(0)

    next_start=None
    # Assuming that the graph is connected, there will always be a node in the path
    # that has not been thoroughly visited until the whole graph is traversed.
    # So each time return a random node in the path that still has children as the starting point in the next iteration.
    for node in path:
        if len(G[node])>0:
            next_start=node
            break
    
    return path, G, next_start

'''
Input: The adjacency list of a directed graph that has an Eulerian path.
Output: An Eulerian path in this graph.
'''
def EulerianPath(G: dict[str, list[str]]) -> list[str]:
    indegdict=IndegreeDict(G)
    keylist=list(indegdict.keys())
   
   # find the startnode and the endnode in the graph and connect them 
   # so that the graph become cyclic.
    for key in keylist :
        # the node with indegree less than outdegree is the startnode
        if indegdict[key]<len(G[key]):
            startnode=key
        # the node with outdegreee more than indegree is the end node.
        elif indegdict[key]>len(G[key]):
            endnode=key 
    G[endnode].append(startnode)

    # Obtain the EulerianCycle in the connected graph
    cycle=EulerianCycle(G)
    # Remove the start to make it a path
    cycle.pop(0)

    # From the path, locate the edge that was added 
    for i in range(len(cycle)):
        if cycle[i]==endnode and cycle[i+1]==startnode:
            index=i
            break
    # swap the path and break the additional edge
    new_path=cycle[index+1:]
    new_path.extend(cycle[:index+1])

    return new_path
'''
Input: An integer k followed by a list of k-mers Patterns.
Output: A string Text with k-mer composition equal to Patterns. 
'''
def StringReconstruction(k: int, Patterns: list[str]) -> str:
    # From a collection of k-mers, first construct the according DeBruijn graph
    DBgraph=DeBruijn(Patterns)
    # Then find a Eulerianpath in the graph
    path=EulerianPath(DBgraph)
    # Finally reconstruct the sequence from the path
    text=PathToGenome(path)

    return text

'''
A k-universal string contains each of the binary k-mers exactly once.
Input: An integer k.
Output: A k-universal circular string.
'''
def KUniversalString(k: int) -> str:
    # First obtain all Binary k-mers, then find a Eulerian cycle in the according Debruijn graph
    k_mers=BinaryStrings(k)
    DBgraph=DeBruijn(k_mers)
    path=EulerianCycle(DBgraph)
    # remove the starting node to avoid repetition
    path.pop(0)
    text=PathToGenome(path)[k-1:]
    return text

'''
Create a list of all the binary k-mer strings.
'''
def BinaryStrings(k: int) -> list[str]:
    KmerList=[]
    for i in range (2**k):
        binary=str(bin(i)).lstrip('0b')
        # fill up the k digits with '0'
        binary='0'*(k-len(binary))+binary
        KmerList.append(binary)

    return KmerList

'''
Input: A collection of k-mers Patterns.
Output: All contigs in DeBruijn(Patterns).
A path in a graph is called non-branching if in(v) = out(v) = 1 for each intermediate node v of this path.
A contig is a maximal non-branching path.
'''
def ContigGeneration(Patterns: list[str]) -> list[str]:
    G=DeBruijn(Patterns)
    contigs=[]

    indict=IndegreeDict(G)
    outdict=OutdegreeDict(G)
    # internodes are the intermediate nodes where in(v) = out(v) = 1
    internodes=Interim(indict,outdict)

    # Start from each internode, find a non-branching path by spreading to both ends as long as possible.
    for currentnode in internodes:
        # In case this internode has been included in another path
        if len(G[currentnode])<1:
            continue
        
        contig=[currentnode]
        parent=GetParent(G,currentnode)
        node=currentnode
        # keep walking left until the node is not an intermediate node (it has no parent or multiple parents)
        while len(parent)==1:
            left=parent[0]
            contig.insert(0,left)
            # Remove the edge once visited
            G[left].remove(node)
            node=left
            parent=GetParent(G,node)
        

        child=G[currentnode]
        node=currentnode
        # keep walking right until the node is not an intermediate node (it has no child or multiple children)
        while len(child)>0:
            contig.append(child[0])  
            if child[0] in internodes:
                node=G[node].pop(0)
                child=G[node]
            else:
                G[node].remove(child[0])
                break
        contigs.append(PathToGenome(contig))

    # append the contigs without internodes
    for key in G.keys():
        for value in G[key]:
            if len(value)>0:
                contigs.append(PathToGenome([key,value]))

    return contigs

'''
Obtain a list of nodes with indegree and outdegree both ==1
'''
def Interim(indict,outdict):
    interim=[]
    for key in outdict.keys():
        if indict[key]==1 and outdict[key]==1:
            interim.append(key)
    return interim

#Creates a dictionary of the keys' outdegrees
def OutdegreeDict(G):
    outdict={}
    for key in G.keys():
        outdict[key]=len(G[key])
    return outdict

#Get a list of the parents of the node the graph
def GetParent(G,childnode) -> str:
    parent=[]
    for key in G.keys():
        for value in G[key]:
            if value==childnode:
                parent.append(key)
    return parent

#Creates a dictionary of the keys' indegrees
def IndegreeDict(G):
    keylist=list(G.keys())

    indegdict={}
    for key in keylist:
        indegdict[key]=0

    for key in keylist:
        for j in G[key]:
            if j not in keylist:
                indegdict[j]=0
                G[j]=[]
            indegdict[j]+=1
    return indegdict

'''
Same as ContigGeneration but need to change datatypes of the sub functions
def MaximalNonBranchingPaths(G: dict[int, list[int]]) :
    contigs=[]
    indict=IndegreeDict(G)
    outdict=OutdegreeDict(G)
    internodes=Interim(indict,outdict)

    for currentnode in internodes:
        if len(G[currentnode])<1:
            continue
        contig=[currentnode]
        parent=GetParent(G,currentnode)
        node=currentnode
        while len(parent)==1:
            left=parent[0]
            contig.insert(0,left)
            G[left].remove(node)
            if left in internodes:
                node=left
                parent=GetParent(G,node)
            else:
                break
        
        child=G[currentnode]
        node=currentnode
        while len(child)>0:
            contig.append(child[0])  
            if child[0] in internodes:
                node=G[node].pop(0)
                child=G[node]
            else:
                G[node].remove(child[0])
                break
        contigs.append(contig)

    for key in G.keys():
        if len(G[key])>0:
            for value in G[key]:
                contigs.append([key,value])

    return contigs
'''