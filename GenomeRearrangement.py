
'''
Implement GreedySorting.

Input: A permutation P.
Output: The sequence of permutations corresponding to applying GreedySorting to P, ending with the identity permutation.
'''
def GreedySorting(P: list[int]) -> list[list[int]]:
    n=len(P)
    Sorting=[]

    for i in range(n):
        if P[i]!=i+1 and P[i]!=-i-1:
            if i+1 in P:
                index=P.index(i+1)
                P=Swap(P,index,i)
                Piece=Copy(P)
                Sorting.append(Piece)
            else:
                index=P.index(-i-1)
                P=Swap(P,index,i)
                Piece=Copy(P)
                Sorting.append(Piece)
        if P[i]==-i-1:
            P[i]=i+1
            Piece=Copy(P)
            Sorting.append(Piece)
    
    return Sorting

'''
A swap operation takes a list of interger (a permutation ) as input and reverses the signs and positions of the intergers from sindex to cindex
and returns the resulting list.
'''
def Swap(P: list[int],sindex,cindex) -> list[list[int]]:
    Pswap=P[:cindex]
    Pmiddle=P[cindex:sindex+1]
    Pend=P[sindex+1:]

    n=len(Pmiddle)
    Pmiddleswapped=[0]*n
    for i in range(n):
        Pmiddleswapped[i]=Pmiddle[n-1-i]*(-1)
    
    Pswap.extend(Pmiddleswapped)
    Pswap.extend(Pend)

    return Pswap

'''
Copies a list of interger.
'''
def Copy(L:list[int]) -> list[int]:
    n=len(L)
    newlist=[0]*n

    for i in range(n):
        newlist[i]=L[i]
    
    return newlist


'''
Number of Breakpoints Problem: Find the number of breakpoints in a permutation.

Input: A permutation.
Output: The number of breakpoints in this permutation.
'''
def BreakpointCount(P: list[int]) -> int:
    n=len(P)
    P.insert(0,0)
    P.append(n+1)

    breakpoint=0
    for i in range(n):
        if P[i+1]-P[i] !=1:
            breakpoint+=1
    
    return breakpoint


'''
Solve the 2-Break Distance Problem.

Input: Genomes P and Q.
Output: The 2-break distance d(P, Q).
'''
def TwoBreakDistance(P: list[list[int]], Q: list[list[int]]) -> int:

    Pedges=ChromosomeToEdges(P)

    numBlock=int(len(Pedges)/2)

    Qedges=ChromosomeToEdges(Q)

    cycles=len(FindCycle(Pedges,Qedges))

    return numBlock-cycles

'''
FindCycle finds every cycle in the breakpoint graph consisted of Pedges and Qedges
and returns a list of an arbituary staring node for each cycle
'''
def FindCycle(Pedges:dict[int,int],Qedges:dict[int,int])->list[int]:
    cyclestart=[]

    for key in Pedges.keys():

        if Pedges[key]==0:
            continue
        nextnode=0
        currentnode=key

        while nextnode!=key:
            node1=Pedges[currentnode]
            #Once en edge is visited, set the value to 0 to mark it as visited.
            Pedges[currentnode]=0
            Pedges[node1]=0
            nextnode=Qedges[node1]
            Qedges[node1]=0
            Qedges[nextnode]=0
            currentnode=nextnode

        else:
            cyclestart.append(key)

    return cyclestart

'''
ChromosomeToEdges takes a chromosome(list of list of integers(blocks)) as input and returns a dictionary of int-int representing edges in this chromosome.
'''
def ChromosomeToEdges(P:list[list[int]])->dict[int,int]:
    edges={}
    for graph in P:
        for i in range(len(graph)):
            if i !=len(graph)-1:
                block=graph[i]
                nextblock=graph[i+1]            
            else:
                block=graph[i]
                nextblock=graph[0]
            
            if block > 0:
                node1=2*block
            else:
                node1=2*(-block)-1


            if nextblock > 0:
                node2=2*nextblock-1
            else:
                node2=-2*nextblock

            
            edges[node1]=node2
            edges[node2]=node1
    return edges


'''
2-Break Sorting Problem: Find a shortest transformation of one genome into another by 2-breaks.

Input: Two genomes with circular chromosomes on the same set of synteny blocks.
Output: The sequence of genomes resulting from applying a shortest sequence of 2-breaks transforming one genome into the other.
'''
def TwoBreakSorting(P: list[list[int]], Q: list[list[int]]) -> list[list[list[int]]]:
    Sorting=[P]
    Pedges=ChromosomeToEdges(P)
    Qedges=ChromosomeToEdges(Q)

    numBlock=int(len(Pedges)/2)
    cycles=FindCycle(CopyDict(Pedges),CopyDict(Qedges))

    while len(cycles)<numBlock:
        node1=0
        for i in range(len(cycles)):
            if Qedges[Pedges[cycles[i]]]!=cycles[i]:
                node1=cycles[i]
        # breaks edges (node1,node2) and (node3, node4)
        node2=Pedges[node1]
        node3=Qedges[node2]
        node4=Pedges[node3]
        # add edges(node2,node3) and (node1, node4)
        Pedges[node2]=node3
        Pedges[node3]=node2
        Pedges[node4]=node1
        Pedges[node1]=node4

        Sorting.append(EdgesToChromosome(CopyDict(Pedges)))
        cycles=FindCycle(CopyDict(Pedges),CopyDict(Qedges))
    return Sorting

'''
Copies a dictionary
'''
def CopyDict(P: dict[int,int])-> dict[int,int]:
    copy={}
    for key in P.keys():
        copy[key]=P[key]

    return copy

'''
EdgesToCHromosome takes a dictionary of int-int(representing edges) as input
and returns a chromosome (consists of graphs of blocks)
'''
def EdgesToChromosome(edges: dict[int])->list[list[int]]:
    cycles=[]

    for key in edges.keys():
        cycle=[]

        if edges[key]==0:
            continue
        if key%2==0:
            block=key//2
        else:
            block=-key//2

        nextblock=0
        nextnode=key
        while nextblock!=block:
            currentnode=edges[nextnode]
            edges[nextnode]=0

            if currentnode%2==1:
                nextblock=currentnode//2 +1
                nextnode=nextblock*2
            else:
                nextblock=-currentnode//2
                nextnode=(-nextblock)*2-1

            edges[currentnode]=0
            cycle.append(nextblock)
        cycles.append(cycle)

    return cycles

'''
ChromosomeToNodeCycle takes a chromosome as input and 
outputs a list of nodes representing the heads and tails of the blocks
'''
def ChromosomeToNodeCycle(P:list[list[int]])->list[list[int]]:
    NodeCycle=[]

    for part in P:
        partnode=[]
        for block in part:
            if block>0:
                partnode.append(2*block-1,2*block)
            else:
                partnode.append(2*(-block),2*(-block)-1)
        NodeCycle.append(partnode)
    
    return NodeCycle

'''
ChromosomeToNodeCycle takes  a list of nodes representing the heads and tails of the blocks as input and 
outputs a chromosome consisting of the blocks
'''
def NodeCycleToChromosome(cycle:list[list[int]])->list[list[int]]:
    chromosome=[]

    for part in cycle:
        chromo=[]
        for i in range(len(part)):
            if i%2==0:
                if part[i]<part[i+1]:
                    block=part[i+1]/2
                else:
                    block=-part[i]/2
                chromo.append(block)
        chromosome.append(chromosome)
    
    return chromosome
    
            
