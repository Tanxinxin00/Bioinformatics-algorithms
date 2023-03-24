'''
This is Week 8 code submission for 02604 Bioinformatics from Tanxin(Thea) Qiao.
Written in Python=3.9.13 
03/24/2023
'''

import math
'''
Distances Between Leaves Problem: Compute the distances between leaves in a weighted tree.

Input:  An integer n followed by the adjacency list of a weighted tree with n leaves.
Output: An n x n matrix (di,j), where di,j is the length of the path between leaves i and j.

'''
def DistofLeaves(numLeaves: int, adjlist: dict[int,dict[int,int]])->dict[int,dict[int,int]]:
    distances={}
    for i in range(numLeaves):
        distances[i]={}
    for leave1 in range(numLeaves):
        for leave2 in range(leave1,numLeaves):
            if leave1==leave2:
                dist=0
            else:
                dist=FindDist(adjlist,set(),leave1,leave2)
                distances[leave1][leave2]=dist


            distances[leave1][leave2]=dist
            distances[leave2][leave1]=dist
    
    return distances       

'''
FindDist performs a depth-first-search-like search to find the distance between two nodes in an acyclic graph.
It takes an adjacency list, a set of visited nodes, the index of the starting node and ending node as input.
And outputs the distance between the starting node and the ending node.
'''
def FindDist(adjlist: dict[int,dict[int,int]],visited,start:int, end:int)-> int:
    visited.add(start)

    for nbr in adjlist[start].keys(): 
        if nbr not in visited:
            edgeweight=adjlist[start][nbr]

            if nbr==end:
                visited.add(end)
                return edgeweight
            else:
                dist=edgeweight+FindDist(adjlist,visited,nbr,end)
                if dist>0:
                    return dist
   
    # if the current path does not reach the end node in any way, return negative infinity and it will be bypassed in the former recursion
    return -math.inf


'''
Solve the Limb Length Problem.

Input: An integer n, followed by an integer j between 0 and n - 1, followed by a space-separated additive distance matrix D (whose elements are integers).
Output: The limb length of the leaf in Tree(D) corresponding to row j of this distance matrix (use 0-based indexing).
'''
def LimbLength(num:int, node:int,distmtx:dict[int,dict[int,int]])-> int:
    #randomly select a node(except from the target node)
    if node!=0:
        i=0
    else:
        i=1
    
    # find the min of D(node,i)+D(node,k)-D(i,k) where k is another node
    minlimb=math.inf
    for k in range(num):
        if k!=i and k!= node:
            limb=distmtx[node][i]+distmtx[node][k]-distmtx[i][k]
            if limb<minlimb:
                minlimb=limb
    
    return minlimb/2

'''
Implement UPGMA.

Input: An integer n followed by a space separated n x n distance matrix.
Output: An adjacency list for the ultrametric tree returned by UPGMA. Edge weights should be accurate to two decimal places (answers in the sample dataset below are provided to three decimal places).
'''

def UPGMA(D, n):
    # Step 1: Initialize clusters and construct a graph T
    Clusters = {i:[i] for i in range(n)}
    # Dnew is a copy of the adjacency list provided that will be updated
    Dnew=CopyDictofDict(D)
    # each node in the tree corresponds to a cluster in clusters
    T = {i: {} for i in range(n)}

    for i in range(n):
        T[i] = {'age': 0, 'edges': {}}

    # Step 2: Merge clusters and update graph T
    newnode=n
    while len(Clusters) > 1:
        # Find the two closest clusters
        min_dist = math.inf
        ci, cj = 0, 0
        for key1 in Dnew.keys():
            for key2 in Dnew[key1]:
                if key1!=key2:
                    dist=Dnew[key1][key2]  
                    if dist<min_dist:   
                        min_dist=dist
                        ci,cj=key1,key2  

        # Merge the two closest clusters
        Cnew = Clusters[ci]+Clusters[cj]
        Clusters[newnode]=Cnew
        del Clusters[ci]
        del Clusters[cj]

        # add a new node and calculate the edge weights linking it to the two children nodes.
        T[newnode] = {'age':min_dist*0.5}

        edgetoi=T[newnode]['age']-T[ci]['age']
        edgetoj=T[newnode]['age']-T[cj]['age']

        T[newnode]['edges']={}
        T[newnode]['edges'][ci]=edgetoi
        T[newnode]['edges'][cj]=edgetoj
        T[ci]['edges'][newnode]=edgetoi
        T[cj]['edges'][newnode]=edgetoj


        # Update the distance matrix
        del Dnew[ci]
        del Dnew[cj]
        for key in Dnew.keys():
            del Dnew[key][ci]
            del Dnew[key][cj]
        Dnew[newnode]={}

        for i in Dnew.keys():
            if i !=newnode:
                dist = sum([D[leave][x] for x in Cnew for leave in Clusters[i]]) / (len(Cnew)*len(Clusters[i]))
                Dnew[newnode][i]=dist
                Dnew[i][newnode]=dist
            else:
                Dnew[i][newnode]=0

        newnode+=1

    return T

#deep copies a dictionary 
def CopyDict(P):
    copy={}
    for key in P.keys():
        copy[key]=P[key]

    return copy

#deep copies a dictionary of dictionary 
def CopyDictofDict(P):
    copy={}
    for key in P.keys():
        copy[key]=CopyDict(P[key])
    return copy