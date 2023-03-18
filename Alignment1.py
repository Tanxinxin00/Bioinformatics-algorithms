# change the maximum recursion limit otherwise the OutputLCS funtion will be suspended easily.
import sys
sys.setrecursionlimit(10000)

'''
Input: An integer money and an array Coins = (coin1, ..., coind).
Output: The minimum number of coins with denominations Coins that changes money.
'''
def DPChange(money: int, Coins: list[int]) -> int:
    # Initializing the array of MinNumCoins[i] with a big value to be potentially overwritten
    MinNumCoins=[money]*(money+1)
    MinNumCoins[0]=0

    # Fill out the array sequentially
    for m in range(1,money+1):
        # range through all coin denominations and get the smallest MinNumCoins[money-denomination]
        for coin in Coins:
            if m >= coin:
                if MinNumCoins[m-coin] + 1<MinNumCoins[m]:
                    MinNumCoins[m]=MinNumCoins[m-coin] +1
    
    return MinNumCoins[money]
'''
Find the length of a longest path in the Manhattan Tourist Problem.
Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right. The two matrices are separated by the "-" symbol.
Output: The length of a longest path from source (0, 0) to sink (n, m) in the rectangular grid whose edges are defined by the matrices Down and Right.
'''
def LongestPathLength(n: int, m: int, Down: list[list[int]], Right:list[list[int]]) -> int:
    # construct a (n+1)x(m+1) matrix
    grid=[]
    for i in range(n+1):
        grid.append([0]*(m+1))
    grid[0][0]=0

    # First fill out the first row and the first column
    for c in range(1,m+1):
        grid[0][c]=grid[0][c-1]+Right[0][c-1]
    for r in range(1,n+1):
        grid[r][0]=grid[r-1][0]+Down[r-1][0]
    
    # Then fill out the other squares row by row
    for r in range(1, n+1):
        for c in range(1,m+1):
            grid[r][c]=grid[r][c-1]+Right[r][c-1]
            # For each square, fill in the max score between two directions according to the DOWN&RIGHT matrix
            if grid[r-1][c]+Down[r-1][c]>grid[r][c]:
                grid[r][c]=grid[r-1][c]+Down[r-1][c]
    
    return grid[n][m]



'''
Input: Two strings s and t.
Output: A longest common subsequence of s and t. 
'''
def LongestCommonSubsequence(s: str, t: str) -> str:
    numcol=len(s)+1
    numrow=len(t)+1

    track=LCSBackTrack(s,t)
    return OutputLCS(track,s,numrow-1,numcol-1)

'''
Given two strings, find the maximum matches while keep track of route in the Manhattan graph
'''
def LCSBackTrack(v: str, w: str) -> list[list[int]]:
    numcol=len(v)+1
    numrow=len(w)+1

    # construct a numcolxnumrow value matrix
    s=[]
    for i in range(numrow):
        row=[0]*(numcol)
        s.append(row)

    # construct a numcolxnumrow track direction matrix
    Backtrack=[]
    for i in range(numrow):
        row=[0]*(numcol)
        Backtrack.append(row)  

    # Fill out the value matrix sequentially, while keep track of the backtrack route
    for i in range(1,numrow):
        for j in range(1,numcol):
            # if the nucleutide is the same, count as match 
            match=0
            if w[i-1]==v[j-1]:
                match=1
            s[i][j]=max(s[i-1][j],s[i][j-1],s[i-1][j-1]+match)

            if s[i][j]==s[i-1][j]:
                Backtrack[i][j]="↓"
            elif s[i][j]==s[i][j-1]:
                Backtrack[i][j]="→"
            else:
                Backtrack[i][j]="↘"

    return Backtrack


'''
Given a recorded track route in a Manhattan graph, return the matches along the route
'''
def OutputLCS(backtrack: list[list[str]], v: str, i: int, j: int) -> str:
    # is i=0 or j=0, it means at least one string has been traversed, no more common nuleutide
    if i == 0 or j == 0:
        return ""
    elif backtrack[i][j] == "↓":
        return OutputLCS(backtrack, v, i - 1, j)
    elif backtrack[i][j] == "→":
        return OutputLCS(backtrack, v, i, j - 1)
    else:
        # only when the track goes in digonal, there is a match
        return OutputLCS(backtrack, v, i - 1, j - 1) + v[j-1]

'''
Input: An integer representing the starting node to consider in a graph, followed by an integer representing the ending node to consider, followed by a list of edges in the graph. The edge notation "0 1 7" indicates that an edge connects node 0 to node 1 with weight 7.  You may assume a given topological order corresponding to nodes in increasing order.
Output: The length of a longest path in the graph, followed by a longest path as a sequence of space-separated node labels. (If multiple longest paths exist, you may return any one.)
'''
def LongestPath(s: int, t: int, E: dict[int, list[tuple[int, int]]]) -> tuple[int, list[int]]:
    # sort the keys in increasing order(also the topological order of the DAG)
    keys=list(E.keys())
    keys.sort()
    # construct a maximum value dictionary
    # each node in the dictionary corresponds to (value, path)
    # where value is the maximum cumulative edge weight from the source node to the node
    # path is the corresponding path from the source node to this node
    maxdict={}
    maxdict[s]=(0,[])

    for node in keys:
        # every node reached starting from 0 will be reached through edges first
        # so if the node is not in keys, means that the node is a source node
        # we need to start from 0 node, so we can simply omit the paths starting from other source nodes.
        if node not in maxdict.keys():
            continue
        # for each edge, compute the cumulative edge value for the childnode, and also keep track of the visited path
        for (nextnode, weight) in E[node]:
            if nextnode not in maxdict.keys():
                path=Copy(maxdict[node][1])
                path.append(node)
                maxdict[nextnode]=(maxdict[node][0]+weight,path)
            # update the maximum value for each node
            elif maxdict[node][0]+weight>maxdict[nextnode][0]:
                path=Copy(maxdict[node][1])
                path.append(node)
                maxdict[nextnode]=(maxdict[node][0]+weight,path)
    
    # add the last node to the path
    maxdict[t][1].append(t)
    return maxdict[t]

'''
Deep copy a list
'''
def Copy(listB:list) :
    listA=[]
    for i in range(len(listB))  :
        listA.append(listB[i])
    return listA