import math
import numpy as np
'''
Solve the Global Alignment Problem.

Input: A match reward, a mismatch penalty, an indel penalty, and two nucleotide strings.
Output: The maximum alignment score of these strings followed by an alignment achieving this maximum score.
'''
def GlobalAlignment(match_reward: int, mismatch_penalty: int, indel_penalty: int, s: str, t: str) -> tuple[int, str, str]:
    numcol=len(s)+1
    numrow=len(t)+1

    # First create the 'score' and 'Backtrack' matrix
    score=[]
    for r in range(numrow):
        row=[-math.inf]*(numcol)
        score.append(row)

    Backtrack=[]
    for r in range(numrow):
        row=[0]*(numcol)
        Backtrack.append(row)
    score[0][0]=0

    # Initialize the score matrix
    for r in range(1,numrow):
        score[r][0]=score[r-1][0]-indel_penalty
    for c in range(1,numcol):
        score[0][c]=score[0][c-1]-indel_penalty

    # Fill in the matrix one by one according to their recursion relationship
    for r in range(1, numrow):
        for c in range(1,numcol):
            if s[c-1]==t[r-1]:
                diagonal_score=match_reward
            else:
                diagonal_score=-mismatch_penalty
            score[r][c]=max(score[r-1][c-1]+diagonal_score,score[r-1][c]-indel_penalty,score[r][c-1]-indel_penalty)
            # keep track of the path while filling matirx
            if score[r][c]==score[r-1][c]-indel_penalty:
                Backtrack[r][c]="↓"
            elif score[r][c]==score[r][c-1]-indel_penalty:
                Backtrack[r][c]="→"
            else:
                Backtrack[r][c]="↘"

    s_align,t_align=AlignBacktrack(Backtrack,s,t,numrow-1,numcol-1)

    return (score[numrow-1][numcol-1],s_align,t_align)

# ALignBacktrack takes the recorded path of an alignment, its location in the graph, the two strings in the alignment as input
# and outputs the aligned strings
def AlignBacktrack(track: list[list[str]], s: str, t: str, r: int, c: int) -> tuple[str,str]:
    if r==0:
        return(s[0:c],'-'*c)
    if c==0:
        return('-'*r,t[0:r])
    if track[r][c]=="↓":
        s_align,t_align=AlignBacktrack(track,s,t,r-1,c)
        return(s_align+'-',t_align+t[r-1])
    elif track[r][c]=="→":
        s_align,t_align=AlignBacktrack(track,s,t,r,c-1)
        return(s_align+s[c-1],t_align+'-')
    else:
        s_align,t_align=AlignBacktrack(track,s,t,r-1,c-1)
        return(s_align+s[c-1],t_align+t[r-1])

# Converts a list of strings to a list of integers
def ConvertStrListToInt(StrList: list[str]) -> list[int]:
    newlist=[]
    for str in StrList:
        newlist.append(int(str))
    return newlist 

# Read the PAM250 matrix
with open('PAM250.txt','r') as f:
    PAM={}
    lines=f.readlines()
    aminos=lines[0].lstrip(' ').rstrip('\n').split('  ')
    for i in range(1,len(lines)):
        lines[i]=lines[i].rstrip('\n').split(' ')
        while '' in lines[i]:
            lines[i].remove('')
        PAM[lines[i][0]]={}
        scores=ConvertStrListToInt(lines[i][1:])
        for c in range(len(aminos)):
            PAM[lines[i][0]][aminos[c]]=scores[c]
    

'''
Solve the Local Alignment Problem.
Input: Two protein strings written in the single-letter amino acid alphabet.
Output: The maximum score of a local alignment of the strings, followed by a local alignment of these strings achieving the maximum score. Use the PAM250 scoring matrix for matches and mismatches as well as the indel penalty σ = 5.
'''
def LocalAlignmentPAM(PAM: dict[str,dict[str,int]], indel_penalty: int, s: str, t: str) -> tuple[int, str, str]:
    numcol=len(s)+1
    numrow=len(t)+1

    # First create the 'score' and 'Backtrack' matrix
    score=[]
    for r in range(numrow):
        row=[0]*(numcol)
        score.append(row)

    Backtrack=[]
    for r in range(numrow):
        row=[0]*(numcol)
        Backtrack.append(row)

    # Initialization
    score[0][0]=0
    maxscore=0

    # Fill out the matrix one by one
    for r in range(1, numrow):
        for c in range(1,numcol):
            diagonal_score=PAM[s[c-1]][t[r-1]]
            # if the maximum socre goes below zero, reset it to zero
            score[r][c]=max(0,score[r-1][c-1]+diagonal_score,score[r-1][c]-indel_penalty,score[r][c-1]-indel_penalty)

            # keep track of the maximum score and its location
            if score[r][c]>maxscore:
                maxscore=score[r][c]
                maxcoordinate=(r,c)
            
            # keep track of the highest-scoring path
            if score[r][c]==score[r-1][c]-indel_penalty:
                Backtrack[r][c]="↓"
            elif score[r][c]==score[r][c-1]-indel_penalty:
                Backtrack[r][c]="→"
            elif score[r][c]==score[r-1][c-1]+diagonal_score:
                Backtrack[r][c]="↘"
            else:
                Backtrack[r][c]="X"

    #backtrack the highest-scoring substring alignment
    s_align,t_align=AlignBacktrackLocal(Backtrack,s,t,maxcoordinate[0],maxcoordinate[1])

    return (maxscore,s_align,t_align)

# ALignBacktrackLocal takes the recorded path of an alignment, its location in the graph, the two strings in the alignment as input
# and outputs the aligned substrings
def AlignBacktrackLocal(track: list[list[str]], s: str, t: str, r: int, c: int) -> tuple[str,str]:
    if r==0 or c==0:
        return('','')
    if track[r][c]=="↓":
        s_align,t_align=AlignBacktrackLocal(track,s,t,r-1,c)
        return(s_align+'-',t_align+t[r-1])
    elif track[r][c]=="→":
        s_align,t_align=AlignBacktrackLocal(track,s,t,r,c-1)
        return(s_align+s[c-1],t_align+'-')
    elif track[r][c]=="X":
        return ('','')
    else:
        s_align,t_align=AlignBacktrackLocal(track,s,t,r-1,c-1)
        return(s_align+s[c-1],t_align+t[r-1])

'''
Solve the Edit Distance Problem.

Input: Two strings.
Output: The edit distance (the lowest amount of insertion, deletion and mismatch) between these strings.
'''
def EditDistance(s: str, t: str) -> int:
    match_reward=0
    mismatch_penalty=1
    indel_penalty=1
    score,s2,t2 =GlobalAlignment(match_reward,mismatch_penalty, indel_penalty,s,t)
    return (-1)*score

'''
Solve the Fitting Alignment Problem.

Input: Two amino acid strings.
Output: A highest-scoring fitting alignment between v and w,
where “Fitting” w to v requires finding a substring v′ of v that maximizes the global alignment score between v′ and w among all substrings of v. 
Use the BLOSUM62 scoring table.
'''
def FittingAlignment(s: str, t: str, BLOSUM: dict[str, dict[str, int]], indel_penalty) -> tuple[int, str, str]:
    numcol=len(s)+1
    numrow=len(t)+1

    #Establish the matrixes
    score=[]
    for r in range(numrow):
        row=[0]*(numcol)
        score.append(row)

    Backtrack=[]
    for r in range(numrow):
        row=[0]*(numcol)
        Backtrack.append(row)
    score[0][0]=0

    #Initialization
    for r in range(1,numrow):
        score[r][0]=score[r-1][0]-indel_penalty

    #Fill in the score matrix according to recursion, keep track of path
    for r in range(1, numrow):
        for c in range(1,numcol):
            diagonal_score=BLOSUM[s[c-1]][t[r-1]]
            score[r][c]=max(score[r-1][c-1]+diagonal_score,score[r-1][c]-indel_penalty,score[r][c-1]-indel_penalty)

            if score[r][c]==score[r-1][c]-indel_penalty:
                Backtrack[r][c]="↓"
            elif score[r][c]==score[r][c-1]-indel_penalty:
                Backtrack[r][c]="→"
            elif score[r][c]==score[r-1][c-1]+diagonal_score:
                Backtrack[r][c]="↘"
            else:
                Backtrack[r][c]="X"
    
    # Find the highet peak in the matrix and trace back
    maxcol=np.argmax(score[numrow-1])
    maxscore=score[numrow-1][maxcol]
    s_align,t_align=AlignBacktrackLocal(Backtrack,s,t,numrow-1,maxcol)

    return (maxscore,s_align,t_align)

'''
Construct a highest-scoring overlap alignment between two strings.
Input: Two strings and a matrix score.
Output: A highest-scoring overlap alignment between the two strings as defined by the scoring matrix score.
An overlap alignment of strings v = v1 ... vn and w = w1 ... wm is a global alignment of a suffix of v with a prefix of w. 
An optimal overlap alignment of strings v and w maximizes the global alignment score 
between an i-suffix of v and a j-prefix of w (i.e., between vi ... vn and w1 ... wj) among all i and j.
'''
def OverlapAlignment( match_reward: int, mismatch_penalty: int, indel_penalty: int,s: str, t: str,) -> tuple[int, str, str]:
    # swapping the row and column
    a=s
    s=t
    t=a
    numcol=len(s)+1
    numrow=len(t)+1

    #Establish the matrixes
    score=[]
    for r in range(numrow):
        row=[0]*(numcol)
        score.append(row)

    Backtrack=[]
    for r in range(numrow):
        row=[0]*(numcol)
        Backtrack.append(row)

    #Initialization
    score[0][0]=0
    for c in range(1,numcol):
        score[0][c]=score[0][c-1]-indel_penalty
    
    for r in range(1, numrow):
        for c in range(1,numcol):
            if s[c-1]==t[r-1]:
                diagonal_score=match_reward
            else:
                diagonal_score=-mismatch_penalty
            score[r][c]=max(score[r-1][c-1]+diagonal_score,score[r-1][c]-indel_penalty,score[r][c-1]-indel_penalty)
            if score[r][c]==score[r-1][c]-indel_penalty:
                Backtrack[r][c]="↓"
            elif score[r][c]==score[r][c-1]-indel_penalty:
                Backtrack[r][c]="→"
            elif score[r][c]==score[r-1][c-1]+diagonal_score:
                Backtrack[r][c]="↘"
            else:
                Backtrack[r][c]="X"
                
    maxcol=np.argmax(score[numrow-1])
    maxscore=score[numrow-1][maxcol]
    s_align,t_align=AlignBacktrackLocal(Backtrack,s,t,numrow-1,maxcol)

    return (maxscore,t_align,s_align)

'''
Solve the Alignment with Affine Gap Penalties Problem.

Input: A match reward, a mismatch penalty, a gap opening penalty, a gap extension penalty, and two nucleotide strings.
Output: The maximum alignment score between v and w, followed by an alignment of v and w achieving this maximum score.
'''
def AffineAlignment(match_reward: int, mismatch_penalty: int, gap_opening_penalty: int, gap_extension_penalty: int, s: str, t: str) -> tuple[int, str, str]:
    numcol=len(s)+1
    numrow=len(t)+1

    #Construct the 3-layer Manhattan graph
    middle=[]
    for r in range(numrow):
        row=[0]*(numcol)
        middle.append(row)

    lower=[]
    for r in range(numrow):
        row=[0]*(numcol)
        lower.append(row)
    for c in range(numcol):
        lower[0][c]=-math.inf

   
    upper=[]
    for r in range(numrow):
        row=[0]*(numcol)
        upper.append(row)
    for r in range(numrow):
        upper[r][0]=-math.inf

    Backtrack=[]
    for r in range(numrow):
        row=[0]*(numcol)
        Backtrack.append(row)
    

    for r in range(numrow):
        for c in range(numcol):
            if r == 0 and c ==0:
                continue

            if c>0:
                upper[r][c]=max(upper[r][c-1]-gap_extension_penalty,middle[r][c-1]-gap_opening_penalty)
            if r>0:
                lower[r][c]=max(lower[r-1][c]-gap_extension_penalty,middle[r-1][c]-gap_opening_penalty)
            if r==0:
                middle[r][c]=upper[r][c]
            elif c==0:
                middle[r][c]=lower[r][c]
            else:
                if s[c-1]==t[r-1]:
                    diagonal_score=match_reward
                else:
                    diagonal_score=-mismatch_penalty

                middle[r][c]=max(middle[r-1][c-1]+diagonal_score,lower[r][c],upper[r][c])

            if middle[r][c]==lower[r][c]:
                Backtrack[r][c]="↓"
            elif middle[r][c]==upper[r][c]:
                Backtrack[r][c]="→"
            else:
                Backtrack[r][c]="↘"

    s_align,t_align=AlignBacktrackAffine(Backtrack,lower,upper,s,t,numrow-1,numcol-1,gap_extension_penalty)

    return (middle[numrow-1][numcol-1],s_align,t_align)


#Affine version of backtrack
def AlignBacktrackAffine(track: list[list[str]],lower: list[list[int]],upper: list[list[int]], s: str, t: str, r: int, c: int,gap_ex: int) -> tuple[str,str]:
    if r==0:
        return(s[0:c],'-'*c)
    if c==0:
        return('-'*r,t[0:r])
    if track[r][c]=="↘":
        s_align,t_align=AlignBacktrackAffine(track,lower,upper,s,t,r-1,c-1,gap_ex)
        return(s_align+s[c-1],t_align+t[r-1])
    
    # When it comes to iindels, check for the length of the gap
    elif track[r][c]=="→": 
        i=0     
        while upper[r][c] == upper[r][c-1]-gap_ex:
            c-=1
            i+=1
        s_align,t_align=AlignBacktrackAffine(track,lower,upper,s,t,r,c-1,gap_ex)
        return(s_align+s[c-1:c+i],t_align+'-'*(i+1))
    else:
        i=0     
        while lower[r][c] == lower[r-1][c]-gap_ex:
            r-=1
            i+=1
        s_align,t_align=AlignBacktrackAffine(track,lower,upper,s,t,r-1,c,gap_ex)
        return(s_align+'-'*(i+1),t_align+t[r-1:r+i])
            
'''
Solve the Multiple Longest Common Subsequence Problem.

Input: Three DNA strings of length at most 10.
Output: The length of a longest common subsequence of these three strings, followed by a multiple alignment of the three strings corresponding to such an alignment.
'''
def MultipleAlignment(s1: str, s2: str, s3: str) -> tuple[int, str, str, str]:
    numd1=len(s1)+1
    numd2=len(s2)+1
    numd3=len(s3)+1

    # establish the matrixes
    score=[]
    for i in range(numd1):
        score2=[]
        score.append(score2)
        for j in range(numd2):
           score3 =[0]*(numd3)
           score2.append(score3)

    track=[]
    for i in range(numd1):
        track2=[]
        track.append(track2)
        for j in range(numd2):
           track3 =[0]*(numd3)
           track2.append(track3)

    # Fill in the 3D matrix
    for i in range(1,numd1):
        for j in range(1,numd2):
            for k in range(1,numd3):
                match=0
                if s1[i-1]==s2[j-1] and s2[j-1]==s3[k-1]:
                    match=1
                score[i][j][k]=max(score[i-1][j][k],score[i][j-1][k],score[i][j][k-1],score[i-1][j-1][k],score[i-1][j][k-1],score[i][j-1][k-1],score[i-1][j-1][k-1]+match)
                
                if score[i][j][k]==score[i-1][j-1][k-1]+match:
                    track[i][j][k]="↓""↓""↓"
                elif score[i][j][k]==score[i-1][j][k]:
                    track[i][j][k]="↓""-""-"
                elif score[i][j][k]==score[i][j-1][k]:
                    track[i][j][k]="-""↓""-"
                elif score[i][j][k]==score[i-1][j-1][k]:
                    track[i][j][k]="↓""↓""-"
                elif score[i][j][k]==score[i-1][j][k-1]:
                    track[i][j][k]="↓""-""↓"
                elif score[i][j][k]==score[i][j-1][k-1]:
                    track[i][j][k]="-""↓""↓"
                else:
                    track[i][j][k]="-""-""↓"

    s1_align, s2_align,s3_align=AlignBacktrack3D(track,s1,s2,s3,numd1-1,numd2-1,numd3-1)

    return (score[numd1-1][numd2-1][numd3-1],s1_align,s2_align,s3_align)

# Get the track of the 3D alignment that yields the highest score
def AlignBacktrack3D(track:list[list[list[str]]],s1:str,s2:str,s3:str,i:int,j:int,k:int)-> tuple[str,str,str]:
    #boundary conditions
    if i ==0:
        if j>k:
            return('-'*j,s2[:j],s3[:k]+'-'*(j-k))
        else:
            return('-'*k,s2[:j]+'-'*(k-j),s3[:k])
    if j ==0:
        if i>k:
            return(s1[:i],'-'*i,s3[:k]+'-'*(i-k))   
        else:
            return(s1[:i]+'-'*(k-i),'-'*k,s3[:k])
    if k ==0:
        if i>j:
            return(s1[:i],s2[:j]+'-'*(i-j),'-'*i)
        else:
            return(s1[:i]+'-'*(j-i),s2[:j],'-'*j)
    
    # For each direction, if the arrow moves in that direction, the corresponding string has a letter aligned, else it's an indel.
    if track[i][j][k][0]=="↓":
        ad1=s1[i-1]
        i2=i-1
    else:
        ad1='-'
        i2=i
    if track[i][j][k][1]=="↓":
        ad2=s2[j-1]
        j2=j-1
    else:
        ad2='-'
        j2=j
    if track[i][j][k][2]=="↓":
        ad3=s3[k-1]
        k2=k-1
    else:
        ad3='-'
        k2=k

    align1,align2,align3=AlignBacktrack3D(track,s1,s2,s3,i2,j2,k2)

    return (align1+ad1,align2+ad2,align3+ad3)

