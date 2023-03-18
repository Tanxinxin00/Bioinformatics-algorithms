'''
This is Week 5 code submission for 02604 Bioinformatics from Tanxin(Thea) Qiao.
Written in Python=3.9.13 
03/03/2023
'''
import numpy as np

'''
Solve the Probability of a Hidden Path Problem.

Input: A hidden path π followed by the states States and transition matrix Transition of an HMM (Σ, States, Transition, Emission).
Output: The probability of this path, Pr(π).
'''
def ProbHiddenPath(path: str, States:list[str], Transition:dict[dict[float]] )->float:
    n=len(States)
    # The initial probability is evenly distributed
    P=1/n

    for i in range (len(path)-1):
        P=P*Transition[path[i]][path[i+1]]

    return P


'''
Solve the Probability of an Outcome Given a Hidden Path Problem.

Input: A string x, followed by the alphabet from which x was constructed, followed by a hidden path π, followed by the states States and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
Output: The conditional probability Pr(x|π) that x will be emitted given that the HMM follows the hidden path π.
'''
def ProbOfEmission(emit: str, Alphabet: list[str], Hiddenpath: str, States: list[str], EmissionMtx: dict[str:dict[str]])-> float:
    P=1.0

    # each emission probability is independent of others
    for i in range(len(Hiddenpath)):
        P*=EmissionMtx[Hiddenpath[i]][emit[i]]
    
    return P


'''
Implement the Viterbi algorithm solving the Decoding Problem.

Input: A string x, followed by the alphabet from which x was constructed, followed by the states States, transition matrix Transition, and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
Output: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.
'''
def Viterbi(emit: str, Alphabet: list[str], States: list[str], Transition: dict[str:dict[str]],EmissionMtx: dict[str:dict[str]])-> str:
    n=len(emit)
    s=len(States)

    # create matrix weights & path
    weights=[]
    for i in range(s):
        row=[1.0]*n
        weights.append(row)
    
    path=[]
    for i in range(s):
        row=['']*n
        path.append(row)

    # Fill in the first column of the matrixes
    for r in range(s):
        weights[r][0]=EmissionMtx[States[r]][emit[0]]
        path[r][0]+=States[r]
    
    # Fill in the rest of the matrixes sequentially
    for c in range(1,n):
        for r in range(s):

            maximum=0
            probemit=EmissionMtx[States[r]][emit[c]]

            # Find the former node with the biggest weight
            for i in range(s):
                weight=weights[i][c-1]*Transition[States[i]][States[r]]*probemit
                if weight>maximum:
                    maximum=weight
                    index=i

            weights[r][c]=maximum
            path[r][c]=path[index][c-1]+States[r]
    
    #locate the maximum value in the last column
    index=np.argmax([weights[r][n-1] for r in range(s)])

    return path[index][n-1]

'''
Solve the Profile HMM Problem.

Input: A threshold θ, followed by an alphabet Σ, followed by a multiple alignment Alignment whose strings are formed from Σ.
Output: The transition matrix followed by the emission matrix of HMM(Alignment, θ).
'''
def ProfileHMM(theta:float, alphabet:list[str],Alignment:list[str]):
    n=len(Alignment)

    ab_index=[]
    for i in range(len(Alignment[0])):
        count=0
        for seq in Alignment:
            if seq[i]=='-':
                count+=1

        if count/n > theta:
            ab_index.append(i)
    
    seq_len=len(Alignment[0])-len(ab_index)

    states=Makestates(seq_len)
    Transition=MakeTransitionMtx(states)
    Emission=MakeEmissionMtx(states,alphabet)

    # count the emissions for each state and symbol
    for seq in Alignment:

        currentstate='S'
        index=0
        for i in range(len(seq)):

            symbol=seq[i]
 
            if i not in ab_index:
                index+=1
                if symbol=='-':
                    nextstate=f'D{index}'
                    Transition[currentstate][nextstate]+=1
                    currentstate=nextstate
                else:
                    nextstate=f'M{index}'
                    Transition[currentstate][nextstate]+=1
                    Emission[nextstate][symbol]+=1
                    currentstate=nextstate
            elif symbol!='-':
                nextstate=f'I{index}'
                Transition[currentstate][nextstate]+=1
                Emission[nextstate][symbol]+=1
                currentstate=nextstate
        Transition[currentstate]['E']+=1   

    return(FreqToProb(Transition),FreqToProb(Emission))

# Create a list of states as needed
def Makestates(n:int)->list[str]:
    states=['S','I0']
    for i in range(n):
        states.append(f'M{i+1}')
        states.append(f'D{i+1}')
        states.append(f'I{i+1}')
    states.append('E')
    return states

# make a fully connected transition matx given a list of states
def MakeTransitionMtx(states:list[str]):
    transition={}
    for state in states:
        transition[state]={}
        for state2 in states:
            transition[state][state2]=0.0
    return transition

# make an emission mtx given a list of states and all possible emissions
def MakeEmissionMtx(states:list[str],alphabet:list[str]):
    emission={}
    for state in states:
        emission[state]={}
        for symbol in alphabet:
            emission[state][symbol]=0.0
    return emission

# Transfer a dict of frequencies into a dict of probabilities
def FreqToProb(matrix:dict[str,dict[str,float]]):

    for key in matrix.keys():
        sum=0
        for key2 in matrix[key].keys():
            sum+=matrix[key][key2]
        for key2 in matrix[key].keys():
            if sum!=0:
                matrix[key][key2]=round(matrix[key][key2]/sum,3)
    
    return matrix

'''
Solve the Profile HMM with Pseudocounts Problem.

Input: A threshold θ and a pseudocount σ, followed by an alphabet Σ, followed by a multiple alignment Alignment whose strings are formed from Σ.
Output: The transition and emission matrices of HMM(Alignment, θ, σ).
'''
def ProfileHMMPseudo(theta:float, alphabet:list[str],Alignment:list[str],pseudo:float):
    n=len(Alignment)

    ab_index=[]
    for i in range(len(Alignment[0])):
        count=0
        for seq in Alignment:
            if seq[i]=='-':
                count+=1
        if count/n > theta:
            ab_index.append(i)
    
    seq_len=len(Alignment[0])-len(ab_index)

    states=Makestates(seq_len)
    Transition=MakeTransitionMtx(states)
    Emission=MakeEmissionMtx(states,alphabet)

    for seq in Alignment:

        currentstate='S'
        index=0
        for i in range(len(seq)):

            symbol=seq[i]
 
            if i not in ab_index:
                index+=1
                if symbol=='-':
                    nextstate=f'D{index}'
                    Transition[currentstate][nextstate]+=1
                    currentstate=nextstate
                else:
                    nextstate=f'M{index}'
                    Transition[currentstate][nextstate]+=1
                    Emission[nextstate][symbol]+=1
                    currentstate=nextstate
            elif symbol!='-':
                nextstate=f'I{index}'
                Transition[currentstate][nextstate]+=1
                Emission[nextstate][symbol]+=1
                currentstate=nextstate
        Transition[currentstate]['E']+=1   
    Transition=FreqToProb(Transition)
    Emission=FreqToProb(Emission)

    Emission=PseudoEmission(Emission,pseudo)
    Transition=PseudoTransition(Transition,seq_len,pseudo)

    return(states, Transition,Emission)

# add pseudo counts to the transition mtx
def PseudoTransition(TransitionMtx,n,pseudo):
    for i in range(n+1):
        if i ==0:
            prev=['S',f'I{i}']
            latter=[f'I{i}',f'M{i+1}',f'D{i+1}']
        elif i ==n:
            prev=[f'I{i}',f'M{i}',f'D{i}']
            latter=['E',f'I{i}']
        else:
            prev=[f'I{i}',f'M{i}',f'D{i}']
            latter=[f'I{i}',f'M{i+1}',f'D{i+1}']

        for key in prev:
            for key2 in latter:
                if TransitionMtx[key][key2]==0:
                    TransitionMtx[key][key2]+=pseudo
            total=sum(TransitionMtx[key].values())
            for key2 in latter:
                TransitionMtx[key][key2]=round(TransitionMtx[key][key2]/total,3)
            
    return TransitionMtx

#add pseudo counts to the emission mtx
def PseudoEmission(Emission,pseudo):
    for key in Emission.keys():
        if key[0]=='I'or key[0]=='M':
            total=sum(Emission[key].values())
            for key2 in Emission[key].keys():
                if Emission[key][key2]==0:
                    Emission[key][key2]+=pseudo
                    total+=pseudo
            for key2 in Emission[key].keys():
                    Emission[key][key2]=round(Emission[key][key2]/total,3)
    return Emission


'''
Solve the Sequence Alignment with Profile HMM Problem.

Input: A string x followed by a threshold θ and a pseudocount σ, followed by an alphabet Σ, followed by a multiple alignment Alignment whose strings are formed from Σ.
Output: An optimal hidden path emitting x in HMM(Alignment, θ, σ).
'''
def SequenceAlignHMM(x:str, theta: float, pseudo:float, alphabet: list[str], Align:list[str]):
    n,states,Transition,Emission= ProfileHMMPseudo(theta,alphabet, Align, pseudo)


    numrow=len(states)-2
    numcol=len(x)

    #create matrix
    weights=[]
    for i in range(numrow):
        row=[1.0]*numcol
        weights.append(row)
    
    path=[]
    for i in range(numrow):
        row=['']*numcol
        path.append(row)

    #initialization
    col_int=[1.0]*(n+1)
    for i in range(n+1):
        if i==1:
            col_int[i]=Transition['S']['D1']
        elif i>1:
            col_int[i]=col_int[i-1]*Transition[f'D{i-1}'][f'D{i}']

    # first fill in the first column
    for i in range(n+1):
        if i==0:
            weights[i][0]=Transition['S'][f'I{i}']*Emission[f'I{i}'][x[0]]
            path[i][0]=f'I{i}'

            weights[i+1][0]=Transition['S'][f'M{i+1}']*Emission[f'M{i+1}'][x[0]]
            path[i+1][0]=f'M{i+1}'

            weights[i+2][0]=weights[i][0]*Transition[f'I{i}'][f'D{i+1}']
            path[i+2][0]=path[i][0]+f' D{i+1}'

        elif i <n:
            k=3*i
            weights[k][0]=col_int[i]*Transition[f'D{i}'][f'I{i}']*Emission[f'I{i}'][x[0]]
            for j in range(1,i+1):
                path[k][0]+=f' D{j}'
            path[k][0]+=f' I{i}'

            weights[k+1][0]=col_int[i]*Transition[f'D{i}'][f'M{i+1}']*Emission[f'M{i+1}'][x[0]]
            for j in range(1,i+1):
                path[k+1][0]+=f' D{j}'
            path[k+1][0]+=f' M{i+1}'

            I=weights[k][0]*Transition[f'I{i}'][f'D{i+1}']
            D=weights[k-1][0]*Transition[f'D{i}'][f'D{i+1}']
            M=weights[k-2][0]*Transition[f'M{i}'][f'D{i+1}']
            weights[k+2][0]=max(I,D,M)
            if weights[k+2][0]==I:
                path[k+2][0]=path[k][0]+f' D{i+1}'
            elif weights[k+2][0]==D:
                path[k+2][0]=path[k-1][0]+f' D{i+1}'
            else:
                path[k+2][0]=path[k-2][0]+f' D{i+1}'

        else:
            k=3*i
            weights[k][0]=col_int[i]*Transition[f'D{i}'][f'I{i}']*Emission[f'I{i}'][x[0]]
            for j in range(1,i+1):
                path[k][0]+=f' D{j}'
            path[k][0]+=f' I{i}'

    # then fill in the rest columns sequentially 
    for c in range(1,numcol):

        for i in range(n+1):
            if i==0:
                weights[i][c]=weights[i][c-1]*Transition[f'I{i}'][f'I{i}']*Emission[f'I{i}'][x[c]]
                path[i][c]=path[i][c-1]+f' I{i}'

                weights[i+1][c]=weights[i][c-1]*Transition[f'I{i}'][f'M{i+1}']*Emission[f'M{i+1}'][x[c]]
                path[i+1][c]=path[i][c-1]+f' M{i+1}'

                weights[i+2][c]=weights[i][c]*Transition[f'I{i}'][f'D{i+1}']
                path[i+2][c]=path[i][c]+f' D{i+1}'

            elif i <n:
                k=3*i
                D=weights[k-1][c-1]*Transition[f'D{i}'][f'I{i}']
                I=weights[k][c-1]*Transition[f'I{i}'][f'I{i}']
                M=weights[k-2][c-1]*Transition[f'M{i}'][f'I{i}']
                weights[k][c]=max(I,D,M)
                if weights[k][c]==I:
                    path[k][c]=path[k][c-1]+f' I{i}'
                elif weights[k][c]==D:
                    path[k][c]=path[k-1][c-1]+f' I{i}'
                else:
                    path[k][c]=path[k-2][c-1]+f' I{i}'
                weights[k][c]*=Emission[f'I{i}'][x[c]]

                I=weights[k][c-1]*Transition[f'I{i}'][f'M{i+1}']
                D=weights[k-1][c-1]*Transition[f'D{i}'][f'M{i+1}']
                M=weights[k-2][c-1]*Transition[f'M{i}'][f'M{i+1}']
                weights[k+1][c]=max(I,D,M)
                if weights[k+1][c]==I:
                    path[k+1][c]=path[k][c-1]+f' M{i+1}'
                elif weights[k+1][c]==D:
                    path[k+1][c]=path[k-1][c-1]+f' M{i+1}'
                else:
                    path[k+1][c]=path[k-2][c-1]+f' M{i+1}'
                weights[k+1][c]*=Emission[f'M{i+1}'][x[c]]

                I=weights[k][c]*Transition[f'I{i}'][f'D{i+1}']
                D=weights[k-1][c]*Transition[f'D{i}'][f'D{i+1}']
                M=weights[k-2][c]*Transition[f'M{i}'][f'D{i+1}']
                weights[k+2][c]=max(I,D,M)
                if weights[k+2][c]==I:
                    path[k+2][c]=path[k][c]+f' D{i+1}'
                elif weights[k+2][c]==D:
                    path[k+2][c]=path[k-1][c]+f' D{i+1}'
                else:
                    path[k+2][c]=path[k-2][c]+f' D{i+1}'

            else:
                k=3*i
                M=weights[k-2][c-1]*Transition[f'M{i}'][f'I{i}']
                D=weights[k-1][c-1]*Transition[f'D{i}'][f'I{i}']
                I=weights[k][c-1]*Transition[f'I{i}'][f'I{i}']
                weights[k][c]=max(I,D,M)
                if weights[k][c]==I:
                    path[k][c]=path[k][c-1]+f' I{i}'
                elif weights[k][c]==D:
                    path[k][c]=path[k-1][c-1]+f' I{i}'
                else:
                    path[k][c]=path[k-2][c-1]+f' I{i}'
                weights[k][c]*=Emission[f'I{i}'][x[c]]

    I=weights[numrow-1][numcol-1]*Transition[f'I{n}']['E']
    D=weights[numrow-2][numcol-1]*Transition[f'D{n}']['E']
    M=weights[numrow-3][numcol-1]*Transition[f'M{n}']['E']
    finalscore=max(I,D,M)
    if finalscore==I:
        return path[numrow-1][numcol-1]
    elif finalscore==D:
        return path[numrow-2][numcol-1]
    else:
        return path[numrow-3][numcol-1]
