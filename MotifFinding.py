import random
import numpy


'''
RandomizedMotifSearch implements an algorithm to find a collection of motifs length of k among t strings of DNA.
By taking a random collection of motifs as a start, the algorithm generates a profile matrix based on the motifs
and relocate the motifs that fits the profile the most in each string of DNA.
This process is run for 1000 times to obtain the best motifs with the lowest score. 
Note: Lower score indicates better motifs.
'''
def RandomizedMotifSearch(Dna: list[str], k: int, t: int) -> list[str]:
    # Use the first set of generated motifs as the initial BestMotifs.
    BestMotifs=RandomizedMotifSearchSingle(Dna,k,t)
    BestScore=Score(BestMotifs)

    for i in range(999):
        new_motifs=RandomizedMotifSearchSingle(Dna,k,t)
        new_score=Score(motifs)
        # If better motifs are generated, update BestMotifs and BestScore
        if new_score<BestScore:
            BestMotifs=new_motifs
            BestScore=new_score
    
    return BestMotifs


'''
RandomizedMotifSearchSingle is a single run of the RandomizedMotifSearch algorithm where
the process of Profile(Motif)-Motif(Profile(Motif),Dna)-Profile(Motif(Profile)) keeps iterating 
until the score of the motifs stops decreasing.
'''
def RandomizedMotifSearchSingle(Dna: list[str], k: int, t: int) -> list[str]:
    # Generate random motifs as a start
    new_motifs=GenerateRandomMotifs(Dna,k,t)
    BestMotifs=new_motifs
    BestScore=Score(BestMotifs)

    forever=True
    # time=0
    # time is the number of times this random collection of motifs gets updated.
    # print(time) before return BestMotifs to see results

    while forever:
        # time+=1
        # Generate the profile of the current motifs.
        new_profile=MakeProfile(new_motifs)
        # Based on the current profile, choose the new motifs. 
        new_motifs=MostProbableMotifs(new_profile,k,Dna)
        # Get the score of the new motifs.
        new_score=Score(new_motifs)

        # if the score is better(lower), update the BestMotifs.
        if new_score<BestScore: 
            BestMotifs=new_motifs
            BestScore=new_score
        else:
        # if the score stops decreasing, end the while loop and return the current BestMotifs.
            #print(time)
            return BestMotifs

'''
This function generates a random collection of motifs length of k from t strings of DNA. 
'''
def GenerateRandomMotifs(Dna: list[str], k: int, t: int) -> list[str]:
    motifs=[] 

    # For each DNA sequence, choose a motif at a random position.
    for i in range(t):
        # n is the length of the DNA sequence
        n=len(Dna[i])
        # There are n-k+1 options of DNA sequence start index.
        index=random.randint(0,n-k)
        mer=Dna[i][index:index+k]
        motifs.append(mer)

    return motifs

'''
MakeProfile generates a frequency profile based on a collection of motifs.
'''
def MakeProfile(Motifs: list[str]) -> list[dict[str,float]]:
    k=len(Motifs[0])
    # Profile is a list of k dictionaries, where k is the length of each motif.
    # Each dictionary contains the frequency of each nucleotide base occurring at the i-th position of the motif.
    profile=[]
    # t is the number of motifs(DNA sequences).
    t=len(Motifs)
    # Each nucleotide base has one pseudo count at each position, thus at each position the total number of nucleotide bases is t+4.
    n=t+4

    for i in range(k):
        freqA=(CountBase('A',Motifs,i)+1)/n
        freqC=(CountBase('C',Motifs,i)+1)/n
        freqG=(CountBase('G',Motifs,i)+1)/n
        freqT=(CountBase('T',Motifs,i)+1)/n
        dict={'A':freqA,'C':freqC,'G':freqG,'T':freqT}
        profile.append(dict)

    return profile

'''
CountBase counts the occurance of Base at the pos-th position of the motifs.
'''
def CountBase(Base:str, Motifs:list[str], pos:int) -> int:
    count=0

    # range through each motif to check if Base occurred.
    for i in range(len(Motifs)):
        if Motifs[i][pos]==Base:
            count+=1
    
    return count

'''
MostProbableMotifs finds the set of motifs length of k in a set of DNA sequences with the highest probability of fitting given profile.
'''
def MostProbableMotifs(Profile: list[dict[str, float]], k: int, Dna: list[str]) -> list[str]:
    motifs=[]

    # range through each DNA string to find the best motif in each sequence
    for i in range(len(Dna)):
        kmer=ProfileMostProbableKmer(Dna[i],k,Profile)
        motifs.append(kmer)

    return motifs

'''
ProfileMostProbableKmer takes a string Text(a DNA sequence), k the length of k-mer and a profile with k dictionaries as input 
and outputs a k-mer in Text that has the higest probability to occur according to the profile.
'''
def ProfileMostProbableKmer(Text: str, k: int, Profile: list[dict[str, float]]) -> str:
    # set the initial highest probability as 0 since all probabilies are equal or greater than 0.
    highest_prob=0.0
    # range through every starting positions of k-mer in the Text string
    for i in range(len(Text)-k+1):
        mer=Text[i:i+k]
        prob=Prob(mer,Profile)
        # update the opt_mer if probability gets higher.
        if prob > highest_prob:
            highest_prob=prob
            opt_mer=mer

    return opt_mer

'''
Prob calculates the probability of a string mer occurring according to the profile.
'''
def Prob(mer: str, Profile: list[dict[str, float]]) -> float:
    score=1.0

    # probability equals to the probability of each base in the mer occuring at each position multiplied together.
    # i.e. P('ATG') = P(mer_1='A') x P(mer_2='T") x P(mer_3='G')
    for i in range(len(mer)):
        score*=Profile[i][mer[i]]

    return score

'''
Score is a scoring function of the motifs used to quantify the goodness of a collection of k-mers.
Score is defined by the number of bases that's different from the most frequently appearred base at each position added together.
'''       
def Score(Motifs: list[str]) -> int:
    # t is the number of motifs.
    t=len(Motifs)
    score=0.0
    # For each position of the mers, count the occurance of each base.
    for i in range(len(Motifs[0])):
        numA=CountBase('A',Motifs,i)
        numT=CountBase('T',Motifs,i)
        numC=CountBase('C',Motifs,i)
        numG=CountBase('G',Motifs,i)

        MostBase=numA
        if numT>MostBase:
            MostBase=numT
        if numC>MostBase:
            MostBase=numC
        if numG>MostBase:
            MostBase=numG
        
        # MostBase is the number of the base that appearred the most.
        # t-MostBase is the number of bases that didn't agree with the Most Base.
        score+=(t-MostBase)
    
    return score

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
'''
GibbsSampler is another algorithm that aims at finding motifs across a number of DNA sequences.
GibbsSampler changes only one k-mer at one step following the probability distribution of each k-mer in a randomly chosed sequence.
GibbsSampler here runs GibbsSampler for 20 times(with 20 different random start) to circumvent the impact of a possibly terrible starting point. 
This number could be altered with respect to the search space.
'''
def GibbsSampler(Dna: list[str], k: int, t: int, N: int) -> list[str] :
    BestMotifs=GibbsSamplerSingle(Dna,k,t,N)
    BestScore=Score(BestMotifs)

    for i in range(19):
        new_motifs=GibbsSamplerSingle(Dna,k,t,N)
        new_score=Score(new_motifs)
        # if score gets lower, update the BestMotifs.
        if new_score<BestScore:
            BestMotifs=new_motifs
            BestScore=new_score
  
    return BestMotifs


def GibbsSamplerSingle(Dna: list[str], k: int, t: int, N: int) -> list[str]:
    new_motifs=GenerateRandomMotifs(Dna,k,t)
    BestMotifs=Copy(new_motifs)
    BestScore=Score(BestMotifs)

    # keep iterating the motifs for N steps.
    for j in range(N): # in each iteration
        # First, randomly choose a k-mer to be removed
        i=random.randint(0,t-1)
        new_motifs.pop(i)
        # Generate the profile based on the remaining k-mers.
        profile=MakeProfile(new_motifs)
        # From the sequence where the corresponding k-mer is deleted, choose a new k-mer that fits the profile according to the probability distribution.
        mer=ProfileRandomKmer(Dna[i],k,profile)
        new_motifs.insert(i,mer)
        # Get the score of the updated motifs.
        new_score=Score(new_motifs)

        # Update the BestMotifs and BestScore if the score gets better
        if new_score<BestScore:
            BestMotifs=Copy(new_motifs)
            BestScore=new_score
    
    return BestMotifs

'''
ProfileRandomKmer selects a k-mer from a string according to the probability distribution of each k-mer given the profile.
'''
def ProfileRandomKmer(Text: str, k: int, Profile: list[dict[str, float]]) -> str:
    # n is the number of start indices that can be chosen from the Text string
    n=len(Text)-k+1
    # Problist is a probability distribution list of all possible k-mers in the Text string
    ProbList=[]
    # sum is the sum of all the probabilities used for normalizing the probability distribution
    sum=0.0

    # for each possible k-mer
    for i in range(n):
        mer=Text[i:i+k]
        # acquire the probability of this k-mer occuring according to the given profile.
        prob=Prob(mer,Profile)
        ProbList.append(prob)
        sum+=prob
    
    # Normalize the probabilities so that the probabilities sum up to 1.00
    for i in range(n):
        ProbList[i]=ProbList[i]/sum
    
    # Following the probability distribution, choose a k-mer 
    index=numpy.random.choice(n,1,p=ProbList)[0]
    return(Text[index:index+k])

'''
Deep copy a list
'''
def Copy(listB: list[str]) -> list[str]:
    listA=[]

    for i in range(len(listB))  :
        listA.append(listB[i])

    return listA


#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------

'''Running GibbsSampler from a given dataset.
with open('dataset_647149_9.txt','r') as f:
    lines=f.readlines()
    Dna=[]
    for i in range(len(lines)):
        lines[i]=lines[i].rstrip('\n')
        if i ==0:
            lines[i]=list(lines[i].split(' '))
        else:
            Dna.append(lines[i])
    
    k=int(lines[0][0])
    t=int(lines[0][1])
    N=int(lines[0][2])
    Motifs=GibbsSampler(Dna,k,t,N)

with open('answer.txt','w') as f:
    for i in range(len(Motifs)):
        Motifs[i]+='\n'
    f.writelines(Motifs)
'''

'''Running RandomizedMotifSearch from a given dataset.
with open('dataset_647148_7.txt','r') as f:
    lines=f.readlines()
    Dna=[]
    for i in range(len(lines)):
        lines[i]=lines[i].rstrip('\n')
        if i ==0:
            lines[i]=list(lines[i].split(' '))
        else:
            Dna.append(lines[i])
    
    k=int(lines[0][0])
    t=int(lines[0][1])
    Motifs=RandomizedMotifSearch(Dna,k,t)

with open('answer.txt','w') as f:
    for i in range(len(Motifs)):
        Motifs[i]+='\n'
    f.writelines(Motifs)
'''




