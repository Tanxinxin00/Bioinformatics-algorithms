# The amino acid encoding codon dictionary
AAdict={'AAA':'K','AAC': 'N','AAG': 'K','AAU': 'N','ACA': 'T','ACC' :'T','ACG': 'T','ACU': 'T','AGA': 'R',
        'AGC':'S','AGG': 'R','AGU':'S','AUA': 'I','AUC': 'I','AUG' :'M','AUU': 'I','CAA': 'Q','CAC': 'H',
        'CAG':'Q','CAU': 'H','CCA': 'P','CCC': 'P','CCG': 'P','CCU':'P','CGA': 'R','CGC':'R','CGG':'R',
        'CGU':'R','CUA': 'L','CUC': 'L','CUG': 'L','CUU': 'L','GAA':'E','GAC':'D','GAG':'E','GAU': 'D',
        'GCA': 'A','GCC': 'A','GCG': 'A','GCU': 'A','GGA': 'G','GGC':'G','GGG':'G','GGU':'G','GUA': 'V',
        'GUC': 'V','GUG': 'V','GUU': 'V','UAA':'*', 'UAC': 'Y','UAG':'*','UAU': 'Y','UCA': 'S','UCC':'S',
        'UCG': 'S','UCU': 'S','UGA':'*' ,'UGC': 'C','UGG' :'W','UGU': 'C','UUA': 'L','UUC': 'F','UUG': 'L',
        'UUU': 'F'}

'''
PeptideEncode takes a string of DNA sequence Text, a string of amino acid sequence Peptide and the codon dictionary as input
and outputs every possible substring of the DNA sequence that could encode the peptide, including reverse complement substrings.
'''
def PeptideEncode(Text: str, Peptide: str, AAdict: dict[str:str]) -> list[str]:
    # Convert the 'T's in DNA sequence to 'U's as in mRNA condons
    Converted=ConvertTtoU(Text)
    # Do the same for the reverse complement string
    RCText=ReverseComplement(Text)  
    RCConverted=ConvertTtoU(RCText)    

    k=3*len(Peptide)
    Substrings=[]
    n=len(Text)

    # range through alll the possible substrings of the converted sequence
    for i in range(n-k+1):
        Substring=Converted[i:i+k]
        Match=True
        # check if the substring could encode the peptide
        for j in range (len(Peptide)):          
            if AAdict[Substring[3*j:3*j+3]]!=Peptide[j]:
                Match = False
                break
        if Match:
            Substrings.append(Text[i:i+k])
    
    # range through alll the possible substrings of the converted reverse complement sequence
    for i in range(n-k+1):
        Substring=RCConverted[i:i+k]
        Match=True
        # check if the substring could encode the peptide
        for j in range (len(Peptide)):          
            if AAdict[Substring[3*j:3*j+3]]!=Peptide[j]:
                Match = False
                break
        if Match:
            Substrings.append(Text[n-i-k:n-i])  

    return Substrings

'''
Converts all the 'T's into 'U's 
'''
def ConvertTtoU(Text: str) -> str:
    ConvertedText=''
    for i in range(len(Text)):
        if Text[i]=='T':
            ConvertedText=ConvertedText+'U'
        else:
            ConvertedText=ConvertedText+Text[i]
    return ConvertedText

'''
Get the reverse complement sequence of a given sequence Text.
'''
def ReverseComplement(Text: str) -> str:
    n=len(Text)
    RCText=''
    for i in range(n):
        if Text[n-1-i]=='T':
            RCText=RCText+'A'
        elif Text[n-1-i]=='A':
            RCText=RCText+'T'
        elif Text[n-1-i]=='C':
            RCText=RCText+'G'
        elif Text[n-1-i]=='G':
            RCText=RCText+'C'

    return RCText

# The amino acid mass dictionary
AAmass={'G': 57,'A': 71,'S': 87,'P': 97,'V': 99,'T': 101,'C': 103,'I': 113,'L': 113,'N': 114,'D': 115,'K': 128,
        'Q' :128,'E': 129,'M' :131,'H' :137,'F': 147,'R': 156,'Y': 163,'W':186}

'''
Cyclospectrum takes a string of amino acid Peptide as input and 
outputs the theoretical spectrum of the peptide, i.e., the masses of all its subpeptides, including 0 and its whole mass
'''
def CycloSpectrum(Peptide:str) -> list[int]:  
    # Get the cyclo subpeptides of the peptide
    cyclosubs=CycloSubPep(Peptide)
    # append itself to it
    cyclosubs.append(Peptide)

    # construct the spectrum of the peptide
    spectrum=[]
    spectrum.append(0)
    for sub in cyclosubs:
        spectrum.append(CalculateMass(sub))
    spectrum.sort()
    return(spectrum)


'''
Linearspectrum takes a string of amino acid Peptide as input and 
outputs the theoretical linear spectrum of the peptide, i.e., the masses of all its linearsubpeptides, including 0 and its whole mass, 
'''
def LinearSpectrum(Peptide:str) -> list[int]:
    linearsubs=LinearSubPep(Peptide)
    linearsubs.append(Peptide)
    
    spectrum=[]
    spectrum.append(0)
    for sub in linearsubs:
        spectrum.append(CalculateMass(sub))
    spectrum.sort()

    return(spectrum)

'''
CycloSubPep takes a string of amino acid Peptide as input
and outputs a list of cyclosubpeptides of Peptide, including the wrap-arounds.
'''
def CycloSubPep(Peptide:str) -> list[str]:
    n=len(Peptide)
    substrings=[]
    for k in range (1, n, 1):
        for i in range (n):
            if i+k-1<n:
                substring=Peptide[i:i+k]
                # When the substring exceeds the end of the string, continue from the start of the string
            else:
                substring=Peptide[i:]+Peptide[:i+k-n]
            substrings.append(substring)

    return substrings

'''
Calculates the mass of a given amino acid string according to the amino acid dictionary
'''
def CalculateMass(AA: str) -> int:
    mass=0

    for amino in AA:
        mass+=AAmass[amino]
    
    return mass

'''LinearSubPep takes a string of amino acid Peptide as input
and outputs a list of linearsubpeptides of Peptide, without the wrap-arounds.
'''
def LinearSubPep(Peptide: str) -> list[str]:
    n=len(Peptide)
    substrings=[]
    for k in range(1,n):
        for i in range(n-k+1):
            substring=Peptide[i:i+k]
            substrings.append(substring)
    return substrings

'''
Scores a peptide against a spectrum according to the matches between the theoretical cyclo spectrum of the peptide and the spectrum
'''
def ScorePeptide(Peptide: str, Spectrum: list[int]) -> int:
    theospec=CycloSpectrum(Peptide)
    n=len(theospec)
    k=len(Spectrum)
    i=0
    j=0
    count=0
    while i<n and j<k :
        if theospec[i]==Spectrum[j]:
            count+=1
            i+=1
            j+=1
        elif theospec[i]<Spectrum[j]:
            i+=1
        else:
            j+=1
    return count

'''
Scores a peptide against a spectrum according to the matches between the theoretical linear spectrum of the peptide and the spectrum
'''
def LinearScorePeptide(Peptide: str, Spectrum: list[int]) -> int:
    theospec=LinearSpectrum(Peptide)
    n=len(theospec)
    k=len(Spectrum)
    i=0
    j=0
    count=0
    while i<n and j<k :
        if theospec[i]==Spectrum[j]:
            count+=1
            i+=1
            j+=1
        elif theospec[i]<Spectrum[j]:
            i+=1
        else:
            j+=1
    return count

'''
CyclopeptideSequencing takes a spectrum (a list of interger masses) as input
and returns a list of mass-strings that could generate this spectrum
'''
def CyclopeptideSequencing(Spectrum : list[int]):
    Wholemass=Spectrum[len(Spectrum)-1]
    # start with single amino acids
    CandidatePeptides=list(AAmass.keys())
    FinalPeptides=[]

    # Keep branching and bounding the candidate peptides until all possibilities are explored
    while len(CandidatePeptides)>0:
        NewCandidates=Copy(CandidatePeptides)

        # range through all the current possible candidate peptides 
        for peptide in CandidatePeptides:
            subspectrum=LinearSpectrum(peptide)
            # if the linear spectrum of the candidate peptide is not consistent with the Spectrum
            # cancel its candidacy
            if Consistent(subspectrum,Spectrum)==False:
                NewCandidates.remove(peptide)
            # if the mass of the candidate is equal or greater than the whole mass indictaed by the Spectrum
            # cancel its candidacy
            elif subspectrum[len(subspectrum)-1]>=Wholemass:
                NewCandidates.remove(peptide)
                # If the spectrum is consistent with the cyclospectrum of the peptide
                # this peptide is a match for the Specturm
                if Consistent(Spectrum,CycloSpectrum(peptide)):
                    FinalPeptides.append(peptide)
        CandidatePeptides=Expand(NewCandidates)

    # Convert the peptides into mass strings
    FinalMasses=[]
    for peptide in FinalPeptides:
        massstr=PeptideToMassString(peptide)
        if massstr not in FinalMasses:
            FinalMasses.append(massstr)

    return FinalMasses   

# Converts a peptide to a mass-string
def PeptideToMassString(Peptide):

    mass_string=str(AAmass[Peptide[0]])
    for i in range(len(Peptide)-1):
        mass_string+='-'+str(AAmass[Peptide[i+1]])

    return mass_string

# Checks if the spectrum of the peptide is consistent with Spec
# A peptide is consistent with Spectrum if every mass in its theoretical spectrum is contained in Spectrum
def Consistent(SubSpec, Spec):
    for mass in SubSpec:
       if mass not in Spec:
           return False

    return True

# Expand takes a list of peptides and 
# expands each peptide with every amino acid possible and returns the expanded peptides
def Expand(peptides):
    newpeptides=[]
    for peptide in peptides:
        for amino in AAmass.keys():
            newpeptides.append(peptide+amino)        
   
    return newpeptides

# Deep copies a list from another
def Copy(listB: list) -> list:
    listA=[]

    for i in range(len(listB))  :
        listA.append(listB[i])

    return listA

'''
Implements an alternative of the CyclopeptideSequencing by keeping a leaderboard, addressing the errors in a real spectrum
'''
def LeaderboardCyclopeptideSequencing(Spectrum, N):
    Wholemass=Spectrum[len(Spectrum)-1]
    # start with single amino acids
    Leaderboard=list(AAmass.keys())
    LeaderPeptide=''

    # keep trimming and expanding the leaderboard until it's empty
    while len(Leaderboard)>0:
        NewLeaderboard=Copy(Leaderboard)

        # range through every peptide in the current leaderboard
        for peptide in Leaderboard:
            mass=CalculateMass(peptide)
            if mass==Wholemass:
                # if the score of the peptide is higher than the leaderpeptide
                # update the leaderpeptdie
                if ScorePeptide(peptide,Spectrum)>ScorePeptide(LeaderPeptide,Spectrum):
                    LeaderPeptide=peptide
            # if the mass of the peptide exceeds the whole mass 
            # cancel its candidacy
            elif mass>Wholemass:
                NewLeaderboard.remove(peptide)
        # trim and expand the current leaderboard to get the new leaderboard for the next iteration
        Leaderboard=Expand(Trim(NewLeaderboard,Spectrum,N))

    return PeptideToMassString(LeaderPeptide)

'''
Trim takes a list of pepetides Peptides, a list of integer masses Spectrum and an interger N as input
and returns the peptides with top N (with ties) scores against the spectrum
'''
def Trim(Peptides,Spectrum,N):
    # construct a score dictionary that stores every score and every peptide corresponding to that score 
    Scoredict={}
    for pep in Peptides:
        score=LinearScorePeptide(pep,Spectrum)
        if score not in Scoredict.keys():
            Scoredict[score]=[pep]
        else:
            Scoredict[score].append(pep)
    Scores=list(Scoredict.keys())

    # Sort the scores in descending order
    Scores.sort(reverse=True)
    board=[]

    # append the peptides in descending order until the number exceeds N
    num=0
    for i in range(len(Scores)):
        board.extend(Scoredict[Scores[i]])
        num+=len(Scoredict[Scores[i]])
        if num>=N:
            break
    
    return board

# Returns the convolution of the given spectrum.
# The convolution of a spectrum is  all positive differences of masses in the spectrum.
def SpectralConvolution(Spectrum: list[int]):
    convolution=[]
    # Assuming the spectrum is in ascending order
    # deduct every mass by all the masses greater than it
    for i in range(0,len(Spectrum)-1):
        for j in range(i+1,len(Spectrum)):
            con=Spectrum[j]-Spectrum[i]
            # do not include 0s
            if con !=0:
                convolution.append(con)

    return convolution

'''
Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.
Output: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements (and ties) of the convolution of Spectrum 
that fall between 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).
'''
def ConvolutionCyclopeptideSequencing(M,N,Spectrum):
    # First construct the convolution dictionary of the spectrum, 
    # storing every mass in convolution along with the number of appearance of that mass
    convolution=SpectralConvolution(Spectrum)
    condict={}
    for mass in convolution:
        if mass >=57 and mass<=200:
            if  mass not in condict.keys():
                condict[mass]=0
            condict[mass]+=1
    
    # Obtain a keylist of the convolution dictionary where the keys are sorted by descending number of appearances
    keylist=list(condict.keys())
    keylist.sort(key=lambda x:condict[x],reverse=True)
    
    # take the elements in the convolution with top M (wth ties) appearances
    value=condict[keylist[0]]
    elements=[]
    for i in range(len(keylist)):
        if condict[keylist[i]]==value:
            elements.append(keylist[i])
        elif i<=M:
            elements.append(keylist[i])
            value=condict[keylist[i]]
        else:
            break

    # run the mass-tide version of the LeaderboardCyclopeptideSequencing algorithm
    Wholemass=Spectrum[len(Spectrum)-1]
    # start with the single elements in elements
    Leaderboard=[]
    for ele in elements:
        Leaderboard.append([ele])

    LeaderPeptide=[]
    while len(Leaderboard)>0:
        NewLeaderboard=Copy(Leaderboard)
        for masstide in Leaderboard:
            mass=Sum(masstide)
            if mass==Wholemass:
                if ScoreMasstide(masstide,Spectrum)>ScoreMasstide(LeaderPeptide,Spectrum):
                    LeaderPeptide=masstide
                NewLeaderboard.remove(masstide)
            elif mass>Wholemass:
                NewLeaderboard.remove(masstide)
        Leaderboard=ElementExpand(TrimMasstides(NewLeaderboard,Spectrum,N),elements)

    return LeaderPeptide

'''
Scores a masstide against a spectrum according to the matches between the theoretical cyclo spectrum of the masstide and the spectrum
'''
def ScoreMasstide(Masstide: list[int], Spectrum: list[int]) -> int:
    theospec=MassSpec(CycloSubMass(Masstide))
    theospec.sort()
    n=len(theospec)
    k=len(Spectrum)
    i=0
    j=0
    count=0
    while i<n and j<k :
        if theospec[i]==Spectrum[j]:
            count+=1
            i+=1
            j+=1
        elif theospec[i]<Spectrum[j]:
            i+=1
        else:
            j+=1
    return count

'''
Scores a masstide against a spectrum according to the matches between the theoretical linear spectrum of the masstide and the spectrum
'''
def LinearScoreMasstide(Masstide: str, Spectrum: list[int]) -> int:
    linearsubmasses=LinearSubMass(Masstide)
    linearsubmasses.append(Masstide)
    theospec=MassSpec(linearsubmasses)
    theospec.insert(0,0)
    theospec.sort()
    n=len(theospec)
    k=len(Spectrum)
    i=0
    j=0
    count=0
    while i<n and j<k :
        if theospec[i]==Spectrum[j]:
            count+=1
            i+=1
            j+=1
        elif theospec[i]<Spectrum[j]:
            i+=1
        else:
            j+=1
    return count

'''LinearSubMass takes a masstide as input
and outputs a list of linearsubmasstides of Masstide, with the wrap-arounds.
'''
def LinearSubMass(Masstide:list[int]) -> list[int]:
    n=len(Masstide)
    substrings=[]
    for k in range(1,n):
        for i in range(n-k+1):
            substring=Masstide[i:i+k]
            substrings.append(substring)
    return substrings

'''
TrimMasstides takes a list of masstides Masstides, a list of integer masses Spectrum and an interger N as input
and returns the masstides with top N (with ties) scores against the spectrum
'''
def TrimMasstides(Masstides,Spectrum,N):
    # construct a score dictionary that stores every score and every peptide corresponding to that score 
    Scoredict={}
    for mass in Masstides:
        score=LinearScoreMasstide(mass,Spectrum)
        if score not in Scoredict.keys():
            Scoredict[score]=[mass]
        else:
            Scoredict[score].append(mass)

    # Sort the scores in descending order
    Scores=list(Scoredict.keys())
    Scores.sort(reverse=True)

    # append the masstides in descending order until the number exceeds N
    board=[]
    num=0
    for i in range(len(Scores)):
        board.extend(Scoredict[Scores[i]])
        num+=len(Scoredict[Scores[i]])
        if num>=N:
            break
    
    return board

# given a list of submasstides, return the mass spectrum
def MassSpec(SubMasses):
    spec=[]
    for submass in SubMasses:
        spec.append(Sum(submass))
    return spec

'''CycloSubMass takes a masstide as input
and outputs a list of cyclo submasstides of Masstide, without the wrap-arounds.
'''
def CycloSubMass(Masstide:list[int]) -> list[str]:
    n=len(Masstide)
    subtides=[[0]]
    for k in range (1, n, 1):
        for i in range (n):
            if i+k-1<n:
                subtide=Masstide[i:i+k]
            else:
                subtide=Masstide[i:]
                subtide.extend(Masstide[:i+k-n])
            subtides.append(subtide)
    subtides.append(Masstide)
    return subtides

# Calculates the mass of a given masstide
def Sum(Masstide:list[int]):
    sum=0
    for mass in Masstide:
        sum+=mass
    return sum

'''
Given a list of list of integers Board,
expand every list of intergers in the Board with every element in the list of Elements
and return the expanded board
'''
def ElementExpand(Board,Elements):
    NewBoard=[]
    
    for pep in Board:
        for ele in Elements:
            newpep=Copy(pep)
            newpep.append(ele)
            NewBoard.append(newpep)

    return NewBoard
