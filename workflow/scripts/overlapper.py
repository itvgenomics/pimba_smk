#!/usr/bin/env python

###    OverlapPER may be used to merge overlapping paired-end reads
###    Copyright (C) 2016  Renato Oliveira
###
###    This program is free software: you can redistribute it and/or modify
###    it under the terms of the GNU General Public License as published by
###    the Free Software Foundation, either version 3 of the License, or
###    any later version.
###
###    This program is distributed in the hope that it will be useful,
###    but WITHOUT ANY WARRANTY; without even the implied warranty of
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###    GNU General Public License for more details.
###
###    You should have received a copy of the GNU General Public License
###    along with this program.  If not, see <http://www.gnu.org/licenses/>.

###Contacts:
###    guilherme.oliveira@itv.org
###    renato.renison@gmail.com

###Usage:
###		python overlapper.py -f <forward_reads.fastq> -r <reverse_reads.fastq> --mo <min_overlap> --ms <min_similarity>


import sys
import getopt

def usage():
	text="\nOverlapPER: Overlapping Paired-End Reads\nUsage: python overlapper.py [options]\n"
	text+="Options:\n"
	text+="\t-h\n\t\tShow this message.\n"
	text+="\t-f\n\t\tPath to forward reads\n"
	text+="\t-r\n\t\tPath to reverse reads\n"
	text+="\t--mo\n\t\tLength of the minimum overlap between the paired\n\t\treads (Default 25).\n"
	text+="\t--ms\n\t\tPercentage of the accepted minimum similarity in an overlap\n\t\tregion oftwo paired reads (default is 0.9).\n"

	print(text)



#creates a list containing the reads sequences and another list 
#with the phred quality information
def createReadsList(arqReads):
    lines = arqReads.readlines()
    readList = []
    phredList = []
    idReadList = []
    i = 0

    while i + 3 < len(lines):
        header = lines[i].strip()
        sequence = lines[i + 1].strip()
        plus = lines[i + 2].strip()
        quality = lines[i + 3].strip()

        # Check if this is a valid FASTQ block
        if header.startswith('@') and plus == '+' and len(sequence) == len(quality):
            sequence = sequence.replace("R", "N").replace("Y", "N").replace("S", "N") \
                               .replace("W", "N").replace("K", "N").replace("M", "N") \
                               .replace("B", "N").replace("D", "N").replace("H", "N") \
                               .replace("V", "N")
            readList.append(sequence)
            phredList.append(quality)
            idReadList.append(header)
            i += 4
        else:
            # Skip malformed entry
            print(f"Warning: Skipped malformed entry at line {i+1}")
            i += 1  # Move forward by one and resync

    return readList, phredList, idReadList

#creates the dictionary able to translate the quality code to numerical code
def createDicPhred():
	dic={'%':4,'*':9,'&':5,')':8,"'":6,',':11,'/':14,'7':22,'M':44,'6':21,'R':49,'9':24,'?':30,'3':18,'0':15,'_':62,'O':46,'A':32,'H':39,'F':37,'P':47,'W':54,'L':43,'4':19,'\\':59,'Y':56,'T':51,'[':58,'U':52,'G':38,'V':53,'S':50,'J':41,'E':36,'>':29,'I':40,'@':31,'=':28,'D':35,';':26,'B':33,'K':42,'C':34,'N':45,'5':20,'8':23,'<':27,'-':12,'+':10,'1':16,':':25,'2':17,'(':7,'Q':48,'.':13,'!':0,'$':3,'"':1,'#':2,'X':55,'Z':57,'^':61,']':60}

	return dic
def createInverseDicPhred():
	dic={4:'%',9:'*',5:'&',8:')',6:"'",11:',',14:'/',22:'7',44:'M',21:'6',49:'R',24:'9',30:'?',18:'3',15:'0',62:'_',46:'O',32:'A',39:'H',37:'F',47:'P',54:'W',43:'L',19:'4',59:'\\',56:'Y',51:'T',58:'[',52:'U',38:'G',53:'V',50:'S',41:'J',36:'E',29:'>',40:'I',31:'@',28:'=',35:'D',26:';',33:'B',42:'K',34:'C',45:'N',20:'5',23:'8',27:'<',12:'-',10:'+',16:'1',25:':',17:'2',7:'(',48:'Q',13:'.',0:'!',3:'$',1:'"',2:'#',55:'X',57:'Z',61:'^',60:']'}
	
	return dic
dic=createDicPhred()
dicInverse = createInverseDicPhred()

#Returns the reverse complement of the read, as well as the respective quality information	
def compReverso(readList,phredList):
	dic={'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g','g':'c','t':'a','n':'n'}
	compReversoList=[]
	qualReversoList=[]
	j=0

	while(j<len(readList)):
		read = readList[j]
		phred = phredList[j]
		i=len(read)-1
		cReverso=""
		phredReverso=""
		while(i>=0):
			cReverso+=dic[read[i]]
			phredReverso+=phred[i]
			i-=1
		compReversoList.append(cReverso)
		qualReversoList.append(phredReverso)
		j+=1

	return (compReversoList,qualReversoList)

#Returns the reverse sequence
def reverse(x):
	reverseX = ""
	i=len(x)-1
	while(i>=0):
		reverseX+=x[i]
		i-=1
	
	return reverseX


def indel(x,y,phredX,phredY,j):
	tol=5
	hit=0
	rest=0
	repeats=0
	result = [0,x,y,phredX,phredY,-1]
	novoX = x[0:j]+"N"+x[j:len(x)]
	novoY = y[0:j]+"N"+y[j:len(y)]
	novoPhredX = phredX[0:j]+"!"+phredX[j:len(phredX)]
	novoPhredY = phredY[0:j]+"!"+phredY[j:len(phredY)]
	while(repeats<5):
		menor = len(y)
		if(len(novoX)<menor): menor=len(novoX)
		if(menor-j <= tol):
			rest = menor
		else:
			rest = j+tol+1 #avoids index out of range error
		for i in range(j+1,rest):
			if(novoX[i]==y[i]):
				hit+=1
		if(hit>=rest-j-1 and hit!=0):
			result[0]=rest-j-1
			result[1]=novoX
			result[2]=y
			result[3]=novoPhredX
			result[4]=phredY
			result[5]=hit
			return result
		
		hit=0
		menor = len(x)
		if(len(novoY)<menor): menor=len(novoY)
		if(menor-j <= tol):
			rest = menor
		else:
			rest = j+tol+1 #avoids index out of range error
		for i in range(j+1,rest):
			if(novoY[i]==x[i]):
				hit+=1
		if(hit>=rest-j-1 and hit!=0):
			result[0]=rest-j-1
			result[1]=x
			result[2]=novoY
			result[3]=phredX
			result[4]=novoPhredY
			result[5]=hit
			return result
		
		hit=0
		repeats+=1
		novoY = novoY[0:j+1]+"N"+novoY[j+1:len(novoY)]
		novoX = novoX[0:j+1]+"N"+novoX[j+1:len(novoX)]
		novoPhredX = novoPhredX[0:j+1]+"!"+novoPhredX[j+1:len(novoPhredX)]
		novoPhredY = novoPhredY[0:j+1]+"!"+novoPhredY[j+1:len(novoPhredY)]
		j+=1
	return result

def mismatch(x,y,j):
	tol=5
	hit=0
	rest=0
	result = [0,x,y,-1]
	novoX = x
	novoY = y

	menor = len(novoY)
	if(len(novoX)<menor): menor=len(novoX)
	if(menor-j<=tol): 
		rest = menor
	else:
		rest = j+tol+1 #avoids index out of range error
	for i in range(j+1,rest): 
		if(novoX[i]==novoY[i]):
			hit+=1
	if(hit>=rest-j-1 and hit!=0):
		result[0]=rest-j-1
		result[1]=novoX
		result[2]=novoY
		result[3]=hit
		return result
		
		
	return result

def decidePhred(phredX,phredY):
	if(dic[phredX]>dic[phredY]):
		return (phredX,0) #se segundo argumento do retorno for 0, entao o escolhido foi o x
	else:
		return (phredY,1) #se for 1, o escolhido foi o y

def sumPhred(phredX, phredY):
	if(dic[phredX]+dic[phredY] >= 62):
		return dicInverse[62]
	else:
		return dicInverse[dic[phredX]+dic[phredY]]

def extendSeed(extSeedX, extSeedY, extSeedPhredX, extSeedPhredY, acerto, erro): #Antigo checkIdentity
	
	length = len(extSeedX)
	if(len(extSeedY)>length): length = len(extSeedY)
	
	consenso = ""
	phredConsenso=""
	i=0
	endOverlap=0
	while(i<=length):
		#print(extSeedX)
		#print(extSeedY)
		#print(consenso)
		#print("\n")
		
		if(i>=len(extSeedY)):
			consenso+=extSeedX[i:len(extSeedX)]
			phredConsenso+=extSeedPhredX[i:len(extSeedPhredX)]
			endOverlap=i
			break
		elif(i>=len(extSeedX)):
			consenso+=extSeedY[i:len(extSeedY)]
			phredConsenso+=extSeedPhredY[i:len(extSeedPhredY)]
			endOverlap=i
			break
		elif(extSeedX[i]==extSeedY[i]):
			consenso+= extSeedX[i]
			phred = sumPhred(extSeedPhredX[i],extSeedPhredY[i])
			phredConsenso+=phred
			acerto+=1
		elif(extSeedX[i]=="N"):
			consenso+= extSeedY[i]
			phredConsenso+=extSeedPhredY[i]
			acerto+=1
		elif(extSeedY[i]=="N"):
			consenso+= extSeedX[i]
			phredConsenso+=extSeedPhredX[i]
			acerto+=1
		else:
			
			
			resMismatch = mismatch(extSeedX,extSeedY,i)
			if(resMismatch[3]>=resMismatch[0]):

				extSeedX = resMismatch[1]
				extSeedY = resMismatch[2]
				erro+=1
				phred,escolhido=decidePhred(extSeedPhredX[i],extSeedPhredY[i])
				if(escolhido):	consenso+=extSeedY[i]
				else:	consenso+=extSeedX[i]
				phredConsenso+=phred

			else:

				resIndel = indel(extSeedX,extSeedY,extSeedPhredX,extSeedPhredY,i)
				if(resIndel[5]>=resIndel[0]):

					extSeedX = resIndel[1]
					extSeedY = resIndel[2]
					extSeedPhredX = resIndel[3]
					extSeedPhredY = resIndel[4]
					i-=1

				else:
					phred,escolhido=decidePhred(extSeedPhredX[i],extSeedPhredY[i])
					if(escolhido):	consenso+=extSeedY[i]
					else:	consenso+=extSeedX[i]
					phredConsenso+=phred
					erro+=1

		length=len(extSeedX)
		i+=1
		if(len(extSeedY)>length): 
			length = len(extSeedY)
		

	#print(acerto, erro)
	return (consenso, phredConsenso, acerto, erro, endOverlap)

def assembleReads(x,y, phredX, phredY, minOverlap,minSimilarity):
	contig=""
	
	phredContig=""
	seedLength=16
	seedStep=12
	yLen = len(y)

	for i in range(yLen,int(yLen/2),-seedStep):
		if(i<seedLength): break
		seed = y[i-seedLength:i]
		iniX=x.rfind(seed)

		if(iniX!=-1):

			acerto=0
			erro=0
			dirSeedX=x[iniX:len(x)]
			dirSeedPhredX=phredX[iniX:len(phredX)]
			dirSeedY=y[i-seedLength:len(y)]
			dirSeedPhredY=phredY[i-seedLength:len(phredY)]

			esqSeedX=x[0:iniX]
			esqSeedPhredX = phredX[0:iniX]
			esqSeedY=y[0:i-seedLength]
			esqSeedPhredY = phredY[0:i-seedLength]
			
			esqSeedX = reverse(esqSeedX)
			esqSeedPhredX = reverse(esqSeedPhredX)
			esqSeedY=reverse(esqSeedY)
			esqSeedPhredY = reverse(esqSeedPhredY)

			dirSeedConsenso,dirPhredSeedConsenso,acerto,erro,dirEndOverlap = extendSeed(dirSeedX,dirSeedY,dirSeedPhredX,dirSeedPhredY,acerto,erro)
			if(dirSeedConsenso==""): break	
			esqSeedConsenso,esqPhredSeedConsenso,acerto,erro,esqEndOverlap = extendSeed(esqSeedX,esqSeedY,esqSeedPhredX,esqSeedPhredY,acerto,erro)
			if(esqSeedConsenso==""): break
			
			overlapLength = dirEndOverlap+esqEndOverlap
			ident = float(acerto)/overlapLength
			#print(overlapLength, ident)

			if(overlapLength>=minOverlap and ident>=minSimilarity):

				esqSeedConsenso = reverse(esqSeedConsenso)
				esqPhredSeedConsenso = reverse(esqPhredSeedConsenso)
				contig+=esqSeedConsenso+dirSeedConsenso
				phredContig+=esqPhredSeedConsenso+dirPhredSeedConsenso

				return (contig,phredContig)

					
	return (contig,phredContig)


try:
	opts, args = getopt.getopt(sys.argv[1:], "f:r:h", ["mo=","ms="])
except getopt.GetoptError as err:
    # print help information and exit:
    print(str(err))  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

minOverlap = 25 #Min acceptable overlap
minSimilarity = 0.9 #Min acceptable similarity in the Overlap

for o, a in opts:
	if o in ("-h", "--help"):
		usage()
		sys.exit()
	elif o in ("-f"):
		nameForward = a #Forward reads
	elif o in ("-r"):
		nameReverse = a #Reverse reads
	elif o in ("--mo"):
		minOverlap = a
	elif o in ("--ms"):
		minSimilarity = a
	else:
		assert False, "unhandled option"


minOverlap = int(minOverlap)
minSimilarity = float(minSimilarity)

arqForward = open(nameForward,"r")
arqReverse = open(nameReverse,"r")

forwardList,phredForward,idReadsForward = createReadsList(arqForward)
reverseList,phredReverse,idReadsReverse = createReadsList(arqReverse)
compReversoList,phredReverso = compReverso(reverseList,phredReverse)

arqContigs = open("overlapped.fastq","w")
textContig = ""
arqNotAssembled1 = open("notAssembled-1.fastq","w")
textNotAssembled1 = ""
arqNotAssembled2 = open("notAssembled-2.fastq","w")
textNotAssembled2 = ""

assembled=0
notAssembled=0
for i in range(0,len(forwardList)):

	contig,phredContig = assembleReads(forwardList[i], compReversoList[i], phredForward[i], phredReverso[i], minOverlap, minSimilarity)
	if(len(contig)>0):
		textContig+=idReadsForward[i]+"\n"
		textContig+=contig+"\n"
		textContig+="+\n"
		textContig+=phredContig+"\n"
		assembled+=1

	else:
		contig,phredContig = assembleReads(compReversoList[i], forwardList[i], phredReverso[i], phredForward[i], minOverlap, minSimilarity)
		if(len(contig)>0):
			textContig+=idReadsForward[i]+"\n"
			textContig+=contig+"\n"
			textContig+="+\n"
			textContig+=phredContig+"\n"
			assembled+=1
		else: 
			notAssembled+=1
			textNotAssembled1+=idReadsForward[i]+"\n"
			textNotAssembled1+=forwardList[i]+"\n"
			textNotAssembled1+="+\n"
			textNotAssembled1+=phredForward[i]+"\n"
			
			textNotAssembled2+=idReadsReverse[i]+"\n"
			textNotAssembled2+=reverseList[i]+"\n"
			textNotAssembled2+="+\n"
			textNotAssembled2+=phredReverse[i]+"\n"
	


print("Merged reads: "+str(assembled))
print("\nNot merged reads: "+str(notAssembled))
arqContigs.write(textContig)
arqNotAssembled1.write(textNotAssembled1)
arqNotAssembled2.write(textNotAssembled2)
arqContigs.close()
arqNotAssembled1.close()
arqNotAssembled2.close()
arqForward.close()
arqReverse.close()



