'''
@author: XiongXiong

@date: 2022/11/16

'''

import Levenshtein
import sys


def BC12(sequence, bait, contaTail, starts, bcN, lens):
	distance = Levenshtein.distance(sequence[(starts+bcN):(starts+bcN+15)], bait)
	if distance > 3:
		return 'None'
	else:
		disList = [15, 15, 15]
		for shifts in range(-1, 2):
			disList[shifts+1] = Levenshtein.distance(sequence[(starts+bcN+shifts):(starts+bcN+shifts+15)], bait)			
		
	posi = 0
	if min(disList) < distance:
		posi = disList.index(min(disList))-1
	
	if sequence[starts+bcN+posi+14] == contaTail:
		return 'Contamination'
	else:
		bc = sequence[(starts+bcN+posi+15):(starts+bcN+posi+15+lens)]
		return(bc)

if __name__ == '__main__':
    
	fastq = sys.argv[1]
	libType = sys.argv[2]
		
	if libType == 'DNA':
		bc1Bait = 'GGCCAGAGCATTCGA'
		contaTail = 'T'
	elif libType == 'RNA':
		bc1Bait = 'GGCCAGAGCATTCGT'
		contaTail = 'A'
	
	with open(fastq, 'r') as fq, open('tmp.ReadType.txt', 'w') as out:
	
		for line in fq:
			names = line.strip('\n').split("\t")[0]
			fastq2 = line.strip('\n').split("\t")[1]

			if len(fastq2) < 90:
				typ = 'TooShort'
				barcode = 'None'
			elif 'GGGGGGGGGGGGGGGGGGGG' in fastq2:
				typ = 'polyG'
				barcode = 'None'
			else:
				disList = [0, 0, 0, 0]
				for numN in range(0, 4):
					disList[numN] = Levenshtein.hamming(fastq2[(17+numN):(17+numN+15)], 'GTGGCCGATGTTTCG')

				bcN = -1
				if min(disList) <= 3:
					bcN = disList.index(min(disList))
				
				bc2 = BC12(fastq2, 'GTGCGAACTCAGACC', contaTail, 33, bcN, 7)
				bc1 = BC12(fastq2, bc1Bait, contaTail, 71, bcN, 4)
								
				if bcN == -1 or bc2 == 'None' or bc1 == 'None':
					typ = 'random'
					barcode = 'None'
				elif bc1 == 'Contamination':
					typ = 'Contamination'
					barcode = 'None'					
				else:
					typ = 'FullyLigated'
					barcode = fastq2[10:17]+bc2+bc1

			out.write(names+'\t'+typ+'\t'+barcode+'\n')

