'''
@author: XiongXiong

@date: 2022/11/23

'''

import sys

if __name__ == '__main__':
	
	countMtx = sys.argv[1]
	
	with open(countMtx, 'r') as Mtx, open('tmp.Matrix.txt', 'w') as MtxOut, open('tmp.Genes.txt', 'w') as GeneOut, open('tmp.Barcode.txt', 'w') as BarcodeOut:
		
		geneCount = 0
		cbCount = 0
		geneDict = dict()
		cbDict = dict()

		for line in Mtx:
			lineSet = line.strip('\n').split('\t')
			geneID = lineSet[0]+'\t'+lineSet[1]
			barcode = lineSet[2]
			nCount = lineSet[3]
			if not(geneID in geneDict):
				geneDict[geneID] = geneCount + 1
				geneCount += 1

			if not(barcode in cbDict):
				cbDict[barcode] = cbCount + 1
				cbCount += 1

			MtxOut.write(str(geneDict.get(geneID))+' '+str(cbDict.get(barcode))+' '+str(nCount)+'\n')
			GeneOut.write(str(geneDict.get(geneID))+'\t'+geneID+'\n')
			BarcodeOut.write(str(cbDict.get(barcode))+'\t'+barcode+'\n')

