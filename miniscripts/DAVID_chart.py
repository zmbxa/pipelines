### 
## get input from cmd
import argparse

parser = argparse.ArgumentParser(description="This script is for DAVID gene enrichment. It take genelist as input and calling for DAVID API, give table-format result.")
parser.add_argument('-e','--email', type=str,default='niuyuxiao@westlake.edu.cn',dest="mail_addr")
parser.add_argument('-idt','--IDType', type=str, help='ID type, one of ENSEMBL_GENE_ID, ENTREZ_GENE_ID, use ENSEMBL if not specified',default='ENSEMBL_GENE_ID',dest="idType")
parser.add_argument('-g', '--geneID', type=str, help='comma split gene ID',dest="gene_list")
parser.add_argument('-o', '--output', type=str, help='output file name, ./DAVIDchartReport.txt if not specified',default='DAVIDchartReport.txt',dest="out")
parser.add_argument('-thd','--threshold', type=float,help='FDR threshold',default=0.05,dest="thd")
parser.add_argument('-ct','--ct', type=int,help='gene Count number threshold, default is 2',default=2,dest="ct")
args = parser.parse_args()

if not args.gene_list:
  print(" Gene list is required!")
  parser.print_usage()
  exit(1)

import datetime
print("\t--------\t",datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),"\t--------\t")
print("DAVID enrichment for ", args.gene_list," with ID type as",args.idType)


## DAVID enrichment table
import sys
sys.path.append('../')

import logging
import traceback
from suds import *
from suds.client import Client
from datetime import datetime

errors = 0

logging.basicConfig(level=logging.DEBUG)

url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'

print('url=%s' % url)

#
# create a service client using the wsdl.
#
client = Client(url)
client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
#
# print the service (introspection)
#
#print(client)

#authenticate user email 
client.service.authenticate(args.mail_addr)

#add a list 
inputIds = args.gene_list
idType = args.idType
listName = 'make_up'
listType = 0
print(client.service.addList(inputIds, idType, listName, listType))

#print(client.service.getDefaultCategoryNames())
categorySting = str(client.service.setCategories('BBID,BIOCARTA,COG_ONTOLOGY,GOTERM_BP_1,GOTERM_BP_2,GOTERM_BP_3,GOTERM_BP_4,GOTERM_BP_5,GOTERM_BP_ALL,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'))
categories=categorySting.split(',')

#getChartReport
chartReport = client.service.getChartReport(args.thd,args.ct)
chartRow = len(chartReport)
print('Total chart records:',chartRow)

resF = args.out
with open(resF, 'w') as fOut:
    fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
    for simpleChartRecord in chartReport:
        categoryName = simpleChartRecord.categoryName
        termName = simpleChartRecord.termName
        listHits = simpleChartRecord.listHits
        percent = simpleChartRecord.percent
        ease = simpleChartRecord.ease
        genes = simpleChartRecord.geneIds
        listTotals = simpleChartRecord.listTotals
        popHits = simpleChartRecord.popHits
        popTotals = simpleChartRecord.popTotals
        foldEnrichment = simpleChartRecord.foldEnrichment
        bonferroni = simpleChartRecord.bonferroni
        benjamini = simpleChartRecord.benjamini
        FDR = simpleChartRecord.afdr
        rowList = [categoryName, termName, str(listHits), str(percent), str(ease), genes, str(listTotals), str(popHits), str(popTotals), str(foldEnrichment), str(bonferroni), str(benjamini), str(FDR)]
        fOut.write('\t'.join(rowList) + '\n')
    print('write file:', resF, 'finished!')

