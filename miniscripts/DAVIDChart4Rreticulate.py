import argparse
import datetime
import logging
import sys
from suds.client import Client
import pandas as pd


def david_gene_enrichment(email='niuyuxiao@westlake.edu.cn', IDType='ENSEMBL_GENE_ID', gene_list=None, output='DAVIDchartReport.txt', threshold=0.05, ct=2):
    logging.basicConfig(level=logging.DEBUG)

    url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'

    client = Client(url)
    client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')

    # Authenticate user email
    client.service.authenticate(email)

    if gene_list is None:
        print(" Gene list is required!")
        return

    print("\t--------\t", datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "\t--------\t")
    print("DAVID enrichment for ", gene_list, " with ID type as", IDType)

    # Add a list
    inputIds = gene_list
    idType = IDType
    listName = 'make_up'
    listType = 0
    print(client.service.addList(inputIds, idType, listName, listType))

    categorySting = str(client.service.setCategories('BBID,BIOCARTA,COG_ONTOLOGY,GOTERM_BP_1,GOTERM_BP_2,GOTERM_BP_3,GOTERM_BP_4,GOTERM_BP_5,GOTERM_BP_ALL,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'))
    categories = categorySting.split(',')

    # Get Chart Report
    chartReport = client.service.getChartReport(threshold, ct)
    chartRow = len(chartReport)
    print('Total chart records:', chartRow)
    
    Chart_df = pd.DataFrame()
    for record in chartReport:
      print(record.termName)
      Chart_df=pd.concat([Chart_df,pd.DataFrame(record).T[1:]])
    
    Chart_df.columns = list(pd.DataFrame(record)[0])
    
    return Chart_df

    # with open(output, 'w') as fOut:
    #     fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
    #     for simpleChartRecord in chartReport:
    #         categoryName = simpleChartRecord.categoryName
    #         termName = simpleChartRecord.termName
    #         listHits = simpleChartRecord.listHits
    #         percent = simpleChartRecord.percent
    #         ease = simpleChartRecord.ease
    #         genes = simpleChartRecord.geneIds
    #         listTotals = simpleChartRecord.listTotals
    #         popHits = simpleChartRecord.popHits
    #         popTotals = simpleChartRecord.popTotals
    #         foldEnrichment = simpleChartRecord.foldEnrichment
    #         bonferroni = simpleChartRecord.bonferroni
    #         benjamini = simpleChartRecord.benjamini
    #         FDR = simpleChartRecord.afdr
    #         rowList = [categoryName, termName, str(listHits), str(percent), str(ease), genes, str(listTotals), str(popHits), str(popTotals), str(foldEnrichment), str(bonferroni), str(benjamini), str(FDR)]
    #         fOut.write('\t'.join(rowList) + '\n')
    #     print('Write file:', output, 'finished!')
