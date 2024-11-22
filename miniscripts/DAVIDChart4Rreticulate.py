import argparse
import datetime
import logging
import sys
from suds.client import Client
import pandas as pd


def david_gene_enrichment(email='niuyuxiao@westlake.edu.cn', IDType='ENSEMBL_GENE_ID', gene_list=None, output='DAVIDchartReport.txt', threshold=0.05, ct=2):
    logging.basicConfig(level=logging.DEBUG)

    #url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
    url = 'https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl'

    client = Client(url,timeout=600)
    #client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
    client.wsdl.services[0].setlocation('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')

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
    
    chartDict = []

    for record in chartReport:
        record_dict = dict(record)
        del record_dict['scores']
        chartDict.append(record_dict)
       
    chart_df = pd.DataFrame(chartDict)
    return chart_df

