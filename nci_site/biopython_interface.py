################################################
#
# Description: Interface between the rest of the site and Biopython/Entrez
#
# Author: George Dietz
#
# Make sure to rerun relearn() if you get an ImportError vis a visa self.committee and NCI_learn
#
###############################################

from Bio import Entrez, Medline
from NCI import NCI_learn
from NCI import NCI_predict


Entrez.email = "george_dietz@brown.edu"
DEFAULT_SEARCH_STRING = '((((((cancer OR neoplasm)) OR cancer[MeSH Terms]) OR neoplasm[MeSH Terms])) AND (lung OR colon OR rectal OR colorectal OR bladder OR breast OR melanoma)) AND ((predictive OR prognostic OR prediction OR prognosis OR prognostication)) AND ((models OR models OR modeling OR modelling OR instrument OR instruments OR tool OR tools OR nomogram) AND (treatment OR treatment[MeSH terms] OR chemotherapy OR diet OR radiation OR surgery)) NOT ("mouse model" OR "animal model")'
predictor = NCI_predict.NCI_predictor()

def relearn():
	NCI_learn.train_NCI_models()

def query_pubmed(date_range=None, relative_date=30, search_term=DEFAULT_SEARCH_STRING):
	records = get_records(date_range=date_range, relative_date=relative_date, search_term=search_term)
	parsed_records = [parse_record(record) for record in records]
	relevant_records = [record for record in parsed_records if predictor.predict(record['title'], record['abstract'], record['keywords'])]
	return relevant_records


def get_records(date_range=None, relative_date=30, search_term=DEFAULT_SEARCH_STRING, retmax=100000):
	''' Yield recods for requested pubmed records '''


	# date_range: given as a tuple: (from, to) in format YYYY/MM/DD, YYYY, or YYYY/MM
	# relative_date: Returns items published within the last n number of days
	# use only one date_range or relative_date, not both


	# esearch limited to 100,000 ids at a tiem, see documentation at
	# http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
	# for more info
	if date_range:
		handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax,
								mindate=date_range[0], maxdate=date_range[1])
	else:
		handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax,
								reldate=relative_date)

	esearch_record = Entrez.read(handle)
	idlist = esearch_record["IdList"]
	handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
	records = Medline.parse(handle)
	return list(records)


def parse_record(record):
	''' Takes a record and extract desired info '''
	# display author, year, pubmed id, publication date, t

	parsed_record = {}
	parsed_record['title'] = get_record_value(record, 'TI')
	parsed_record['abstract'] = get_record_value(record, 'AB')
	parsed_record['keywords'] = get_record_value(record, 'MH') # Actually MeSH terms (!)

	parsed_record['pubdate']  = get_record_value(record, 'DP')
	parsed_record['pmid']     = get_record_value(record, 'PMID')
	return parsed_record
	


def get_record_value(record, field, missing=''):
	try:
		return record[field]
	except KeyError:
		return missing
