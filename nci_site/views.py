#from django.http import HttpResponse
from django.shortcuts import render

from biopython_interface import query_pubmed

def index(request):
    context = {}
    return render(request, 'index.html', context)

def get_data_for_display(request, from_date, to_date):
    # dates given in YYYYMMDD format
    # need to convert to YYYY/MM/DD format
    noslashes_to_slashes = lambda noslash: '/'.join([noslash[0:4],noslash[4:6],noslash[6:8]])
    date_range = (noslashes_to_slashes(from_date), noslashes_to_slashes(to_date))
    results = query_pubmed(date_range=date_range)
    return parsed_results_to_json(results)

def parsed_results_to_json(results):
    ''' results is a list of dictionaries containing fields and values from relevant records in the query to pubmed '''
    pass
    # continue later