from django.http import HttpResponse
from django.shortcuts import render

import json

from biopython_interface import query_pubmed

def index(request):
    context = {}
    return render(request, 'index.html', context)

def get_data(request, from_date, to_date):
    ''' Takes range of dates given in YYYYMMDD format and returns a json object
    to be sent to the template renderer for displaying the requested info in the
    table '''
    
    
    
    print("from_date, to_date: %s, %s" % (from_date, to_date))
    
    # need to convert dates to YYYY/MM/DD format
    noslashes_to_slashes = lambda noslash: '/'.join([noslash[0:4],noslash[4:6],noslash[6:8]])
    date_range = (noslashes_to_slashes(from_date), noslashes_to_slashes(to_date))
    print("date range: %s,%s" % date_range)
    
    records = query_pubmed(date_range=date_range)
    json_tmp = {'aaData':[],}
#                'aoColumns':[{"sTitle":"PMID"},
#                             {"sTitle":"Title"}
#                            ]
#                }
    json_tmp['aaData']=[[record['pubdate'],
                         "<a href=\"http://www.ncbi.nlm.nih.gov/pubmed/" + record['pmid']+"\"" + ">" +record['pmid']+ "</a>",
                         record['author'],
                         record['title']] for record in records]

    return HttpResponse(json.dumps(json_tmp), content_type="application/json")
    #return render(request, 'table_content.html', context)