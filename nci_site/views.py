from django.http import HttpResponse
from django.shortcuts import render

import json
import urlparse
import csv

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
    

def export_pmids(request):
    query_string = request.META['QUERY_STRING']
    parsed_qs = urlparse.parse_qs(query_string)
    fmt = parsed_qs['format'][0]
    
    # convert crazy parsed url format to a sane python array of arrays
    data = []
    data_keys = [key for key in parsed_qs.keys() if key[0:4]=='data']
    nrows = len(data_keys)
    for i in range(nrows):
        key = 'data[%d][]' % i
        data.append(parsed_qs[key])
    for i,row in enumerate(data):
        for j, item in enumerate(row):
            if item == "&nbsp;":
                data[i][j]=''
                
    # extract just pmids
    pmid_rows = [[row[1]] for row in data]
        
    
    print("format")
    print(fmt)
    
    if fmt=="csv":
        return csv_response(pmid_rows)
    elif fmt == "txt":
        return txt_response(pmid_rows)
    #elif format == "xlsx":
    #    return xlsx_response(pmid_rows)

def csv_response(data):
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="pmidlist.csv"'
    
    writer = csv.writer(response)
    writer.writerows(data)
    
    return response

def txt_response(data):
    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="pmidlist.txt"'
    #import pdb; pdb.set_trace();
    for row in data:
        response.write(row[0]+"\r\n")
    
    
    
    return response