#from django.http import HttpResponse
from django.shortcuts import render

from biopython_interface import query_pubmed

def index(request):
    context = {}
    return render(request, 'index.html', context)