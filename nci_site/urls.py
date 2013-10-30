# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 13:43:22 2013

@author: gdietz
"""

from django.conf.urls import patterns, url

from nci_site import views

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    url(r'^$', views.index, name="index"),
    url(r'^(?P<from_date>\d{8})/(?P<to_date>\d{8})/$', views.get_data),
    url(r'^export/$', views.export_pmids),
)