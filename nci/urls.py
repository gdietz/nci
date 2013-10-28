from django.conf.urls import patterns, include, url

from nci_site import views

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'nci.views.home', name='home'),
    # url(r'^nci/', include('nci.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),
    
    url(r'^$', views.index, name="index"),
    url(r'^(?P<from_date>\d{8})/(?P<to_date>\d{8})/$', views.get_data),
    #url(r'^20110901/20120901/', views.get_data),
)
