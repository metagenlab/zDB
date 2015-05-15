from django.conf.urls import patterns, url
from django.conf.urls import include


urlpatterns = patterns('chlamdb.views',
                       url(r'^home/([a-zA-Z0-9_]+)$', 'home', name="home"),
                       url(r'^homology/([a-zA-Z0-9_]+)$', 'homology', name="homology"),
                       url(r'^search_taxonomy/([a-zA-Z0-9_]+)$', 'search_taxonomy', name="search_taxonomy"),
                       url(r'^locusx/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)$', 'locusx', name="locusx"),
                       url(r'^locusx/([a-zA-Z0-9_]+)/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.]+)$', 'locusx', name="locusx"),
                       url(r'^blastnr/([a-zA-Z0-9_]+)/([a-zA-Z0-9_]+)$', 'blastnr', name="blastnr"),
                       url(r'^plot_region/([a-zA-Z0-9_]+)$', 'plot_region', name="plot_region"),
                       url(r'^orthogroups/$', 'orthogroups', name="orthogroups"),
                       url(r'^circos/([a-zA-Z0-9_]+)$', 'circos', name="circos"),
                       url(r'^extract_region/([a-zA-Z0-9_]+)$', 'extract_region', name="extract_region"),
                       url(r'^circos_homology/([a-zA-Z0-9_]+)$', 'circos_homology', name="circos_homology"),
                       url(r'^search/([a-zA-Z0-9_]+)$', 'search', name="search"),
                       url(r'^interpro/([a-zA-Z0-9_]+)$', 'interpro', name="interpro"),
                       url(r'^blast/([a-zA-Z0-9_]+)$$', 'blast', name="blast"),
                       url(r'^mummer/([a-zA-Z0-9_]+)$$', 'mummer', name="mummer"),
                       url(r'^extract/([a-zA-Z0-9_]+)$$', 'extract', name="extract"),
                       url(r'^motif_search/([a-zA-Z0-9_]+)$$', 'motif_search', name="motif_search"),
                       url(r'^primer_search/([a-zA-Z0-9_]+)$$', 'primer_search', name="primer_search"),
                       url(r'^choose_db/$', 'choose_db', name="choose_db"),
                       url(r'^circos2genomes/([a-zA-Z0-9_]+)/$', 'circos2genomes', name="circos2genomes"),
                       url(r'^crossplot/$', 'crossplot', name="crossplot"),
                       url(r'^login/$', 'chlamdb_login', name='chlamdb_login'),
                       #url(r'^chaining/', include('smart_selects.urls')),
                       url(r'^autocomplete/', include('autocomplete_light.urls')),
                        url(r'^logout/$', 'logout_view', name="logout_view"),

)


urlpatterns += patterns('',
  (r'^logout/$', 'django.contrib.auth.views.logout', {'next_page': '/'}),
)
