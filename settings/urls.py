from django.conf.urls import include, url
#import autocomplete_light
#autocomplete_light.autodiscover()

from django.contrib import admin
admin.autodiscover()

urlpatterns = [
    url(r'^admin/', include(admin.site.urls)),
    url(r'^chlamdb/', include('chlamdb.urls')),
]
