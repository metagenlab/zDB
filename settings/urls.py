from django.conf.urls import include, url
#import autocomplete_light
#autocomplete_light.autodiscover()

from django.contrib import admin
admin.autodiscover()
import chlamdb
from django.views.generic.base import RedirectView

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^chlamdb/', include('chlamdb.urls')),
    url(r'^', include('chlamdb.urls')),
]
