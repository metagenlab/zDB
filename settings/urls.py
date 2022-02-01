from django.urls import include, path, re_path
#import autocomplete_light
#autocomplete_light.autodiscover()

from django.contrib import admin
admin.autodiscover()
import chlamdb
from django.views.generic.base import RedirectView

urlpatterns = [
    re_path(r'^admin/', admin.site.urls),
    re_path(r'^chlamdb/', include('chlamdb.urls')),
    re_path(r'^', include('chlamdb.urls')),
]
