from django.contrib import admin
from django.urls import include, re_path

admin.autodiscover()

urlpatterns = [
    re_path(r'^admin/', admin.site.urls),
    re_path(r'^chlamdb/', include('chlamdb.urls')),
    re_path(r'^', include('chlamdb.urls')),
]
