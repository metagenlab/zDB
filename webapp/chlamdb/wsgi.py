"""
WSGI dispatcher that mounts the ETE web app under /ete alongside Django.
It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/howto/deployment/wsgi/
"""

import os

from django.core.wsgi import get_wsgi_application
from ete4.smartview import explorer as _explorer
from werkzeug.middleware.dispatcher import DispatcherMiddleware

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings.settings")
django_app = get_wsgi_application()

# `default_app` returns the Bottle WSGI application instance
ete_app = _explorer.default_app()

application = DispatcherMiddleware(django_app, {
    '/ete': ete_app,
})


class RefererRedirectMiddleware:
    """WSGI middleware that redirects requests originating from ete to /ete.

    Many ETE frontend assets and API calls use absolute paths like
    `/static/...` or `/api/...`. When the explorer is mounted at `/ete`, the
    browser still requests the absolute paths. This middleware detects such
    requests (by checking `HTTP_REFERER`) and redirects them to include
    the `/ete` prefix so the DispatcherMiddleware routes them to the
    ETE app. We cannot simply reroute them as we need to also detect subsequent
    requests coming from js scripts downloaded by ete.
    """

    def __init__(self, app, ete_app, prefixes=None):
        self.app = app
        self.ete_app = ete_app
        self.prefixes = prefixes or ('/static/', '/api', '/trees', '/help')
        self.mount_prefix = "/ete"

    def __call__(self, environ, start_response):
        path = environ.get('PATH_INFO', '') or ''
        # If the request was referred from any ETE page/resource, treat it
        # as belonging to the mounted explorer and dispatch to the ETE app.
        referer = (environ.get('HTTP_REFERER') or '').lower()
        if 'ete' in referer and not path.startswith(self.mount_prefix):
            # Respond with a redirect to the mounted explorer path so the
            # browser will request the correct /ete-prefixed URL instead of
            # attempting absolute-root requests that miss the mount point.
            location = self.mount_prefix + path
            qs = environ.get('QUERY_STRING', '')
            if qs:
                location = f"{location}?{qs}"
            start_response('302 Found', [('Location', location)])
            return [b'']
    
        return self.app(environ, start_response)


# Wrap the composed application so requests to ETE prefixes are handled by
# the ETE app without modifying site-packages.
application = RefererRedirectMiddleware(application, ete_app)
