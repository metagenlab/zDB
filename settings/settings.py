# -*- coding: utf-8 -*-

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
from django.forms.utils import flatatt

BASE_DIR = os.path.dirname(os.path.dirname(__file__))


DEBUG = True

SECURE_CONTENT_TYPE_NOSNIFF = False
SECURE_BROWSER_XSS_FILTER = False
SECURE_SSL_REDIRECT = False
SESSION_COOKIE_SECURE = False
CSRF_COOKIE_SECURE = False
X_FRAME_OPTIONS = "SAMEORIGIN"


LOGIN_URL = '/chlamdb/login/'
LOGOUT_URL = '/chlamdb/logout/'
LOGIN_REDIRECT_URL = ''


DOCS_ROOT="/home/chlamdb/docs/_build/"

SECRET_KEY   = os.environ["RUN_NAME"]
RUN_NAME     = os.environ["RUN_NAME"]
NEXTFLOW_DIR = os.environ["NEXTFLOW_DIR"]
hosts = os.environ.get("ALLOWED_HOSTS", "")

BIODB_DB_PATH = NEXTFLOW_DIR+"/db/"+RUN_NAME
SEARCH_INDEX  = NEXTFLOW_DIR+"/search_index/"+RUN_NAME
BLAST_DB_PATH = NEXTFLOW_DIR + "/blast_DB/"+RUN_NAME

ALLOWED_HOSTS = hosts.split(",")


BIODB_CONF = {
    "chlamdb.db_type" : "sqlite",
    "chlamdb.db_name" : "George",
    "chlamdb.db_psswd" : ""
}


INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'chlamdb',
    'gunicorn',
    'templatetags',
    'django.contrib.admin',
    'django.contrib.sitemaps',
    'crispy_forms'
)


MIDDLEWARE = [
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.security.SecurityMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'chlamdb.middleware.MyMiddleware',
]

ROOT_URLCONF = 'settings.urls'


LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'Europe/Paris'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.6/howto/static-files/

STATIC_URL = '/assets/'


STATICFILES_DIRS = (
    os.path.join(BASE_DIR, "assets"),
)


APPEND_SLASH = True  # Ajoute un slash en fin d'URL

GOOGLE_ANALYTICS_DOMAIN = 'chlamdb.ch'
GOOGLE_ANALYTICS_PROPERTY_ID = 'UA-125948409-1'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',

        'DIRS': [os.path.join(BASE_DIR, 'templates')],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

INTERNAL_IPS = (
    '127.0.0.1',
)


BROKER_URL = 'django://'
CELERY_BROKER_URL = 'amqp://chlamdb:estrella3@localhost:5672/my_vhost'
CELERYD_STATE_DB = "/home/tpillone/celery_worker_state"
CELERY_CACHE_BACKEND = 'django-cache'
CELERY_ACCEPT_CONTENT = ['json']
CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'
