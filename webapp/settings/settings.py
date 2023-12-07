# -*- coding: utf-8 -*-

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
import secrets

from django.forms.utils import flatatt

BASE_DIR = os.path.dirname(os.path.dirname(__file__))

SECURE_CONTENT_TYPE_NOSNIFF = False
SECURE_BROWSER_XSS_FILTER = False
SECURE_SSL_REDIRECT = False
SESSION_COOKIE_SECURE = False
CSRF_COOKIE_SECURE = False
X_FRAME_OPTIONS = "SAMEORIGIN"


SECRET_KEY = secrets.token_urlsafe()

DEBUG = int(os.environ.get("DEBUG", 0))
RUN_NAME = os.environ["RUN_NAME"]
hosts = os.environ.get("ALLOWED_HOSTS", "localhost")

PREFIX = "assets"

BIODB_DB_PATH = PREFIX+"/db/"+RUN_NAME
SEARCH_INDEX = PREFIX+"/search_index/"+RUN_NAME
BLAST_DB_PATH = PREFIX+"/blast_DB/"+RUN_NAME

ALLOWED_HOSTS = hosts.split(",")

BIODB_CONF = {
    "zdb.db_type": "sqlite",
    "zdb.db_name": "George",
    "zdb.psswd": ""
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
