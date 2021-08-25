# -*- coding: utf-8 -*-

"""
Django settings for test_1 project.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.6/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
import os
from django.forms.utils import flatatt

BASE_DIR = os.path.dirname(os.path.dirname(__file__))

DOCS_ROOT = "/scratch/hdd/bmarquis/github-reps/chlamdb/docs/_build/html"
# DOCS_ROOT = '/home/tpillone/work/dev/metagenlab/chlamdb/docs/_build/html'

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.6/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!

import os
SECRET_KEY = "flurb-o-matic_2000" # os.environ['SECRET_KEY']

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

SECURE_CONTENT_TYPE_NOSNIFF = False
SECURE_BROWSER_XSS_FILTER = False
SECURE_SSL_REDIRECT = False
SESSION_COOKIE_SECURE = False
CSRF_COOKIE_SECURE = False

ALLOWED_HOSTS = ['0.0.0.0', '155.105.138.104']

LOGIN_URL = '/chlamdb/login/'
LOGOUT_URL = '/chlamdb/logout/'
LOGIN_REDIRECT_URL = ''
#import memcache
import simplejson
import dill


class SimplejsonWrapper(object):
    def __init__(self, file, protocol=None):
        self.file = file


    def dump(self, value):
        dill.dumps(value, self.file)

    def load(self):
        dill.loads(value, self.file)



#CACHES = {
#    'default': {
#        'BACKEND': 'django.core.cache.backends.memcached.MemcachedCache',

#        'pickler': 'SimplejsonWrapper',
#        'unpickler': 'SimplejsonWrapper',
#    }
#}

#CACHES = {
#    'default': {
#        'BACKEND': 'django_pylibmc.memcached.PyLibMCCache',
#        'LOCATION': '127.0.0.1:1211',
#        'TIMEOUT': 500,
#        'BINARY': True,
#        'OPTIONS': {  # Maps to pylibmc "behaviors"
#            'tcp_nodelay': True,
#            'ketama': True
#        }
#    }
#}

#CACHES = {'default':
#            {'BACKEND': 'django.core.cache.backends.memcached.PyLibMCCache',
#             'LOCATION': '127.0.0.1:11211',
#             'TIMEOUT': 5000,
#             'BINARY': True,
#        }
#}


CACHES = {
        "default": {
                    "BACKEND": "django_redis.cache.RedisCache",
                    "LOCATION": "redis://127.0.0.1:6379/1",
                    "OPTIONS": {
                                    "CLIENT_CLASS": "django_redis.client.DefaultClient",
                                }
                }
    }


#TEMPLATE_DEBUG


BIODB = "db_gloomy_meucci" # '2019_06_PVC' #'chlamydia_04_16' #'2017_06_29b_motile_chlamydiae' #db_sleepy_thompson
BIODB_DB_PATH = "/scratch/hdd2/acarrara/ChlamDB_pipeline/Comparative_db/annotation_pipeline_3_genomes_MENU/db/" + BIODB  + ".db" #july_26_ann_pipeline

FOLDER_PATH= "/scratch/hdd2/acarrara/ChlamDB_pipeline/Comparative_db/annotation_pipeline_3_genomes_MENU/" #useful to have access to all the subfolders of the main one where the Nf pipeline runs


BIODB_CONF = {
        "chlamdb.db_type" : "sqlite",
        "chlamdb.db_name" : "George",
        "chlamdb.db_psswd" : ""
}

# Application definition
# 'django.contrib.addn',
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
    'crispy_forms',
    'django_celery_results',
)
#    'bootstrap3',
#     'django.middleware.security.SecurityMiddleware',
#    'django.middleware.csrf.CsrfViewMiddleware',
#
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

#WSGI_APPLICATION = 'test_1.wsgi.application'


# Database
# https://docs.djangoproject.com/en/1.6/ref/settings/#databases

# TODO : will need to remove this for cleaner solution

DB_DRIVER = "sqlite"
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': "/scratch/hdd2/acarrara/ChlamDB_pipeline/Comparative_db/annotation_pipeline_3_genomes_MENU/db/db_gloomy_meucci" #os.path.join(BASE_DIR, 'db.sqlite3'),
    }
}

# Internationalization
# https://docs.djangoproject.com/en/1.6/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'Europe/Paris'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.6/howto/static-files/

STATIC_URL = '/assets/'

"""
STATICFILES_DIRS = (
    "/assets/",
)
"""



STATICFILES_DIRS = (
    os.path.join(BASE_DIR, "assets"),
)

PROJECT_DIR = "/scratch/hdd2/acarrara/bin/chlamdb/" #"/home/tpillone/work/dev/metagenlab/chlamdb/"

# STATIC_ROOT = os.path.join(PROJECT_DIR, 'assets')

#STATIC_ROOT = '/static/'

#TEMPLATE_DIRS = (
#    "/home/trestan/work/dev/django/chlamydia/templates/",
#)

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

#HOME PAGE DATABASE
INTRO="Hi,"
TITLE= "your comparative database is ready"
SUBTITLE= "A customizable comparative genomic databas (changeable title)"
