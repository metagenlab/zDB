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
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.6/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 's+#sfx7#1n+bi!00jze0b^5l4s-j86-%&ki$lkp@lj5#--(kk5'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

TEMPLATE_DEBUG = True

ALLOWED_HOSTS = ['*']

'''
ADMINS = (
    ('Trestan Pillonel', 'trestan.pillonel@gmail.com')
)
EMAIL_USE_TLS = False
EMAIL_HOST = 'smtp.gmail.com'
EMAIL_HOST_USER = 'trest_p@hotmail.com'
EMAIL_HOST_PASSWORD = 'agnathe3'
EMAIL_PORT = 25
'''

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




CACHES = {'default':
            {'BACKEND': 'django_pylibmc.memcached.PyLibMCCache',
             'LOCATION': '127.0.0.1',
             'TIMEOUT': 5000,
             'BINARY': True,

        }

}




TEMPLATE_DEBUG



# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'chlamdb',
    'autocomplete_light',
    'gunicorn',
    'sorting_bootstrap',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'settings.urls'

#WSGI_APPLICATION = 'test_1.wsgi.application'


# Database
# https://docs.djangoproject.com/en/1.6/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
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

PROJECT_DIR = "/webapps/chlamdb_django/chlamdb/"

STATIC_ROOT = os.path.join(PROJECT_DIR, 'assets')

#STATIC_ROOT = '/static/'

TEMPLATE_DIRS = (
    "/webapps/chlamdb_django/chlamdb/templates/",
)

APPEND_SLASH = True  # Ajoute un slash en fin d'URL
