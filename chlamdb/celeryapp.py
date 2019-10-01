from __future__ import absolute_import, unicode_literals
import os
from celery import Celery

# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings.settings')

app = Celery('chlamdb_tasks', backend='amqp', broker='amqp://chlamdb:estrella3@localhost:5672/my_vhost', result_backend='django-db')

# Using a string here means the worker doesn't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')

# Load task modules from all registered Django app configs.
app.autodiscover_tasks()

CELERYD_MAX_TASKS_PER_CHILD = 20
BROKER_TRANSPORT_OPTIONS = {'confirm_publish': True}

@app.task(bind=True)
def debug_task(self):
    print('Request: {0!r}'.format(self.request))
