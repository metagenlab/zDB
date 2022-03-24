wsgi_app = "chlamdb.wsgi:application"

# The granularity of Error log outputs
loglevel = "debug"

# The number of worker processes for handling requests
workers = 2

# The socket to bind
bind = "0.0.0.0:8000"

# Write access and error info to /var/log
accesslog = errorlog = "/usr/local/gunicorn/gunicorn.log"

# Redirect stdout/stderr to log file
capture_output = True

# PID file so you can easily fetch process ID
pidfile = "/usr/local/gunicorn/gunicorn.pid"

# Daemonize the Gunicorn process (detach & enter background)
daemon = True
