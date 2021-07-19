#! /usr/bin/python3.6
activate_this = '/home/rein/megaconstellations/webapp/venv/bin/activate_this.py'
with open(activate_this) as file_:
    exec(file_.read(), dict(__file__=activate_this))

import logging
import sys
logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, '/home/rein/megaconstellations/webapp/')
from app import app as application
application.secret_key = 'torontoapikey'
