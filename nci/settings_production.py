'''
Author: George Dietz
Description: Settings for production version of nci webiste
copy this file onto settings.py for production, otherwise for developement
copy settings_dev.py on to settings.py
'''

from settings_dev import *

DEBUG = TEMPLATE_DEBUG = False
ALLOWED_HOSTS = ['predictivemodels.cebm.brown.edu']

