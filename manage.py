#!/usr/bin/env python
import os
import sys

import nci_site.biopython_interface

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "nci.settings")

    from django.core.management import execute_from_command_line

    execute_from_command_line(sys.argv)
