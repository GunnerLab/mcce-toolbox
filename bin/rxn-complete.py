#!/usr/bin/env python
"""
Complete a mcce-tpl RADIUS records with C6 and C12 parameters
"""

import sys
import logging

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(levelname)-10s: %(message)s')

    if len(sys.argv) < 3:
        logging.info("Requires two command line arguments for the ftpl file and dielectric constant.")
        print("rxn-complete.py  <ftpl_file>  <dielectric constant id>")
        print("Example: rxn-complete.py  tmp.ftpl 04")
        sys.exit()
