#!/usr/bin/env python
"""
Replace "rxn=" with "rxn__=" in CONFORMER record inline. __ is the dielectric constant in format %02d.
epsilon = 2 -> rxn02
epsilon = 4 -> rxn04
epsilon = 12 -> rxn12
epsilon = 100 -> rxn100
"""

import sys
import logging

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(levelname)-10s: %(message)s')

    if len(sys.argv) < 3:
        logging.info("Requires two command line arguments for the ftpl file and dielectric constant.")
        print("rxn-complete.py  <ftpl_file>  <dielectric constant>")
        print("Example: rxn-complete.py  tmp.ftpl 4")
        sys.exit()

    lines = open(sys.argv[1]).readlines()

    for line in lines:
        