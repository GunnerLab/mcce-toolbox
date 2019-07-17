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

    to_be_replaced = "rxn="
    replaced_by = "rxn%02d=" % int(sys.argv[2])
    lines = open(sys.argv[1]).readlines()

    for line in lines:
        line = line.strip("\n")
        line_str = line.strip().split("#", 1)[0].strip()
        fields = line_str.split(":")
        if len(fields) == 2:
           key = fields[0].split(",")[0]
           if key == "CONFORMER":
               line = line.replace(to_be_replaced, replaced_by)

        print(line)