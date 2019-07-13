#!/usr/bin/env python
"""Convert mcce format tpl to free format."""

# bug, single atom residue doesn't have CONNECT in mcce tpl, but should have a CONNECT record in free tpl

import sys
import logging

mccedb = {}  # parameter database in mcce format
extra_records = ["EXTRA", "SCALING"] # records that can be directly translated


def atom_consistency(conf):
    passed = False
    natom = int(mccedb["NATOM", conf, "    "])
    for i in range(natom):
        try:
            key = ("ATOMNAME", conf, "%4d" % i)
            atomname = "{:<4}".format(mccedb[key][:4])
        except:
            logging.debug("# Error in fetching number %d atom. Check ATOMNAME record of conformer %s" % (i, conf))
            return passed
        try:
            key = ("IATOM", conf, atomname)
            iatom = int(mccedb[key].strip())
        except:
            logging.debug("# Error in finding index for atom \"%s\" of conformer %s" % (atomname, conf))
            return passed
        if iatom == i:
            passed = True

    return passed


def make_atom(conf):
    natom = int(mccedb["NATOM", conf, "    "])
    lines = []
    for i in range(natom):
        key = ("ATOMNAME", conf, "%4d" % i)
        atomname = "{:<4}".format(mccedb[key][:4])
        key = ("CONNECT", conf, atomname)
        connect = mccedb[key].rstrip()
        orbital_type = connect[:9].strip()
        if len(connect[10:]) < 1:
            nconnected = 0
        else:
            nconnected = len(connect[10:])/10+1
        connected_atoms = []
        for j in range(1, nconnected+1):
            if connect[j*10:j*10+5].strip() == "LIG":
                catomname = '%4s' % ("{:<4}".format(connect[j*10+5: j*10+9]))
            else:
                serial_str = connect[j*10:j*10+5]
                serial = int(serial_str)
                catomname = '%4s' % ("{:<4}".format(connect[j*10+5: j*10+9]))
                if serial != 0:
                    catomname = " ?  "
            connected_atoms.append(catomname)
        quoted = ['"%s"' % x for x in connected_atoms]
        str_value = ", ".join(quoted).rstrip(",")
        line = "CONNECT, \"%s\", %s: %s, %s\n" % (atomname, conf, orbital_type, str_value)
        lines.append(line)

    return lines

def make_charge(conf):
    natom = int(mccedb["NATOM", conf, "    "])
    lines = []
    for i in range(natom):
        key = ("ATOMNAME", conf, "%4d" % i)
        atomname = "{:<4}".format(mccedb[key][:4])
        key = ("CHARGE", conf, atomname)
        if mccedb.has_key(key):
            charge = float(mccedb[key])
        else:
            charge = 0.0
        line = "CHARGE, %s, \"%4s\": %6.3f\n" % (conf, atomname, charge)
        lines.append(line)

    return lines

def make_radius(conf):
    natom = int(mccedb["NATOM", conf, "    "])
    lines = []
    for i in range(natom):
        key = ("ATOMNAME", conf, "%4d" % i)
        atomname = "{:<4}".format(mccedb[key][:4])
        key = ("RADIUS", conf[:3], atomname)
        if mccedb.has_key(key):
            radius = float(mccedb[key])
        else:
            radius = 0.0
        line = "RADIUS, %s, \"%4s\": %6.3f, to_be_filled, to_be_filled\n" % (conf, atomname, radius)
        lines.append(line)
        #lines = list(set(lines))

    return lines


def make_confparm(conformers):
    lines = []
    for conf in conformers:
        if conf[-2:] == "BK" or conf[-2:] == "DM": continue
        key = ("PROTON", conf, "    ")
        nH = int(mccedb[key])
        key = ("ELECTRON", conf, "    ")
        ne = int(mccedb[key])
        key = ("PKA", conf, "    ")
        pKa0 = float(mccedb[key])
        key = ("EM", conf, "    ")
        Em0 = float(mccedb[key])
        key = ("RXN", conf, "    ")
        rxn = float(mccedb[key])

        line = "CONFORMER, %s: Em0=%6.1f, pKa0=%6.2f, ne=%2d, nH=%2d, rxn=%7.3f\n" % (conf, Em0, pKa0, ne, nH, rxn)
        lines.append(line)
    return lines


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(levelname)-s: %(message)s')
    filename = sys.argv[1]
    lines = open(filename).readlines()

    extralines = []
    for line in lines:
        end = line.find("#")
        line = line[:end]
        if len(line) < 20:
            continue
        key1 = line[:9].strip()
        key2 = line[9:15].strip()
        key3 = line[15:19]
        value = line[20:]

        # Direct translate: This part directly translate mcce to free format
        if key1 in extra_records:
            line = "%s, %s, %s: %s\n" % (key1, key2, key3, value)
            extralines.append(line)
        else:
            mccedb[(key1, key2, key3)] = value

    # Collect all conformer names from the read file
    conformers = []
    for k in mccedb.keys():
        if k[0] == "CONFLIST":
            conformers += mccedb[k].split()
    logging.debug("# Detected these conformers: [%s]" % ', '.join(map(str, conformers)))

    # check consistency between ATOMNAME and IATOM
    for conf in conformers:
        if atom_consistency(conf):      # pased
            logging.debug("# Consistency test passed for ATOM records of conformer %s." % conf)
        else:
            logging.debug("# There are discrepancies in ATOM records of conformer %s shown above." % conf)

    # Make conflist
    tplout = []
    tplout.append("\n# Values of the same key are appended and separated by \",\"\n")
    residue_names = [x[:3] for x in conformers]
    residues = list(set(residue_names))
    for residue in residues:
        line = "CONFLIST, %s: " % residue
        conflist = []
        for conf in conformers:
            if conf[:3] == residue:
                conflist.append(conf)
        line += ", ".join(conflist)
        line += "\n"
        tplout.append(line)

    # Make atom records
    tplout.append("\n# Atom definition\n")
    for conf in conformers:
        tplout += make_atom(conf)

    # Make charge records
    tplout.append("\n# Atom charges\n")
    for conf in conformers:
        tplout += make_charge(conf)

    # Make radius records
    tplout.append("\n# Atom radius, dielelctric boundary radius, van der Waals radius, and energy well depth\n")
    for conf in conformers:
        tplout += make_radius(conf)

    # Make conformer parameters
    tplout.append("\n# Conformer parameters that appear in head3.lst: ne, Em0, nH, pKa0, rxn\n")
    tplout += make_confparm(conformers)

    # Make rotatable bonds
    tplout.append("\n# Rotatable bonds. The atoms extended in the bond direction will all be rotated.\n")
    lines = []
    for key in mccedb.keys():
        if key[0] == "ROTAMER":
            residue = key[1]
            value = mccedb[key]
            atom1 = "\"%s\"" % value[:4]
            atom2 = "\"%s\"" % value[5:9]
            bond = "%s - %s" % (atom1, atom2)
            line = "ROTATE, %s: %s\n" % (residue, bond)
            lines.append(line)
    tplout += lines

    # Direct translate: This part directly translate mcce to free format
    lines = []
    extra_records = ["EXTRA", "SCALING"]
    for key in mccedb.keys():
        if key[0] in extra_records:
            line = "%s, %s: %s\n" % (key[0], key[1], mccedb[key])
            lines.append(line)
    tplout += lines

    tplout += extralines

    sys.stdout.writelines(tplout)
