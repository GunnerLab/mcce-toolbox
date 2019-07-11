#!/usr/bin/python
# Complete a mcce-tpl RADIUS records with C6 and C12 parameters

import sys

vdw_parm = {"C": (2.000, 0.150),
            "H": (1.000, 0.020),
            "O": (1.600, 0.200),
            "N": (1.750, 0.160),
            "S": (2.000, 0.200),
            "X": (2.000, 0.173)
}

def complete_radius(keys, values):
    atomname = keys[2].strip().strip("\"")
    if atomname[0] == "H": # H atom
        element = "H"
    else:
        element = atomname[1]
    if not (element in vdw_parm.keys()):
        element = "X"
    vdw = vdw_parm[element]
    line = "RADIUS, %s, %s: %s, %7.3f, %7.3f\n" % (keys[1], keys[2], values[0], vdw[0], vdw[1])
    return line


if __name__ == "__main__":
    try:
        fp = open(sys.argv[1])
        lines = fp.readlines()
        fp.close()
    except IndexError:
        print("Usage:")
        print("    vdw-complete.py free_format_tpl_file")
        sys.exit()

    newlines = []
    for line in lines:
        i = line.find("#")
        nline = line[:i]
        fields = nline.split(":")
        if len(fields) != 2:
            newlines.append(line)
            continue
        keys = fields[0].strip().split(",")
        values = fields[1].strip().split(",")
        record_name = keys[0].strip().upper()
        if record_name == "RADIUS":
            nline = complete_radius(keys, values)
            newlines.append(nline)
        else:
            newlines.append(line)

    sys.stdout.writelines(newlines)