CONFLIST OEH        OEHBK OEH-1 OEH-2

NATOM    OEHBK      0
NATOM    OEH-1      2
NATOM    OEH-2      1

IATOM    OEH-1  O   0
IATOM    OEH-1  H   1
IATOM    OEH-2  O   0

ATOMNAME OEH-1    0  O
ATOMNAME OEH-1    1  H
ATOMNAME OEH-2    0  O

RXN      OEH-1       -25.54
RXN      OEH-2       -103.5
PKA      OEH-1        15.7
PKA      OEH-2        21.9
PROTON   OEH-1        -1
PROTON   OEH-2        -2


#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  OEH-1  O   sp3       0     H
CONNECT  OEH-1  H   s         0     O
CONNECT  OEH-2  O   ion

CHARGE   OEH-1  O    -1.70
CHARGE   OEH-1  H     0.70
CHARGE   OEH-2  O    -2.00

RADIUS   OEH    O     1.52
RADIUS   OEH    H     1.00

VDW_RAD  OEH-1  O   1.5
VDW_EPS  OEH-1  O   0.21
VDW_RAD  OEH-1  H   0.0
VDW_EPS  OEH-1  H   0.0
VDW_RAD  OEH-2  O   1.5
VDW_EPS  OEH-2  O   0.21

#-------|-----|----|--------------
#ParaNam|Res |Atom|Param/toggle
TRANS    OEH         t

