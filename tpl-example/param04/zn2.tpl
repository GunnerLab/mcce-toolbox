CONFLIST ZN2        ZN2BK ZN2+2 ZN2DM

NATOM    ZN2BK      0
NATOM    ZN2+2      1
NATOM    ZN2DM      0

IATOM    ZN2+2 ZN2  0

ATOMNAME ZN2+2   0  ZN2 


#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   ZN2+2      0
PKA      ZN2+2      0.0
ELECTRON ZN2+2      0
EM       ZN2+2      0.0
RXN      ZN2+2      -113.568

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  ZN2+2 ZN2  ion

#3.Atom Parameters: Partial Charges and Radii
# Radii(vdw) were collected from http://www.webelements.com/
RADIUS   ZN2   ZN2  1.39

CHARGE   ZN2+2 ZN2  2.00
#-------|-----|----|--------------
#ParaNam|Res |Atom|Param/toggle
TRANS    ZN2         t

