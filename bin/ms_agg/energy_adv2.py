#!/usr/bin/env python
import pandas as pd
import numpy as np
from data_adv import *
##change pairwise output into matrix form



PW_PRINT_CUT = 0.001

def analyze_state_energy(state, ph=7.0, eh=0.0, T=ROOMT, cutoff = PW_PRINT_CUT):
    print("Environment: pH = %.2f  eh = %.f  Temperature = %.2f K" % (ph, eh, T))
    print("Microstate: %s" % (",".join(["%d" % x for x in state])))
    print("Self energy in kCal/mol:")
    print("iConf  CONFORMER    FL  occ eh0-Em0 pKa-pKa0    vdw0    vdw1    epol    tors   dsolv   extra entropy "
          "E_self*occ")

    moving_occ = [x.occ for x in head3lst]
    # Self energy
    E_self = []

    #store all self energy details
    E_self_data=[]
    for ic in state:
        conf = head3lst[ic]
        E_ph = T / ROOMT * conf.nh * (ph - conf.pk0) * PH2KCAL
        E_eh = T / ROOMT * conf.ne * (eh - conf.em0) * PH2KCAL / 58.0
        if conf.flag.upper() == "T":
            occ = moving_occ[ic] = conf.occ
        else:
            occ = moving_occ[ic] = 1.0

        E_self_ic = occ * (E_eh + E_ph + conf.vdw0 + conf.vdw1 + conf.epol + conf.tors + conf.dsolv +
                                   conf.extra +
                           conf.entropy)
        E_self.append(E_self_ic)
        E_self_data.append([ic + 1, conf.confname,conf.flag, occ, E_eh,E_ph,conf.vdw0,conf.vdw1,conf.epol,\
         conf.tors,conf.dsolv,conf.extra,conf.entropy,E_self_ic])
                                                                                                                                                                                          
        print("%05d %s %c %4.2f %7.3f  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f   %8.3f" % (ic + 1,
                                                                                                   conf.confname,
                                                                                                   conf.flag, occ, E_eh,
                                                                                                   E_ph,
                                                                                                   conf.vdw0, conf.vdw1,
                                                                                                   conf.epol, conf.tors,
                                                                                                   conf.dsolv,
                                                                                                   conf.extra,
                                                                                                   conf.entropy,
                                                                                                   E_self_ic))
    E_self_total = sum(E_self)
    E_self_data=pd.DataFrame(E_self_data, columns=["iConf",  "CONFORMER","FL","occ","eh0-Em0",\
        "pKa-pKa0", "vdw0", "vdw1","epol","tors","dsolv","extra", "entropy", "E_self*occ"])
    print("Total E_self   %96.3f" % (E_self_total))
    print("\n")

    # Pairwise energy
    E_pw = 0.0
    E_pw_matrix=np.zeros((len(state), len(state)))
    col1=[]
    for i, ic in enumerate(state):
        col1.append(head3lst[ic].confname)
        for j, jc in enumerate(state):
            E_pw_icjc = pairwise[ic][jc] * moving_occ[ic] * moving_occ[jc]
            E_pw_matrix[i][j]=E_pw_icjc
            E_pw += E_pw_icjc * 0.5    # This is because the interaction will be counted twice A <- B, B <- A
            #if abs(E_pw_icjc) > 0:
               # print("%s <- %s: %5.2f" % (head3lst[ic].confname, head3lst[jc].confname, E_pw_icjc))
    pd.set_option('display.expand_frame_repr', False)
    E_pw_matrix= pd.DataFrame(E_pw_matrix,index=col1, columns=col1)
    print E_pw_matrix

    print("Total pairwise interaction: %6.2f\n" % E_pw)

    # Mfe energy
    E_mfe = 0.0
    n_conf=len(head3lst)
    state_res = [head3lst[x].confname[5:10] for x in state]
    #print state_res
    E_mfe_matrix=[]
    #print n_conf
    header=[] 
    for i, ic in enumerate(state):
        col2=[]
        header.append(head3lst[ic].confname)
        E_mfe_matrix.append([])
        for jc in range(n_conf):
            if head3lst[jc].confname[5:10] in state_res: continue
            col2.append(head3lst[jc].confname)
            E_mfe_icjc = pairwise[ic][jc] * moving_occ[ic] * moving_occ[jc]
            E_mfe_matrix[i].append(E_mfe_icjc)
            E_mfe += E_mfe_icjc * 0.5    # This is because the interaction will be counted twice A <- B, B <- A
            if abs(E_mfe_icjc) > cutoff:
                print("%s <- %s: %5.2f" % (head3lst[ic].confname, head3lst[jc].confname, E_mfe_icjc))
    pd.set_option('display.expand_frame_repr', False)
    #print the mfe matrix
    #print pd.DataFrame(np.transpose(E_mfe_matrix),index=col2, columns=header)
    state_res_mfe = np.sum(E_mfe_matrix, axis=1)*0.5  #This is because the interaction will be counted twice A <- B, B <- A, so half the pairwise
    state_res_mfe = pd.DataFrame([state_res_mfe], index=['mfe_tot'],columns=header)
    print state_res_mfe
    print("Total mfe interaction: %6.2f" % E_mfe)



    print("Cluster energy = E_self + E_pairwise + E_mfe = %.2f + %.2f + %.2f = %.2f" %(E_self_total, E_pw, E_mfe, E_self_total+E_pw+ E_mfe))

    print("\n")
    return E_self_total, E_pw, E_mfe, E_self_data, E_pw_matrix, state_res_mfe
