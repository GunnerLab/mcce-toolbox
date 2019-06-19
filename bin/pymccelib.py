#!/usr/bin/env python
#import numpypy
import sys
import os
import numpy as np
import shutil
import random
import gzip
import math

Delta_PW_warning = 0.1
ROOMT = 298.15
PH2KCAL = 1.364
KCAL2KT = 1.688
KJ2KCAL = 0.239

class Env:
    def __init__(self):
        # Hard coded values
        self.runprm = "run.prm"
        self.version = "PyMCCE 1.0"
        self.fn_conflist1 = "head1.lst"
        self.fn_conflist2 = "head2.lst"
        self.fn_conflist3 = "head3.lst"
        self.energy_table = "energies"
        self.mc_states = "microstates"
        self.prm = self.load_runprm()
        self.tpl = {}
        self.read_extra()
        return

    def load_runprm(self):
        float_values = ["EPSILON_PROT", "TITR_PH0", "TITR_PHD", "TITR_EH0", "TITR_EHD", "CLASH_DISTANCE",
                        "BIG_PAIRWISE", "MONTE_T", "MONTE_REDUCE"]
        int_values = ["TITR_STEPS", "MONTE_RUNS", "MONTE_TRACE", "MONTE_NITER", "MONTE_NEQ",
                      "MONTE_NSTART", "MONTE_FLIPS", "NSTATE_MAX", "MONTE_NEQ"]
        prm = {}
        print("   Loading %s" % self.runprm)
        lines = open(self.runprm).readlines()
        # Sample line: "t        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)"
        for line in lines:
            line = line.strip()
            line = line.split("#")[0]  # This cuts off everything after #
            left_p = line.rfind("(")
            right_p = line.rfind(")")
            if left_p > 0 and right_p > left_p + 1:
                key = line[left_p + 1:right_p]
                fields = line[:left_p].split()
                if len(fields) >= 1:
                    value = fields[0]
                    if key in float_values:
                        prm[key] = float(value)
                    elif key in int_values:
                        prm[key] = int(value)
                    else:
                        prm[key] = value
        return prm

    def print_runprm(self):
        for key in self.prm.keys():
            print("%-25s:%s" % (key, str(self.prm[key])))
        return

    def load_ftpl(self, file):
        """Load a tpl file."""
        float_values = ["EXTRA", "SCALING"]
        int_values = []

        print("   Loading ftpl file %s" % file)
        lines = open(file).readlines()
        for line in lines:
            line = line.split("#")[0]
            fields = line.split(":")
            if len(fields) != 2:
                continue

            key_string = fields[0].strip()
            keys = key_string.split(",")
            keys = [x.strip().strip("\"") for x in keys]
            keys = [x for x in keys if x]
            keys = tuple(keys)

            value_string = fields[1].strip()
            if keys[0] in float_values:
                self.tpl[keys] = float(value_string)
            elif keys[0] in int_values:
                self.tpl[keys] = int(value_string)
            else:
                self.tpl[keys] = value_string

        return


    def load_tpl(self, file):
        """Load a tpl file."""
        print("   Loading tpl file %s" % file)
        float_values = ["EXTRA", "SCALING"]
        int_values = []

        lines = open(file).readlines()
        for line in lines:
            line = line.split("#")[0]
            if len(line) < 21:
                continue
            keys = [line[:9], line[9:14], line[15:19]]
            value_string = line[20:].strip()

            keys = [x for x in keys if x]
            keys = tuple(keys)

            if keys[0] in float_values:
                self.tpl[keys] = float(value_string)
            elif keys[0] in int_values:
                self.tpl[keys] = int(value_string)
            else:
                self.tpl[keys] = value_string

        return


    def read_extra(self):
        """Read extra.tpl."""
        fname = self.prm["EXTRA"]

        print("   Extra tpl parameters in file %s" % fname)
        if os.path.isfile(fname):
            if fname[-5:] == ".ftpl":
                self.load_ftpl(fname)
            elif fname[-4:] == ".tpl ":
                self.load_tpl(fname)

        default_values_keys = [("SCALING", "VDW0"),
                               ("SCALING", "VDW1"),
                               ("SCALING", "VDW"),
                               ("SCALING", "TORS"),
                               ("SCALING", "ELE"),
                               ("SCALING", "DSOLV")]
        for element in default_values_keys:
            if element not in self.tpl:
                print("      Set to default: %s = 1.0" % ",".join(element))
                self.tpl[element] = 1.0

        return

    def print_scaling(self):
        """Print scaling factors."""
        # print self.param
        print("      Scaling factors:")
        print("      VDW0  = %.3f" % self.tpl[("SCALING", "VDW0")])
        print("      VDW1  = %.3f" % self.tpl[("SCALING", "VDW1")])
        print("      VDW   = %.3f" % self.tpl[("SCALING", "VDW")])
        print("      TORS  = %.3f" % self.tpl[("SCALING", "TORS")])
        print("      ELE   = %.3f" % self.tpl[("SCALING", "ELE")])
        print("      DSOLV = %.3f" % self.tpl[("SCALING", "DSOLV")])
        return


class Conformer:
    def __init__(self, fields):
        # directly from head3.lst
        self.iConf = int(fields[0])
        self.confname = fields[1]
        self.flag = fields[2].lower()
        self.on = False
        self.occ = float(fields[3])
        self.crg = float(fields[4])
        self.em0 = float(fields[5])
        self.pk0 = float(fields[6])
        self.ne = int(fields[7])
        self.nh = int(fields[8])
        self.vdw0 = float(fields[9]) * env.tpl[("SCALING", "VDW0")]
        self.vdw1 = float(fields[10]) * env.tpl[("SCALING", "VDW1")]
        self.tors = float(fields[11]) * env.tpl[("SCALING", "TORS")]
        self.epol = float(fields[12]) * env.tpl[("SCALING", "ELE")]
        self.dsolv = float(fields[13]) * env.tpl[("SCALING", "DSOLV")]
        self.extra = float(fields[14])
        self.history = fields[15]
        # needed by MC process
        self.E_self = 0.0  # self energy in head3.lst
        self.E_self_mfe = 0.0  # self energy including pairwise contribution from fixed residues
        return

    def printme(self):
        print("%05d %s %c %4.2f %6.3f %5d %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %s" % (self.iConf,
                                                                                                   self.confname,
                                                                                                   self.flag,
                                                                                                   self.occ,
                                                                                                   self.crg,
                                                                                                   self.em0,
                                                                                                   self.pk0,
                                                                                                   self.ne,
                                                                                                   self.nh,
                                                                                                   self.vdw0,
                                                                                                   self.vdw1,
                                                                                                   self.tors,
                                                                                                   self.epol,
                                                                                                   self.dsolv,
                                                                                                   self.extra,
                                                                                                   self.history))


class MC_Protein:
    """Monte Carlo Protein data structure."""

    def __init__(self):
        print("\n   Reading and interpreting input energy and conformer list.")
        self.head3list, self.confnames = self.read_head3list()
        self.pairwise = self.read_pairwise()
        self.fixed_conformers, self.free_residues, self.biglist = self.group_conformers()
        self.report_residues()
        return

    def read_head3list(self):
        head3list = []
        fname = env.fn_conflist3
        print("      Loading confomer self energy from %s" % fname)

        lines = open(fname).readlines()
        lines.pop(0)
        for line in lines:
            fields = line.split()
            if len(fields) >= 16:
                conf = Conformer(fields)
                if conf.flag == "t":
                    conf.on = False
                else:
                    conf.on = True
                head3list.append(conf)

        # validate
        confnames = [x.confname for x in head3list]
        for name in confnames:
            if len(name) != 14:
                print("      ERROR: %s is not a conformer name.")
                sys.exit()
            occurrence = confnames.count(name)
            if occurrence > 1:
                print("      ERROR: Conformer %s occurred %d times" % (name, occurrence))
                sys.exit()
        return head3list, confnames

    def print_headlist(self):
        for conf in self.head3list:
            print("%05d %s %c %4.2f %6.3f %5d %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %s" % (conf.iConf,
                                                                                                       conf.confname,
                                                                                                       conf.flag,
                                                                                                       conf.occ,
                                                                                                       conf.crg,
                                                                                                       conf.em0,
                                                                                                       conf.pk0,
                                                                                                       conf.ne,
                                                                                                       conf.nh,
                                                                                                       conf.vdw0,
                                                                                                       conf.vdw1,
                                                                                                       conf.tors,
                                                                                                       conf.epol,
                                                                                                       conf.dsolv,
                                                                                                       conf.extra,
                                                                                                       conf.history))
        return

    def read_pairwise(self):
        """Read pairwise interactions from opp files in folder."""
        folder = env.energy_table
        print("      Loading pairwise interactions from opp files in folder %s ..." % folder)
        n_size = len(self.confnames)
        pairwise = np.zeros((n_size, n_size))
        for i in range(n_size):
            conf = self.head3list[i]
            oppfile = "%s/%s.opp" % (folder, conf.confname)
            resid_i = conf.confname[:3] + conf.confname[5:11]
            if os.path.isfile(oppfile):
                lines = open(oppfile)
                for line in lines:
                    fields = line.split()
                    if len(fields) < 6:
                        continue
                    confname = fields[1]
                    j = self.confnames.index(confname)
                    if j < 0:
                        print("      Warning: %s in file %s is not a conformer" % (confname, oppfile))
                        continue

                    resid_j = confname[:3] + confname[5:11]
                    if resid_i != resid_j:  # not within a residue
                        ele = float(fields[2])
                        vdw = float(fields[3])
                        pw = ele * env.tpl[("SCALING", "ELE")] + vdw * env.tpl[("SCALING", "VDW")]
                        pairwise[i][j] = pw

        # Average pairwise after loading
        for i in range(n_size - 1):
            for j in range(i + 1, n_size):
                pw1 = pairwise[i][j]
                pw2 = pairwise[j][i]
                if abs(pw1 - pw2) > Delta_PW_warning:
                    print("         Warning: big pairwise difference between %s: %.3f and %s: %.3f" % (
                        self.confnames[i],
                        pw1,
                        self.confnames[j],
                        pw2))
                pairwise[i][j] = pairwise[j][i] = (pw1 + pw2) / 2

        return pairwise

    def group_conformers(self):
        fixed_conformers = []
        free_residues = []
        residue_ids = []
        for confname in self.confnames:
            resid = confname[:3] + confname[5:11]
            if resid not in residue_ids:
                residue_ids.append(resid)
        self.residues = [[] for i in range(len(residue_ids))]  # residue stores indices to conformers
        for i in range(len(self.confnames)):
            confname = self.confnames[i]
            resid = confname[:3] + confname[5:11]
            index = residue_ids.index(resid)
            self.residues[index].append(i)

        # Verify head3list flag and occ; Find free and fixed residues
        # if total occ of "t" flagged conformers is 1:
        # the rest conformers will be set to "t 0.00", and this residue is "fixed"
        # else if total occ of "t" flagged conformers is 0:
        # if only one conformer is left:
        # the lone conformer is set to "t 1.00", and this conformer and this residue will be "fixed"
        # else:
        # this residue is "free" and occ of conformers is 0.
        # otherwise:
        #    partial fixed occupancy not allowed

        print("      Grouping conformers ...")
        # Verify flags
        # Group conformers
        for res in self.residues:
            socc = 0.0
            n_freeconf = len(res)
            for i in res:
                if not self.head3list[i].on:
                    socc += self.head3list[i].occ
                    n_freeconf -= 1
                elif abs(self.head3list[i].occ) > 0.001:  # free residue has non-0 occ
                    print("         %s %c %4.2f -> %s f  0.00 (free conformer initial occ = 0)" % (
                        self.head3list[i].confname,
                        self.head3list[i].flag,
                        self.head3list[i].occ, self.head3list[i].confname))
                    self.head3list[i].occ = 0.0
            if abs(socc - 1.0) < 0.001:  # total occ of fixed conformers are 1.0
                for i in res:
                    fixed_conformers.append(i)
                    if self.head3list[i].on:
                        print("         %s %c %4.2f -> %s t  0.00 (fixed conformers already have occ 1.0)" % (
                            self.head3list[i].confname,
                            self.head3list[i].flag, self.head3list[i].occ, self.head3list[i].confname))
                        self.head3list[i].occ = 0.0
                        self.head3list[i].on = False
                        self.head3list[i].flag = "t"
            elif abs(socc) < 0.001:  # total occ is 0
                if n_freeconf == 1:
                    for i in res:
                        if self.head3list[i].on:
                            print("         %s %c %4.2f -> %s t  1.00 (single conformer of the residue)" % (
                                self.head3list[
                                    i].confname,
                                self.head3list[
                                    i].flag,
                                self.head3list[i].occ,
                                self.head3list[
                                    i].confname))
                            self.head3list[i].on = False
                            self.head3list[i].occ = 1.0
                            self.head3list[i].flag = "t"
                            fixed_conformers.append(i)
                            break  # because only one "f"
                else:
                    free_conformers = []
                    for i in res:
                        if not self.head3list[i].on:
                            fixed_conformers.append(i)
                        else:
                            free_conformers.append(i)
                    free_residues.append(free_conformers)
            else:  # total occ is neither 0 or 1
                print("      Error: Total residue occupancy is %.2f, 0.00 or 1.00 expected." % socc)
                for i in res:
                    self.head3list[i].printme()
                print("      Exiting ...")
                sys.exit()

        # Make big list. A big list is the size of free residues. It contains other free residue index numbers that
        # have big interactions
        bigpw = env.prm["BIG_PAIRWISE"]
        biglist = [[] for i in range(len(free_residues))]
        for ir in range(len(free_residues)):
            for jr in range(ir + 1, len(free_residues)):
                next_jr = False
                for ic in free_residues[ir]:
                    if next_jr:
                        break
                    for jc in free_residues[jr]:
                        if next_jr:
                            break
                        pw = self.pairwise[ic][jc]
                        if abs(pw) >bigpw:
                            biglist[ir].append(jr)
                            biglist[jr].append(ir)
                            next_jr = True

        return fixed_conformers, free_residues, biglist

    def update_energy(self, T=298.15, ph=7.0, eh=0.0):
        # get self energy
        for ic in range(len(self.head3list)):
            conf = self.head3list[ic]
            E_ph = T / ROOMT * conf.nh * (ph - conf.pk0) * PH2KCAL
            E_eh = T / ROOMT * conf.ne * (eh - conf.em0) * PH2KCAL / 58.0
            self.head3list[
                ic].E_self = conf.vdw0 + conf.vdw1 + conf.epol + conf.tors + conf.dsolv + conf.extra + E_ph + E_eh

            # mfe from fixed conformer
            mfe = 0.0
            for jc in self.fixed_conformers:
                mfe += self.pairwise[ic][jc] * self.head3list[jc].occ

            self.head3list[ic].E_self_mfe = self.head3list[ic].E_self + mfe

    def report_biglist(self):
        fname = "biglist.info"
        lines = ["iRes iRes_with_big_interactions\n"]
        for ires in range(len(self.biglist)):
            biglist = ",".join(["%d" % x for x in self.biglist[ires]])
            if biglist:
                lines.append("%4d %s\n" % (ires, biglist))
        open(fname, "w").writelines(lines)
        return

    def report_residues(self):
        fname = "fixed_conformers.info"
        lines = ["iConf CONFORMER     FL  occ    crg ne nH\n"]
        for ic in self.fixed_conformers:
            conf = self.head3list[ic]
            lines.append("%5d %s %s %4.2f %6.3f %2d %2d\n" % (ic, conf.confname, conf.flag, conf.occ,
                                                               conf.crg, conf.ne, conf.nh))
        open(fname, "w").writelines(lines)

        fname = "free_residues.info"
        lines = ["iRes iConf CONFORMER     FL    crg ne nH\n"]
        ires = 0
        for res in self.free_residues:
            for ic in res:
                conf = self.head3list[ic]
                lines.append("%4d %5d %s %s %6.3f %2d %2d\n" % (ires, ic, conf.confname, conf.flag,
                                                               conf.crg, conf.ne, conf.nh))
            lines.append("%s\n" % ("."*35))
            ires += 1

        open(fname, "w").writelines(lines)
        return

def mc_prepdir():
    # prepare mc folder
    if os.path.exists(env.mc_states):
        if os.path.isdir(env.mc_states):
            shutil.rmtree(env.mc_states)
        else:
            os.remove(env.mc_states)

    os.mkdir(env.mc_states)

    return


def mc_sample(prot, T=298.15, ph=7.0, eh=0.0):
    print("   Titration at T = %.2f, ph = %5.2f and eh = %.0f mv" % (T, ph, eh))

    b = -KCAL2KT / (T / ROOMT)
    n_free = len(prot.free_residues)
    nflips = env.prm["MONTE_FLIPS"]

    # get ph and eh patched self energy
    prot.update_energy(T=T, ph=ph, eh=eh)

    # loop independent runs
    n_conf = sum([len(x) for x in prot.free_residues])
    runs = env.prm["MONTE_RUNS"]
    for i in range(runs):
        fname = "ph%.1f-eh%.0f-run%02d.ms" % (ph, eh, i)
        #fh = open(fname, "w")
        fh = gzip.open("%s.gz" % fname, "wb")

        # randomize a state
        state = [random.choice(x) for x in prot.free_residues]

        # obtain a complete state
        line = "T=%f, ph=%f, eh=%f\n" % (T, ph, eh)
        fh.write(line.encode())
        E = get_state_energy(prot, state)
        line = "%.3f: %s\n" % (E, ",".join(["%d" % x for x in state]))
        fh.write(line.encode())

        # MC sampling

        for iterations in range((env.prm["MONTE_NITER"])*n_conf):
            old_state = list(state)

            # choose new state
            ires = random.randrange(n_free)
            #ires = np.random.randint(n_free)
            while True:
                new_conf = random.choice(prot.free_residues[ires])
                if new_conf != state[ires]:
                    break

            old_conf = state[ires]
            state[ires] = new_conf

            dE = prot.head3list[new_conf].E_self_mfe - prot.head3list[old_conf].E_self_mfe
            for j in range(n_free):
                dE += prot.pairwise[new_conf][state[j]] - prot.pairwise[old_conf][state[j]]

            # multiflip
            if prot.biglist[ires]:
                flip_probablity = 0.5
                flip_counter = nflips
                while flip_counter > 0:
                    if random.random() < flip_probablity:
                        iflip = random.choice(prot.biglist[ires])
                        old_conf = state[iflip]
                        new_conf = random.choice(prot.free_residues[iflip])
                        state[iflip] = new_conf

                        dE += prot.head3list[new_conf].E_self_mfe - prot.head3list[old_conf].E_self_mfe
                        for j in range(n_free):
                            dE += prot.pairwise[new_conf][state[j]] - prot.pairwise[old_conf][state[j]]

                    flip_counter -= 1
                    flip_probablity = flip_probablity / 2.0

            # evaluate
            if dE < 0.0:
                flip = True
            elif random.random() < math.exp(b*dE):
                flip = True
            else:
                flip = False

            if flip:
                new = set(state)
                old = set(old_state)
                on_confs = new - old
                off_confs = old - new
                E += dE
                line = "%.3f:" % E + ",".join(["-%d"%x for x in off_confs])+","+ ",".join(["%d"%x for x in
                                                                                          on_confs])+"\n"
                fh.write(line.encode())
            else:
                state = old_state
                fh.write("\n".encode())

        fh.close()

    return

def validate_state(prot, state):
    # each conf in state is in free_residues
    # each res in free_residues has one and only one conf in state
    # This makes sure the state is free residues only and garauntees correct energy
    counters = [0 for x in prot.free_residues]  # on conf occurance in each residue, should be all 1
    outsiders = []  # on conf not in free residues, should be empty

    for ic in state:
        found = False
        for ir in range(len(prot.free_residues)):
            if ic in prot.free_residues[ir]:
                counters[ir] += 1
                found = True
        if not found:
            outsiders.append(ic)

    matched = True
    for ir in range(len(counters)):
        if counters[ir] == 0:
            i_1stconf = prot.free_residues[ir][0]
            confname = prot.confnames[i_1stconf]
            resid = confname[:3] + confname[5:11]
            print("Free residue %s doesn't have any on-conformer in microstate." % resid )
            matched = False
        elif counters[ir] > 1:
            i_1stconf = prot.free_residues[ir][0]
            confname = prot.confnames[i_1stconf]
            resid = confname[:3] + confname[5:11]
            print("Free residue %s has multiple on-conformers in microstate." % resid )
            matched = False

    return matched


def get_state_energy(prot, state):
    E = 0.0

    # all fixed self energy
    for ic in prot.fixed_conformers:
        E += prot.head3list[ic].E_self_mfe * prot.head3list[ic].occ

    # minus one side of pw fixed to fixed
    n_fixed_conformers = len(prot.fixed_conformers)
    for i in range(n_fixed_conformers -1):
        ic = prot.fixed_conformers[i]
        for j in range(i+1, n_fixed_conformers):
            jc = prot.fixed_conformers[j]
            E -= prot.pairwise[ic][jc]*prot.head3list[ic].occ*prot.head3list[jc].occ

    # plus self on-conformers
    for ic in state:
        E += prot.head3list[ic].E_self_mfe

    # plus pw on-conformer to on-conformer
    for kc in range(len(state) - 1):
        ic = state[kc]
        for lc in range(kc+1, len(state)):
            jc = state[lc]
            E += prot.pairwise[ic][jc]

    return E


def get_state_energy_details(prot, state):
    E = 0.0

    #print("Microstate: %s" % ",".join(["%d" % x for x in state]))
    # all fixed self energy
    for ic in prot.fixed_conformers:
        E += prot.head3list[ic].E_self_mfe * prot.head3list[ic].occ
        #print("%s %.3f" % (prot.head3list[ic].confname, prot.head3list[ic].occ))

    # minus one side of pw fixed to fixed
    n_fixed_conformers = len(prot.fixed_conformers)
    for i in range(n_fixed_conformers -1):
        ic = prot.fixed_conformers[i]
        for j in range(i+1, n_fixed_conformers):
            jc = prot.fixed_conformers[j]
            E -= prot.pairwise[ic][jc]*prot.head3list[ic].occ*prot.head3list[jc].occ

    #print(state, prot.fixed_conformers)
    state = list(set(state) - set(prot.fixed_conformers))
    state.sort()
    print(state)

    # plus self on-conformers
    for ic in state:
        E += prot.head3list[ic].E_self_mfe

    # plus pw on-conformer to on-conformer
    for kc in range(len(state) - 1):
        ic = state[kc]
        for lc in range(kc+1, len(state)):
            jc = state[lc]
            E += prot.pairwise[ic][jc]

    return E


def deltaE(prot, state, off_confs, on_confs):
    """
    Calculate delta E based on conformer difference, state is not altered
    """
    dE = 0.0
    for ic in off_confs:
        dE -= prot.head3list[ic].E_self_mfe
        state = state - {ic}
        for jc in list(state):
            dE -= prot.pairwise[ic][jc]
    for ic in on_confs:
        dE += prot.head3list[ic].E_self_mfe
        for jc in list(state):
            dE += prot.pairwise[ic][jc]
        state = state.add(ic)

    return dE


env = Env()

if __name__ == "__main__":
    print("This is pymcce module.")
    print("Use pymonte.py to run mcce step 4.")
    prot = MC_Protein()
    prot.report_biglist()