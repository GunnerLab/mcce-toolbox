#!/usr/bin/python

import sys
import string
import operator
from itertools import groupby
import math
from collections import defaultdict
import struct
import os.path
import pandas as pd
from openpyxl import load_workbook
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import re
from mfe_adv import *
from energy_adv2 import *
import copy

####READ NEW FORMAT OF MS.DAT AND RE_MS.DAT (/home/cai/mcce3.5_enum_ms)


class MSRECORD:
   def __init__(self):
      self.counter = []
      self.stdev=0.0
      self.occ=[]

conformers=[]
grouped_ms_states=[]
ms_states=[]
grouped_crgs = []
enum_flag = 0
grouped_crgseqs=[]
E_min=0.0
E_max=0.0

def read_re_ms():
	global ms_states
	global conformers, enum_flag
	lines = read_file_respect_comments("re_ms.dat")
	del lines[:3]   # remove the two lines
	if len(ms_states) > 0:
		print "WARNING: adding to non empty microstate list."
	
	#read run.prm
	first_ph()
	#read head3.lst to read charges for each conformer
	conformers = read_headlst()
	count=0
	#judge if it counter is occupancy instead of count
	if (lines[1].split()[-2] == "occ:"): 
		enum_flag =1
		print "re_ms.dat is obtained from enumeration MCCE."
	else: enum_flag =0
	
	for i in range(0, len(lines), 2):
		ms_state =MSRECORD()
		confid_list = lines[i].split()
		energies = lines[i+1].split()
#		print energies
#		print confid_list
		ms_state.confseq = confid_list
		ms_state.cumE = float(energies[2])
		ms_state.E = float(energies[5])
		ms_state.counter = float(energies[-1])
		ms_state.occ= float(energies[-1])
		ms_state.id = count
		#assign charge state for the microstate
		ms_state.crg = 0
		ms_state.confidseq=[]
		ms_state.crgseq=""
		for index in ms_state.confseq:
			ms_state.crg += conformers[int(index)].crg
			ms_state.crgseq += str(int(conformers[int(index)].crg))
			confid = conformers[int(index)].id
			ms_state.confidseq.append(confid)
		count += 1
		ms_states.append(ms_state)
	return

def read_ms():
        global ms_states
        global conformers, enum_flag

        with open("ms.dat", "rb") as f:
		#read run.prm
		first_ph()
		#read head3.lst to read charges for each conformer
		conformers = read_headlst()

		byte = f.read(4)
		n_spe = struct.unpack('1i', byte)[0]
		print "There are %i residues." % (n_spe)
		spe_list=""
		for i in range(n_spe):
			byte = f.read(8)
			res_spe= struct.unpack('8s', byte)[0]
			spe_list = str(spe_list) + str(res_spe) + ","
		print "Microstates on "+str(spe_list)
		byte = f.read(9)
		method = struct.unpack('9s', byte)[0]
		print "Microstate is obtained from %s." %(method)
		if (method == "ENUMERATE"):
			enum_flag=1
		else: 
			enum_flag=0
		count=0
		byte = f.read(2)
		while (byte != ''):
			ms_state =MSRECORD()
			confid_list=[]
			res_confid = struct.unpack('1H', byte)[0]
                        confid_list.append(res_confid)
			for i in range(n_spe-1):
				byte=f.read(2)
				res_confid = struct.unpack('1H', byte)[0]
				confid_list.append(res_confid)
			ms_state.confseq = confid_list
			byte = f.read(8)
			ms_state.cumE = struct.unpack('1d', byte)[0]
			byte = f.read(8)
                        ms_state.E = struct.unpack('1d', byte)[0]
			if (enum_flag):
				byte=f.read(8)
				#ms_state.counter=struct.unpack('1d', byte)[0]
				ms_state.occ= struct.unpack('1d', byte)[0]
			else:
				byte=f.read(4)
				ms_state.counter=struct.unpack('1i', byte)[0]
				#ms_state.occ= struct.unpack('1i', byte)[0]
			ms_state.id = count
			#assign charge state for the microstate
			ms_state.crg = 0
			ms_state.confidseq=[]
			ms_state.crgseq=""
			for index in ms_state.confseq:
				ms_state.crg += conformers[int(index)].crg
				ms_state.crgseq += str(int(conformers[int(index)].crg))
				confid = conformers[int(index)].id
				ms_state.confidseq.append(confid)
			count += 1
			ms_states.append(ms_state)
			byte=f.read(2)
		return



def group_ms():
	global grouped_ms_states
	global E_min, E_max 
	if len(ms_states) < 1:
		print "WARNING: no microstate to evaluate."
		return     # no microstates
	ms_states.sort(key=lambda x: x.confseq)
	if not enum_flag: 
		#acumulate the duplicated microstate, dulicated microstate judging standard: same ms_state.confseq

		Es=[[x.E, x.counter] for x in ms_states]
		#checkpoint
		#print(Es[0])
		E_min, E_max=min(Es, key=lambda x: x[0]), max(Es, key=lambda x: x[0])
		tot_count = sum(x[1] for x in Es)
		print "Grouping the microstates ..."

		#group the microstate by its conformer sequence
		for k, g in groupby(ms_states, key=lambda item: item.confseq):
			#new_state=copy.deepcopy(list(g)[0])
			#print(new_state.counter, new_state.confseq)
			#print(list(g)[0].counter)
			counter_list=[]
			E_tmp=[]
			E_list=[]
			i=0
			for item in g:
				if i==0:
					new_state=copy.deepcopy(item)
				i +=1
				counter_list.append(item.counter)
				E_tmp.append(item.counter * item.E)
				E_list.extend([item.E] * item.counter)
			
			#new_state=copy.deepcopy(list(g)[0])
			new_state.counter=sum(counter_list)
			new_state.E=sum(E_tmp)/new_state.counter
			new_state.stdev=np.std(E_list)
			new_state.occ=float(new_state.counter)/tot_count
			#new_state.counter= sum( item.counter for item in list(g))
			#new_state.E=sum( item.counter * item.E for item in g)/new_state.counter
			#new_state.occ=new_state.counter/tot_count
			#E_list=[]
			#new_state.stdev=np.std(E_list.extend([item.E] * item.counter) for item in g)
			grouped_ms_states.append(new_state)
	

		#coumpute standard deviation 
		#for item in seen_ms_state:
		#	item.stdev=0.0
		#	if item.counter < 1: 
		#		print "Error: counter in re_ms.dat should be larger than 1"
		#		return
		#	t=0.0			                        
		#	for ms_state in ms_states:
		#		if (ms_state.confseq == item.confseq):
		#			t += (ms_state.E - item.E) * (ms_state.E - item.E) * ms_state.counter
		#	item.stdev = math.sqrt(t/item.counter)


		

	else:  # for enumeration case, pass the ms_states to grouped_ms_states
		print "The microstates obtained from enumaratioin."
		Es=[x.E for x in ms_states]
		E_min, E_max=min(Es), max(Es)
		grouped_ms_states= ms_states  

	#grouped_ms_states.sort(key=lambda x: x.E)
	print "Microstates are grouped."
	return

def group_crg():
	global grouped_crgs
	if len(ms_states) < 1:
		print "WARNING: no grouped microstate to evaluate."
		return     # no grouped microstates
	print "grouping microstates by the total charge..."
	#group microstates with same ionization state together
	grouped_ms_states.sort(key=lambda x: (x.crg, x.E))
	for k, g in groupby(grouped_ms_states, key=lambda item: item.crg):
			grouped_crgs.append([k, list(g)])
			

	print "Ionization grouping is done."
	return grouped_crgs	


def group_crgseq():
        global grouped_crgseqs
        if len(ms_states) < 1:
                print "WARNING: no grouped microstate to evaluate."
                return     # no grouped microstates
        print "grouping microstates by ionization sequence..."
        #group microstates with same ionization sequence together
	grouped_ms_states.sort(key= lambda x:(x.crgseq, x.E))
        for k, g in groupby(grouped_ms_states, key=lambda item: item.crgseq):
			grouped_crgs.append([k, list(g)])

        print "Ionization sequence grouping is done."
        return grouped_crgseqs

def save_excel(file, sheetname):
	group_crg()

	#save summamary into file
	#save self energy and pairwise of the microstates for each ionization state
	ph = 7.0
	eh = 0.0
	delta_E=0.0
	delta_E_self=0.0
	delta_E_pw=0.0
	delta_E_mfe=0.0
	delta_E_mfe1=0.0
	
	all_states = []
	detailed_data=defaultdict(dict)
	for index, i in enumerate(grouped_crgs):
		crg=i[0]
		crg_ms=i[1][0]
		E_self_total, E_pw, E_mfe, E_self_data, E_pw_matrix, state_res_mfe=analyze_state_energy(crg_ms.confseq, ph=7.0, eh=0.0, cutoff=cutoff)
		detailed_data[crg]={"E_self": E_self_data, "E_pw": E_pw_matrix, "E_mfe": state_res_mfe}
		all_states.append([index, i[0] ,i[1][0].crgseq,i[1][0].E, E_self_total, E_pw, E_mfe, i[1][0].occ, \
			i[1][0].confseq,i[1][0].confidseq])

	all_states=pd.DataFrame(all_states,columns=["iState","Tot_crg","Crg_seq","Min_E(Kcal)","E_self(Kcal)",\
		"E_pw(Kcal)","E_mfe(Kcal)","Occ","Conf_ids","Conf_Names"])

	##wrtie data into excel file
	writer=pd.ExcelWriter(file, engine="openpyxl")
	if os.path.isfile(file):
		writer.book = load_workbook(file)

	all_states.to_excel(writer, sheet_name=sheetname)
	if detailed_data[-3]:
		written_rows=len(all_states)+3
		detailed_data[-3]["E_self"].to_excel(writer, sheet_name=sheetname, startrow=written_rows)
		written_rows +=len(detailed_data[-3]["E_self"])+3
		detailed_data[-3]["E_pw"].to_excel(writer, sheet_name=sheetname, startrow=written_rows)
		written_rows +=len(detailed_data[-3]["E_pw"])+3
		detailed_data[-3]["E_mfe"].to_excel(writer, sheet_name=sheetname, startrow=written_rows)
	if detailed_data[-4]:
		written_rows +=len(detailed_data[-3]["E_mfe"])+3
		#print(written_rows)
        	detailed_data[-4]["E_self"].to_excel(writer, sheet_name=sheetname, startrow=written_rows)
        	written_rows +=len(detailed_data[-4]["E_self"])+3
        	detailed_data[-4]["E_pw"].to_excel(writer, sheet_name=sheetname, startrow=written_rows)
        	written_rows +=len(detailed_data[-4]["E_pw"])+3
        	detailed_data[-4]["E_mfe"].to_excel(writer, sheet_name=sheetname, startrow=written_rows)
        	written_rows +=len(detailed_data[-4]["E_mfe"])+3



	writer.save()





def print_group_crg():
        group_crg()
        print "\n"
        print "================================="
        print "iState  Tot_crg  Crg_seq  Min_E(Kcal)  Occ  Conf_ids  Conf_Names"
        all_states = []
        for index, i in enumerate(grouped_crgs):
                crg=i[0]
                crg_ms=i[1]
                sorted_crg_ms=sorted(crg_ms, key=operator.attrgetter('E'))
                print "%6d%9.2f  %s%13.2f%5.2f %s %s" % (index, crg,sorted_crg_ms[0].crgseq,sorted_crg_ms[0].E, sorted_crg_ms[0].occ, sorted_crg_ms[0].confseq,sorted_crg_ms[0].confidseq)
                #draw histogram for each charge state
                Es=[x.E for x in crg_ms]
                if crg==-3: Es_3=Es
                if crg==-4: Es_4=Es
                n,bins,patches=plt.hist(Es, bins='auto',color='#0504aa', alpha=0.7, rwidth=0.85)
                plt.grid(axis='y', alpha=0.75)
                plt.xlabel('State Energy')
                plt.ylabel('Frequency')
                plt.title(struct_state+'Energy histogram at Crg '+str(int(crg)))
                plt.xlim(E_min-(E_max-E_min)/10, E_max+10+(E_max-E_min)/10)
                plt.savefig(outputpath+'/'+struct_state+'hist_'+str(int(crg))+'.png',dpi=300)
                plt.clf()

                each_state={}
                each_state["iState"]=index
                each_state["Tot_crg"]=crg
                each_state["min_E_ms"]=sorted_crg_ms[0]
                all_states.append(each_state)
        print "\n"

        #draw hists of -3 and -4 together
        plt.hist(Es_3, bins='auto',label="crg=-3", alpha=0.7, rwidth=0.85)
        plt.hist(Es_4, bins='auto',label='crg=-4', alpha=0.7, rwidth=0.85)
        plt.xlabel('State Energy')
        plt.ylabel('Frequency')
        plt.title(struct_state+'Energy histogram at Crg -3 and -4')
        plt.grid(axis='y', alpha=0.75)
        plt.legend(loc='upper right')
        plt.savefig(outputpath+'/'+struct_state+'hist_-3_-4.png',dpi=300)
        plt.clf()

        #draw hist for all microstate 
        E_tots=[x.E for x in grouped_ms_states]
        plt.hist(E_tots, bins='auto',color='#0504aa', alpha=0.7, rwidth=0.85)
        plt.xlabel('State Energy')
        plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.75)
        plt.title(struct_state+'Energy histogram for all MS')
        plt.savefig(outputpath+'/'+struct_state+'hist.png',dpi=300)
        plt.clf()


        #print out self energy and pairwise of the microstates for each ionization state
        ph = 7.0
        eh = 0.0
        delta_E=0.0
        delta_E_self=0.0
        delta_E_pw=0.0
        delta_E_mfe=0.0
        delta_E_mfe1=0.0
        for i in all_states:
                if i["Tot_crg"] != -4.0 and i["Tot_crg"] != -3.0:
                        continue
                print "iState: %d, Tot_crg: %5.2f, Occ: %5.2f" %(i["iState"],i["Tot_crg"],i["min_E_ms"].occ)
                state=i["min_E_ms"].confseq

                E_self_total, E_pw, E_mfe =analyze_state_energy(state, ph=7.0, eh=0.0, cutoff=cutoff)
                if i["Tot_crg"] == -4.0:
                        delta_E -= i["min_E_ms"].E
                        delta_E_self -= E_self_total
                        delta_E_pw -= E_pw
                        delta_E_mfe1 -=E_mfe
                elif i["Tot_crg"] == -3.0:
                        delta_E += i["min_E_ms"].E
                        delta_E_self += E_self_total
                        delta_E_pw += E_pw
                        delta_E_mfe1 += E_mfe
        delta_E_mfe=delta_E - delta_E_self - delta_E_pw
        print "Delta_E[iS(-3)-iS(-4)]: %.2f, Delta_cluster_E_self: %.2f, Delta_cluster_E_pw: %.2f, Delta_E_mfe: %.2f, Delta_E_mfe1: %.2f" \
        %(delta_E, delta_E_self, delta_E_pw, delta_E_mfe, delta_E_mfe1)
        print("\n")
	return



def print_group_crgseq(given_flag=0, given_crgseqs=[]):
	group_crgseq()
        all_crgseq_states = []
        for index, i in enumerate(grouped_crgseqs):
                crgseq=i[0]
                crgseq_ms=i[1]
                sorted_crgseq_ms=sorted(crgseq_ms, key=operator.attrgetter('E'))
                each_crgseq_state={}
                each_crgseq_state["iState"]=index
                each_crgseq_state["Crg_seq"]=crgseq
                each_crgseq_state["crg"]=sorted_crgseq_ms[0].crg
		each_crgseq_state["E"]=sorted_crgseq_ms[0].E
		#print each_crgseq_state["crg"]
                each_crgseq_state["min_E_ms"]=sorted_crgseq_ms[0]
                all_crgseq_states.append(each_crgseq_state)
        print "\n"
        print "================================="
        print "iState  Tot_crg  Crg_seq  Min_E(Kcal)  Occ  Conf_ids  Conf_Names"
	sorted_all_crgseq_states=sorted(all_crgseq_states, key=lambda i:(i["crg"], i["E"]))
	for i in sorted_all_crgseq_states:
		if i["crg"] != -4.0 and i["crg"] != -3.0: continue
		if given_flag and i["Crg_seq"] not in given_crgseqs: continue
        	print "%6d%9.2f  %s%13.2f%5.2f %s %s" % (i["iState"], i["crg"], i["Crg_seq"],i["min_E_ms"].E, i["min_E_ms"].occ, i["min_E_ms"].confseq,i["min_E_ms"].confidseq)
			
		
	print "\n"

        #print out self energy and pairwise of the microstates for each ionization state
        ph = 7.0
        eh = 0.0
        delta_E=0.0
        delta_E_self=0.0
        delta_E_pw=0.0
        delta_E_mfe=0.0
        delta_E_mfe1=0.0
        for i in sorted_all_crgseq_states:
		tot_crg=i["min_E_ms"].crg
                if tot_crg != -4.0 and tot_crg != -3.0:
                        continue
		if given_flag and i["Crg_seq"] not in given_crgseqs: continue
                print "iState: %d, Tot_crg: %5.2f, Crg_seq: %s, Occ: %5.2f" %(i["iState"], tot_crg, i["Crg_seq"], i["min_E_ms"].occ)
                state=i["min_E_ms"].confseq
                E_self_total, E_pw, E_mfe=analyze_state_energy(state, ph=7.0, eh=0.0, cutoff=cutoff)

		if given_flag:
                	if tot_crg == -4.0:
				delta_E -= i["min_E_ms"].E
				delta_E_self -= E_self_total
				delta_E_pw -= E_pw
                                delta_E_mfe1 -= E_mfe
			elif tot_crg == -3.0:
                                delta_E += i["min_E_ms"].E
                                delta_E_self += E_self_total
                                delta_E_pw += E_pw
                                delta_E_mfe1 += E_mfe
        if given_flag:
		delta_E_mfe=delta_E - delta_E_self - delta_E_pw
        	print "Delta_E[iS(-3)-iS(-4)]: %.2f, Delta_cluster_E_self: %.2f, Delta_cluster_E_pw: %.2f, Delta_E_mfe: %.2f, Delta_E_mfe1: %.2f" \
        	%(delta_E, delta_E_self, delta_E_pw, delta_E_mfe, delta_E_mfe1)
        	print("\n")
	return




if __name__ == '__main__':
	if len(sys.argv) < 2:
		print("./ms_agg.py sheet_name")
		sys.exit()
	else:
		sheetname=sys.argv[1]

	#read microstates
	if os.path.exists("./ms.dat"):
		read_ms()
	else: 
		read_re_ms()
	#print ms_states[-1].counter

	cutoff = 2.0
	print "Microstate movements number: "+str(ms_states[-1].id +1)
	outputfile='ms_E.xlsx'
	struct_state=sheetname	
	group_ms()
	save_excel(outputfile, sheetname)
	#print("Grouping ionization seqence for the cluster...")
	#print_group_crgseq(1, given_crgseqs=[unloaded_crgseq, loaded_crgseq])


