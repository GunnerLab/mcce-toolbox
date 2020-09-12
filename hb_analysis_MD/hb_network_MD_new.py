import MDAnalysis
import MDAnalysis.analysis.hbonds
import pandas as pd
import numpy as np
import os
from collections import defaultdict 
import networkx as nx
import matplotlib.pyplot as plt
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
#logger.addHandler(logging.FileHandler('test.log', 'a'))
print = logger.info

sys.setrecursionlimit(1000)
print(sys.getrecursionlimit())


class HB_MD:
	def __init__(self, frame):
		self.direct_connection = pd.DataFrame(columns = ['donor_residue', 'acceptor_residue'])
		self.one_water_connection = pd.DataFrame(columns = ['donor_residue', 'acceptor_residue'])
		self.two_water_connection = pd.DataFrame(columns = ['donor_residue', 'acceptor_residue'])
		self.three_water_connection = pd.DataFrame(columns = ['donor_residue', 'acceptor_residue'])
		self.four_water_connection = pd.DataFrame(columns = ['donor_residue', 'acceptor_residue'])
		self.hb_analysis(frame, self.direct_connection, self.one_water_connection, self.two_water_connection, self.three_water_connection, self.four_water_connection)
		return 


	def addEdge(self, graph,u,v): 
	    graph[u].append(v) 
	  
	def generate_edges(self, graph): 
	    edges = [] 
	    for node in graph: 
	        for neighbour in graph[node]: 
	            edges.append((node, neighbour)) 
	    return edges 


	def find_path(self, graph, start, end, path =[]): 
		path = path + [start] 
		if start == end: 
			return path 
		for node in graph[start]: 
			if node not in path: 
				newpath = find_path(graph, node, end, path) 
			if newpath:  
				return newpath 
		return None


	def find_all_path(self, graph, start, path, paths):
		if len(path) == 6:
			return paths.append(list(path))
		if len(graph[start]) == 0:
			return paths.append(list(path))

		for node in graph[start]:
			if node in path:
				continue
			path.append(node)
			self.find_all_path(graph, node, path, paths)
			path.pop()

	def get_chain(self, frame, chain):
		i = 0
		pdb = open(frame, 'r')
		#os.system('sed -i "s/1H / H1/" hoh.pdb')
		for line in pdb:
			#line.replace('HOH', 'TIP3')
			if line[0:4] != 'ATOM':
				continue
			chain[i] = line[21:22]
			i += 1
		return 

	def MDtraj(self, pdb):
		#print('Getting coordinate')
		h3 = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(pdb, 'not resname ALA and not resname GLN and not resname GLY and not resname ILE and not resname LEU and not resname PHE and not resname PRO and not resname VAL', 
			                                                  'not resname ALA and not resname GLN and not resname GLY and not resname ILE and not resname LEU and not resname PHE and not resname PRO and not resname VAL', distance=3.5, angle=90.0, acceptors = {'O1', 'O2'})
		#print('Analyzing')
		h3.run()
		#print('Generating table')
		h3.generate_table()
		#print('Generating form')
		df3 = pd.DataFrame.from_records(h3.table)
		h3.generate_table() 
		df3 = pd.DataFrame.from_records(h3.table)
		return df3

	def get_all_connection(self, df3, chain, index_donor, index_accept):
		for index2, row2 in df3.iterrows():
			if row2['donor_resnm'] == 'TIP3'and row2['acceptor_resnm'] != 'TIP3':
				if row2['donor_atom'] == 'H1':
					index_donor.append(row2['donor_resnm'] + '_' + str(row2['donor_index']-1))
					index_accept.append(row2['acceptor_resnm'] + '_' + chain[row2['acceptor_index']] + '_' + str(row2['acceptor_resid']))
				if row2['donor_atom'] == 'H2':
					index_donor.append(row2['donor_resnm'] + '_' + str(row2['donor_index']-2))
					index_accept.append(row2['acceptor_resnm'] + '_' + chain[row2['acceptor_index']] + '_' + str(row2['acceptor_resid'])) 

			elif row2['acceptor_resnm'] == 'TIP3' and row2['donor_resnm'] != 'TIP3':
				index_accept.append(row2['acceptor_resnm'] + '_' + str(row2['acceptor_index']))
				index_donor.append(row2['donor_resnm'] + '_' + chain[row2['donor_index']] + '_' + str(row2['donor_resid']))
			elif row2['acceptor_resnm'] == 'TIP3' and row2['donor_resnm'] == 'TIP3':
				if row2['donor_atom'] == 'H1':
					index_donor.append(row2['donor_resnm'] + '_' + str(row2['donor_index']-1))
					index_accept.append(row2['acceptor_resnm'] + '_' + str(row2['acceptor_index']))
				if row2['donor_atom'] == 'H2':
					index_donor.append(row2['donor_resnm'] + '_' + str(row2['donor_index']-2))
					index_accept.append(row2['acceptor_resnm'] + '_' + str(row2['acceptor_index']))
			else:
				index_donor.append(row2['donor_resnm'] + '_' + chain[row2['donor_index']] + '_' + str(row2['donor_resid']))
				index_accept.append(row2['acceptor_resnm'] + '_' + chain[row2['acceptor_index']] + '_' + str(row2['acceptor_resid']))

		return

	def divide_networks(self, hb_two, donor_residue, acceptor_residue, donor_residue2, acceptor_residue2):
		#print('Divide networks')
		for row in range(len(hb_two)):
			if hb_two['donor_residue'][row][0:3] != 'TIP' and hb_two['acceptor_residue'][row][0:3] != 'TIP': 
				if hb_two['donor_residue'][row] == hb_two['acceptor_residue'][row]:
					continue
				else:
					donor_residue.append(hb_two['donor_residue'][row])
					acceptor_residue.append(hb_two['acceptor_residue'][row])
			else:
				if hb_two['donor_residue'][row] == hb_two['acceptor_residue'][row]:
					continue
				else:
					donor_residue2.append(hb_two['donor_residue'][row])
					acceptor_residue2.append(hb_two['acceptor_residue'][row])
		return

	def count_water_num(self, path, donor, accept, wat_num):
		#print('Count number of water in paths')
		for item in path:
			donor_column = [item[0]]
			accpt_column = []
			count = 0
			for r in range(1, len(item)):
				if item[r][0:3] != 'TIP':
					donor_column.append(item[r])
					accpt_column.append(item[r])
					wat_num.append(count)
					count = 0
				else:

					count += 1
			if len(donor_column) > len(accpt_column):
				donor_column.pop()
			else:
				accpt_column.pop()
			donor.extend(donor_column)
			accept.extend(accpt_column)
		return

	#c = u.select_atoms("protein and prop z > 85 or around 3.0 protein and prop z > 85 ")
	#c.write('/Users/zhangyingying/Dropbox (City College)/Yingying/large_file/new_trajectories_PSII_wt/cut_frame32_50_test.pdb')

	def hb_analysis(self, frame, direct_connection, one_water_connection, two_water_connection, three_water_connection, four_water_connection):
		chain = {}
		graph = defaultdict(list) 
		pdb = MDAnalysis.Universe(frame)
		self.get_chain(frame, chain)

		df3 = self.MDtraj(pdb)

		index_donor = []
		index_accept = []
		self.get_all_connection(df3, chain, index_donor, index_accept)
		df3['donor_residue'] = index_donor
		df3['acceptor_residue'] = index_accept

		dic_hdonnor = {'ASP':['HD1', 'HD2'], 'ARG': ['HH11', 'HH12', 'HH21', 'HH22', 'HE'], 'GLU':['HE1', 'HE2'], 'HIS':['HD1', 'HE2'], 'HSD':['HD1', 'HE2'], 'HSE':['HD1', 'HE2'], 'HSP':['HD1', 'HE2'],
		                'SER':['HG'], 'THR':['HG1'], 'ASN':['HD21', 'HD22'], 'GLN':['HE21', 'HE22'], 'CYS':['HG'], 'TYR':['HH'], 'TRP':['HE1'], 'LYS':['HZ1', 'HZ2', 'HZ3'], 'TIP3':['H1', 'H2'], 'HOH':['1H', '2H']}
		dic_accept = {'ASP':['OD1', 'OD2'], 'HCO': ['OC1', 'OC2'], 'ARG': ['NE', 'NH1', 'NH2'], 'GLU':['OE1', 'OE2'], 'HSD':['ND1', 'NE2'], 'HSE':['ND1', 'NE2'], 'HSP':['ND1', 'NE2'], 'HIS':['ND1', 'NE2'],
		                'SER':['OG'], 'THR':['OG1'], 'ASN':['OD1'], 'GLN':['OE1'], 'CYS':['SG'], 'TYR':['OH'], 'LYS':['NZ'], 'MET':['SD'], 'CLX':['CLX'], 'CLA':['CLA'], 'OX2':['OX2'], 'PL9':['O1', 'O2'], 'FX':['FX'], 'TIP3':['OH2'], 'HOH':['O'], 'MQ8':['O1', 'O2']}

		donor_residue_pick = []
		acceptor_residue_pick = []
		for index, row in df3.iterrows():
		    if row['donor_resnm'] in dic_hdonnor.keys() and row['acceptor_resnm'] in dic_accept.keys():
		        if row['donor_atom'] in dic_hdonnor[row['donor_resnm']] and row['acceptor_atom'] in dic_accept[row['acceptor_resnm']]:
		            donor_residue_pick.append(row['donor_residue'])
		            acceptor_residue_pick.append(row['acceptor_residue'])
		    else:
		        continue

		hb_two = pd.DataFrame({'donor_residue':donor_residue_pick, 'acceptor_residue':acceptor_residue_pick})
		
		donor_residue = []
		acceptor_residue = []
		donor_residue2 = []
		acceptor_residue2 = []
		self.divide_networks(hb_two, donor_residue, acceptor_residue, donor_residue2, acceptor_residue2)
		dire_con = pd.DataFrame({'donor_residue': donor_residue, 'acceptor_residue': acceptor_residue, 'wat_num': [0]*len(donor_residue)})
		wat_con = pd.DataFrame({'donor_residue': donor_residue2, 'acceptor_residue': acceptor_residue2})

		# connection via water
		wat_con = wat_con.drop_duplicates()
		wat_con.index = range(0, len(wat_con))

		# direct connection
		dire_con = dire_con.drop_duplicates()
		dire_con.index = range(0, len(dire_con))


		#wat_con.to_csv('/Users/zhangyingying/Dropbox (City College)/Yingying/PSII/quinone/hb_network/conncetion_hoh_frame32_50.csv')

		#print('Generating graph')
		for i in range(len(wat_con)):
			self.addEdge(graph, wat_con['donor_residue'][i], wat_con['acceptor_residue'][i])

		visited = []
		path = []
		#print('Finding all paths through water')
		for res in range(len(wat_con)):
			results = []
			if wat_con['donor_residue'][res] not in visited and wat_con['donor_residue'][res][0:3] != 'TIP':
				self.find_all_path(graph, wat_con['donor_residue'][res], [wat_con['donor_residue'][res]], results)
				path = path + results
				visited.append(wat_con['donor_residue'][res])
			else:
				continue
					
		donor = []
		accept = []
		wat_num = []
		self.count_water_num(path, donor, accept, wat_num)

		# put all the connection together get the network
		res_wat_res = pd.DataFrame({'donor_residue': donor, 'acceptor_residue': accept, 'wat_num': wat_num})
		res_wat_res = res_wat_res.drop_duplicates()
		hb_network = pd.concat([dire_con, res_wat_res])
		hb_network.index = range(0, len(hb_network))
		visited_1 = [] 
		visited_2 = []
		visited_3 = []
		visited_4 = []
		for i in range(0, len(hb_network)):
			if hb_network['wat_num'][i] == 0:
				new_row = pd.Series({'donor_residue': hb_network['donor_residue'][i], 'acceptor_residue': hb_network['acceptor_residue'][i]})
				direct_connection = direct_connection.append(new_row, ignore_index=True)

			if hb_network['wat_num'][i] <= 1 and [hb_network['donor_residue'][i], hb_network['acceptor_residue'][i]] not in visited_1:
				visited_1.append([hb_network['donor_residue'][i], hb_network['acceptor_residue'][i]])
				new_row = pd.Series({'donor_residue': hb_network['donor_residue'][i], 'acceptor_residue': hb_network['acceptor_residue'][i]})
				one_water_connection = one_water_connection.append(new_row, ignore_index=True)

			if hb_network['wat_num'][i] <= 2 and [hb_network['donor_residue'][i], hb_network['acceptor_residue'][i]] not in visited_2:
				visited_2.append([hb_network['donor_residue'][i], hb_network['acceptor_residue'][i]])
				new_row = pd.Series({'donor_residue': hb_network['donor_residue'][i], 'acceptor_residue': hb_network['acceptor_residue'][i]})
				two_water_connection = two_water_connection.append(new_row, ignore_index=True)

			if hb_network['wat_num'][i] <= 3 and [hb_network['donor_residue'][i], hb_network['acceptor_residue'][i]] not in visited_3:
				visited_3.append([hb_network['donor_residue'][i], hb_network['acceptor_residue'][i]])
				new_row = pd.Series({'donor_residue': hb_network['donor_residue'][i], 'acceptor_residue': hb_network['acceptor_residue'][i]})
				three_water_connection = three_water_connection.append(new_row, ignore_index=True)

			if hb_network['wat_num'][i] <= 4 and [hb_network['donor_residue'][i], hb_network['acceptor_residue'][i]] not in visited_4:
				visited_4.append([hb_network['donor_residue'][i], hb_network['acceptor_residue'][i]])
				new_row = pd.Series({'donor_residue': hb_network['donor_residue'][i], 'acceptor_residue': hb_network['acceptor_residue'][i]})
				four_water_connection = four_water_connection.append(new_row, ignore_index=True)

		self.direct_connection = direct_connection
		self.one_water_connection = one_water_connection
		self.two_water_connection = two_water_connection
		self.three_water_connection = three_water_connection
		self.four_water_connection = four_water_connection

		return 

if __name__ == "__main__":
	traj_index = 100

	Direct = pd.DataFrame(columns = ['donor_residue', 'acceptor_residue'])
	One_water = pd.DataFrame(columns = ['donor_residue', 'acceptor_residue'])
	Two_water = pd.DataFrame(columns = ['donor_residue', 'acceptor_residue'])
	Three_water = pd.DataFrame(columns = ['donor_residue', 'acceptor_residue'])
	Four_water = pd.DataFrame(columns = ['donor_residue', 'acceptor_residue'])

	for f in range(0, 100, 2):
		os.system('sed -i "s/HOH /TIP3/" "/archive/yzhang/PS2/quinone/md_traj/traj_%s/frame%s_%s-stripped.pdb"' % (traj_index, f, traj_index))
		os.system('sed -i "s/O   TIP/OH2 TIP/" "/archive/yzhang/PS2/quinone/md_traj/traj_%s/frame%s_%s-stripped.pdb"' % (traj_index, f, traj_index))
		frame = '/archive/yzhang/PS2/quinone/md_traj/traj_%s/frame%s_%s-stripped.pdb' % (traj_index, f, traj_index)
		hb_info = HB_MD(frame)
		Direct = pd.concat([Direct, hb_info.direct_connection])
		One_water = pd.concat([One_water, hb_info.one_water_connection])
		Two_water = pd.concat([Two_water, hb_info.two_water_connection])
		Three_water = pd.concat([Three_water, hb_info.three_water_connection])
		Four_water = pd.concat([Four_water, hb_info.four_water_connection])
		#os.remove('frame%s_%s.pdb' % (f, traj_index))
		print('Frame%s finished' % f)		

	Direct = Direct.groupby(['donor_residue', 'acceptor_residue']).size().reset_index(name='Frequency')
	One_water = One_water.groupby(['donor_residue', 'acceptor_residue']).size().reset_index(name='Frequency')
	Two_water = Two_water.groupby(['donor_residue', 'acceptor_residue']).size().reset_index(name='Frequency')
	Three_water = Three_water.groupby(['donor_residue', 'acceptor_residue']).size().reset_index(name='Frequency')
	Four_water = Four_water.groupby(['donor_residue', 'acceptor_residue']).size().reset_index(name='Frequency')

	Direct.to_csv('/archive/yzhang/PS2/quinone/md_traj/hb_csv/direct_connection_traj%s.csv' % traj_index)
	One_water.to_csv('/archive/yzhang/PS2/quinone/md_traj/hb_csv/onewater_connection_traj%s.csv' % traj_index)
	Two_water.to_csv('/archive/yzhang/PS2/quinone/md_traj/hb_csv/twowater_connection_traj%s.csv' % traj_index)
	Three_water.to_csv('/archive/yzhang/PS2/quinone/md_traj/hb_csv/threewater_connection_traj%s.csv' % traj_index)
	Four_water.to_csv('/archive/yzhang/PS2/quinone/md_traj/hb_csv/fourwater_connection_traj%s.csv' % traj_index)
	#print(Direct, len(Direct))
	#print(One_water, len(One_water))
	#print(Two_water, len(Two_water))
	#print(Three_water, len(Three_water))
	#print(Four_water, len(Four_water))
