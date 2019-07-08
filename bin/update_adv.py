#!/usr/bin/python

import sys



if len(sys.argv) < 3:
	print "Usage: %s fort.38 head_list [exclusion_list]\n" % sys.argv[0]
	sys.exit(0)
exclusion_file=""
if len(sys.argv) > 3:
	exclusion_file = sys.argv[3]
conformer = {}

fort = sys.argv[1]
headlist = sys.argv[2]
exclusion_list=[]

if exclusion_file:
	for line in open(exclusion_file).readlines():
		line = line.strip()
		exclusion_list.append(line)

for line in open(fort).readlines():
	line = line.strip()
	if line:
		(key, value) = line.split()
		#if value in ['0.000', '1.000']:
		#	conformer[key] = float(value)

		conformer[key] = float(value)

for line in open(headlist).readlines():
	line = line.strip()
	if line:
		tmp = line.split()
		conf_name = tmp[1]
		res_name = conf_name[0:3]+conf_name[5:10]
		if res_name in exclusion_list:
			print line
			continue
		if conformer.has_key(conf_name):
			new_str = '%s t %3.2f ' % (conf_name, conformer[conf_name])	
			old_str = '%s f %s ' % (conf_name, tmp[3])	
			new_line = line.replace(old_str, new_str)
			print new_line
		else:
			print line

