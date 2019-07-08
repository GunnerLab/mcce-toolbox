import os, re
def water_mc(input_file, output_file):
	with open(output_file,mode='w',encoding='utf-8') as out_file:
		with open(input_file,mode='r',encoding='utf-8') as in_file:
			lines=[l.split() for l in in_file.readlines()]
			num_w=0
			for l in lines:
				if 'HOHDM' in l[0]:
					resn=l[0][0:3]
					resid=l[0][5:10]
					dm_w=l[1]
					occ_w=1-float(dm_w)
					num_w=num_w+occ_w
					outputstr=resn+resid+'  '+str(occ_w)+'\n'
					out_file.write(outputstr)
			return str(num_w)

inputpathway='/home/cai/jsc/water_count/xray_aa3_fort38'
outputpathway='/home/cai/jsc/water_count/xray_aa3_water_mc'
#mainpathway='/home/cai/jsc/water_count/xray_aa3_fort38'

with open(outputpathway+'/water_mc', mode='w', encoding='utf-8') as wat_file:
	water_num={}
	for root,dirs,files in os.walk(inputpathway):
		for subdir in dirs:
			subdirpath=os.path.join(root,subdir)
			outputsubdir=re.sub(inputpathway,outputpathway,subdirpath) 
			os.mkdir(outputsubdir)
		for name in files:
			snap=name[8:]
			print(snap)
			inputfilepath=os.path.join(root,name)
			outputfilepath_1=re.sub(inputpathway,outputpathway,inputfilepath)
			outputfilepath=re.sub('fort.38','occ.wat',outputfilepath_1)
			print("creating %s." % (str(outputfilepath)))
			water_num[snap]=water_mc(inputfilepath,outputfilepath)
			outstr=snap+'  '+ water_num[snap]+'\n'
			wat_file.write(outstr)




