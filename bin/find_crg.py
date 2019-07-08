
# coding: utf-8

# In[1]:

import os, re


# In[19]:

def find_sign(resname):
#find out the sign of the standard charge of resname, decide if it is neutral, positive or negative
    
    #define the class list
    neutral_list=["HIS", "TYR", "CYS"]
    positive_list=["ARG", "LYS"]
    negative_list=["ASP", "GLU"]
    
    sign=[str, float]
    
    if resname in neutral_list:
        sign[0]="NEU"
        sign[1]=0.0
    elif resname in positive_list:
        sign[0]="POS"
        sign[1]=1.0
    elif resname in negative_list:
        sign[0]="NEG"
        sign[1]=-1.0
    else:
        sign[0]="NULL"
        sign[1]=0.0
    return sign


# In[40]:

def define_crg(resname, rescrg):
#define if the charge of resname is consitent with the standard charge, if yes, return 1, else return 0.

    #crg_sign judge if the charge of resname is consitent or not, or not defined
    crg_sign=True

    #define the error of charge permitted
    error=0.00
    max_range=find_sign(resname)[1]+error
    min_range=find_sign(resname)[1]-error
    
    if find_sign(resname)[0] == "NULL":
        #if find residue that is not defined, skip judging
        return crg_sign
    else:
        if rescrg >= min_range and rescrg <= max_range:
            crg_sign=True
        else:
            crg_sign=False
    
    return crg_sign
    


# In[73]:

def aberr_res(input_file,sub_dir):
#find out residues that have aberrant charge in input_file
        with open(input_file,mode='r',encoding='utf-8') as in_file:
            lines=[l.split() for l in in_file.readlines()]
            aberr_res_list=[]
            for l in lines[1:-4]:
                #print(l)
                if len(l)==2:
                    #print(l)
                    resn=l[0][0:3]
                    reschain=l[0][4]
                    resid=int(l[0][5:9])
                    #correct residue id in f2
                    if (sub_dir[0:2] == "f2" and  reschain == "A"):
                        resid=resid+11
                    res=resn+reschain+str(resid)
                    res_crg=float(l[1])
                    if ( not (define_crg(resn,res_crg))):
                        aberr_res=[str,str]
                        aberr_res[0]=res
                        aberr_res[1]=l[1]
                        aberr_res_list.append(aberr_res)
            return aberr_res_list


# In[77]:

# eg: d={"1":[ 'johnny':crg, 'surname':crg, 'smith':crg, 'age':crg], "2":['name', 'johnny', 'surname', 'ryan', 'age'],
# "3":['name', 'jakob', 'surname', 'smith', 'age']}
# Iterate through and find out how many times each key value occurs, return the key value: [feq, [key]]
def find_freq(d):
    vals={}
    for i in d.values():
        for j in set(i):                         
            vals[j]=[0,0,[]]                
    for i in d.values():
        for j in set(i):                
            vals[j][0] +=1
                                  
    for j in vals.keys():
        for h in d.keys():
            if j in d[h]:
                vals[j][2].append(h)
                vals[j][1] +=float(d[h][j])/vals[j][0]
    return vals



# inputpathway to search sum_crg.out, generate outfile at outputpathway
inputpathway='/archive/cai/aa3/xray/quick_run/step4/'
outputpathway='/archive/cai/aa3/xray/quick_run/step4/'

with open(outputpathway+'aberr_res_detail', mode='w', encoding='utf-8') as f1:
    #output file name: 
    
    with open(outputpathway+'aberr_res_sum', mode='w', encoding='utf-8') as f2:
        f2.write("path:  "+inputpathway+'\n')
        #f2.write("#snapshot:  residue lists which have abnormal charge after MC sampling"+ "\n")
        aberr_data_dict={}
        freq_data={}
        for root,dirs,files in os.walk(inputpathway):
            for name in files:
                if name=="sum_crg.out":
                   # aberr_data_set=set()
                    aberr_data_list=[]
                    aberr_res_dict={}
                    inputfilepath=os.path.join(root,name)
                    print("processing %s." % (str(inputfilepath)))
                    subdirpath=os.path.dirname(inputfilepath)
                    filename, subdir_name=os.path.split(subdirpath)
                    aberr_data_list=aberr_res(inputfilepath,subdir_name)
                    f1.write("path: "+inputfilepath+"\n")
                    f1.write("#residue sum_crg"+"\n")
                    for i in range(len(aberr_data_list)):
                        f1.write(aberr_data_list[i][0]+"  "+aberr_data_list[i][1]+"\n")
                        aberr_res_dict[aberr_data_list[i][0]]=aberr_data_list[i][1]
                        #aberr_data_set.add(aberr_data_list[i][0])
                    #aberr_data_dict[subdir_name]=list(aberr_data_set)
                    aberr_data_dict[subdir_name]=aberr_res_dict
                    #f2.write(subdir_name+":  "+str(list(aberr_data_set))+"\n")
                    #f2.write(subdir_name+":  "+str(aberr_data_dict)+"\n")
        f2.write("\n")
        f2.write("SUMMARY:"+"\n")
        f2.write("RESIDUE\tFREQUENCY\tAVE_CRG\tNOR_CRG\tCHANGE_RATIO\t[SNAP LIST]"+"\n")
        freq_data=find_freq(aberr_data_dict)
        for i in freq_data.keys():
            aberr_resn=str(i)[0:3]
            aberr_crg=find_sign(aberr_resn)[1]
            change_ratio=freq_data[i][1]-aberr_crg
            f2.write(str(i)+"\t"+ str(freq_data[i][0])+ "\t"+"%.3f"% freq_data[i][1]+"\t"+"%.2f"% aberr_crg+"\t"+"%.3f" % change_ratio+"\t"+str(freq_data[i][2])+"\n")


            
     


# In[ ]:





