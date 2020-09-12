# Hydrogen bond analysis for MD frames
This folder is used for analysing the hrdrogen bond network in MD frames. Three python scripts are included.   
- ```grep_frame.py``` is used for extracting the frame from MD trajectory.   
- ```stripoff_solvent.py``` is used for removing all the water molecules outside of the protein.  
- ```hb_network_MD_new.py``` is used for the hb analysis for the frame you choose.   

## Howto use
### 1. First step is using ```grep_frame.py```.
Input file: MD trajectory  (traj.dcd)  
Outout file: MD snapshot (frameXX.pdb)   

### 2. ```stripoff_solvent.py```
Input file: MD snapshot
Output file: pdb without lipid outside protein

### 3. ```hb_network_MD_new.py```
Inputfile: pdb file from step2
Outputfile: Direct.to_csv One_water.to_csv Two_water.to_csv Three_water.to_csv Four_water.to_csv    
(You could use the ```HB_Analysis_MD_traj.ipynb``` to use ```hb_network_MD_new.py``` step by step)   
            
