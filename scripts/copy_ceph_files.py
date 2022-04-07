import os
import glob


in_dir 	= "/ceph/cms/store/user/fsetti/c++_looper_ul_output/30Mar2022_local_data_skim/Data/"
out_dir	= "/ceph/cms/store/user/fsetti/c++_looper_ul_output/23Mar2022/Data/"

files = glob.glob(in_dir+'/*.root')
for file in files[200:]:
	file_name	= file.split("/")[-1]
	file_name = file_name.split("_201")[0] + '_missingdata_201' + file_name.split("_201")[1] 
	print("Copying %s into %s/%s"%(file,out_dir,file_name))
	os.system("cp %s %s/%s"%(file,out_dir,file_name))
	#print ("cp %s %s/%s"%(file,out_dir,file_name))
