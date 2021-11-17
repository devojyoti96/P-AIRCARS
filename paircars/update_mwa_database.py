import os,pip

def install(package):
    if hasattr(pip, 'main'):
        pip.main(['install', package])
    else:
        pip._internal.main(['install', package])

packages=["extension-helpers","numpy==1.19.0"]
for package in packages:
	install(package)
datadir = os.path.dirname(__file__)
import numpy as np,datetime as dtt,json,urllib.request,copy,time
def update_mwa_obsids(obsid_file='',verbose=False,force=False):
	'''
	Function to update MWA OBSIDs 

	Parameters
	----------
	obsid_file : str 
		Name of the file to save MWA OBSIDs
	verbose : bool
		Verbose output
	force :  bool
		Update forcefully
	Returns
	-------
	str
		OBSID file name
	int
		Update success or failure code (0 or 1)
	'''
	if verbose==True:
		print ('Updating local MWA OBSid file......\n')
	if obsid_file=='':
		obsid_file=datadir+'/MWA_OBSids'
	BASEURL='http://ws.mwatelescope.org/'
	temp_array=np.empty(0,dtype='int')
	if os.path.isfile(obsid_file+'.npy')==True and force==False:
		try:
			obsids=np.load(obsid_file+'.npy',allow_pickle=True)
			start_obsid=np.max(obsids)
			temp_array=np.append(temp_array,obsids)
		except:
			start_obsid=972654120
	else:
		start_obsid=972654120
	end_obsid=3786480018  # Till 2100-01-01
	searchurl=BASEURL+'metadata/find?maxtime='+str(end_obsid)+'&page=20000000000000'
	future_search_url=BASEURL+'metadata/find?maxtime='+str(end_obsid)+'&future=on'
	try:
		end_obsid=json.load(urllib.request.urlopen(future_search_url,timeout=10))[0][0]
		if verbose==True:
			print ('Last OBSID in MWA metadata server : '+str(end_obsid)+'\n')
		while True:
			searchurl=BASEURL+'metadata/find?mintime='+str(start_obsid)+'&maxtime='+str(start_obsid+432000)
			try:
				OBSid=json.load(urllib.request.urlopen(searchurl,timeout=150))
				OBSid=np.array(OBSid)[:,0].astype('int')
				start_obsid=np.max(OBSid)+235
			except:
				OBSid=np.empty(0,dtype='int')
				start_obsid=start_obsid+3600
			if len(OBSid)!=0:
				temp_array=np.append(temp_array,OBSid)
			if start_obsid>=end_obsid:
				break
		np.save(obsid_file,temp_array)
		if verbose==True:
			print ('Updated successfully.\n')
		os.system('rm -rf casa*log')
		return obsid_file+'.npy',0
	except Exception as e:
		if verbose==True:
			print ('Error in update : '+str(e)+'\n')
			print ('Update not successful.\n')
		os.system('rm -rf casa*log')
		return obsid_file+'.npy',1
