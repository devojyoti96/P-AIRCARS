import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from . import inputs

CLEARSCRATCH	=	inputs.CLEARSCRATCH											#	Clear the scratch directory
ukey			=	inputs.ukey
vkey			=	inputs.vkey
wkey			=	inputs.wkey
ugrids			=	inputs.ugrids
vgrids			=	inputs.vgrids
plotuv			=	inputs.plotuv
CONVERTFITS		=	inputs.CONVERTFITS											#	Convert FITS to binary ?	
DOFLAG			=	inputs.DOFLAG											#	Do flagging ?
FLAGMODE		=	inputs.FLAGMODE 									#	'baseline' / 'uvbin' 
BASEFLAGMEAN	=	inputs.BASEFLAGMEAN		#	[FLAGON,	tolerance_mean,	tolearnce_rms,	min fraction]]			ONLY for 'baseline'
READBACK		=	inputs.READBACK											#	Read back baselines ?
SHOWBASE		=	inputs.SHOWBASE											#	SHOW baseline stats ?
SHOWTF			=	inputs.SHOWTF											#	Show time-frequency plots ?
WRITEOUT		=	inputs.WRITEOUT											#	Write output ?
BLOCKPOW		=	inputs.BLOCKPOW											#	Power low for Block non-Gaussianity (DON'T CHANGE UNLESS YOU KNOW WHAT IT IS !)

#	----------------------------------------------------------------------------------------
#
#	Function to convert UVFITS file to binary files and write in scratch direcotory
#
#	----------------------------------------------------------------------------------------

def uvfitstobinary(data,scratchdir,ugrids,vgrids,plotuv,npols,verbose):
	
	ngroups	=	len(data)
	recids	=	np.arange(ngroups)
	
	#print('Parameters	'+str(data.parnames))
	if verbose:
		print('Total no of records	=	%d'%ngroups)
	
	uvwarr		=	np.transpose(np.array([data.par(ukey),data.par(vkey),data.par(wkey)]))

	usorted		=	np.argsort(uvwarr[:,0])
	uvwarrus	=	uvwarr[usorted]
	recdataus	=	data[usorted]
	recidus		=	recids[usorted]
	
	uvbins		=	[]

	ustart		=	0
	uglen	=	1 + int(ngroups/ugrids)
	
	for ug in range (0,ugrids):
	
		ustop		=	ustart + uglen
		
		if(ustop>(ngroups+1)):
			ustop	=	ngroups+1
			
		ugridarr	=	uvwarrus[ustart:ustop]
		ugsize		=	len(ugridarr)
		#print('%d		%d		%d	%d'%(ug,ugsize,ustart,ustop))
		vsorted		=	np.argsort(ugridarr[:,1])
		ugarrvs		=	ugridarr[vsorted]
	
		urecdata	=	recdataus[ustart:ustop]
		urecvs		=	urecdata[vsorted]
		
		urecid		=	recidus[ustart:ustop]
		urecidvs	=	urecid[vsorted]
	
		vstart		=	0
		vglen		=	1 + int(ugsize/vgrids)
	
		for vg in range (0,vgrids):
			vstop		=	vstart + vglen
			if(vstop>(ugsize+1)):
				vstop	=	ugsize+1
			vgridarr	=	ugarrvs[vstart:vstop]	
			vgsize		=	len(vgridarr)
			#print('%d	%d	%d		%d	%d'%(ug,vg,vgsize,vstart,vstop))
			if verbose:
				print('%d	%d		writing...'%(ug,vg))
			
			vrecdata	=	urecvs[vstart:vstop]
			vrecid		=	urecidvs[vstart:vstop]
			uvbins.append([ug,vg,vrecdata,vrecid])
		
			#print(vgridarr[100],vrecdata.par(ukey)[100],vrecdata.par(vkey)[100],vrecdata.par(wkey)[100],vrecdata.par('BASELINE')[100])
			if(plotuv):
				colsel		=	'b.'
				if((ug+vg)%2):
					colsel	=	'r.'
				plt.plot(vgridarr[:,0],vgridarr[:,1],colsel,markersize=0.5)
		
			vstart		=	vstop
		
		ustart		=	ustop
	
	if(plotuv):
		plt.xlim([-0.0001,0.0001])
		plt.ylim([-0.0001,0.0001])
		plt.show()
	
	uvbinstats	=	[]
		
	for ug in range (0,ugrids):
		for vg in range (0,vgrids):
			
			datauvg		=	uvbins[ug*vgrids+vg]
			
			datarray	=	datauvg[2].data
			#print(np.shape(datarray))
			
			uvbinstats.append([ug,vg,1,len(datarray),len(datarray[0,0,0,0])])
			
			bindetails	=	[datauvg[3],datauvg[2].par(ukey),datauvg[2].par(vkey),datauvg[2].par(wkey),datauvg[2].par('BASELINE')]
			bindetails	=	np.array(bindetails)
			bindetails	=	np.transpose(bindetails)
			
			np.savetxt(scratchdir+'uvbindetails_%d_%d.txt'%(ug,vg),bindetails,fmt='%d	%e	%e	%e	%d')
			
			fname		=	open(scratchdir+'uvbin_%d_%d.array'%(ug,vg),'wb')		
			for p in range (0,npols):
				for d in range (0,3):					
					fname.write(datarray[:,0,0,0,:,p,d].astype('float32').tobytes())
			fname.close()	
	
	uvbinstats		=	np.array(uvbinstats)
	np.savetxt(scratchdir+'uvbin_status.txt', uvbinstats, fmt='%4d	%4d	%4d	%4d	%4d')
		
	return

#	------------------------------------------------------------------------------------------------------







#	----------------------------------------------------------------------------------------
#
#	Function to convert binary files from the scratch directory to UVFITS file
#
#	----------------------------------------------------------------------------------------

def uvfitsfrombinary(data,scratchdir,ugrids,vgrids,npols,verbose):
		
	for ug in range (0,ugrids):
		for vg in range (0,vgrids):
		
			bindetails	=	np.loadtxt(scratchdir+'uvbindetails_%d_%d.txt'%(ug,vg))
			recids		=	bindetails[:,0].astype('int32')
			baselines	=	bindetails[:,4]
			
			temparr	=	np.fromfile(scratchdir+'uvbin_%d_%d_f.array'%(ug,vg),dtype='float32',count=-1,sep="")
			
			if((len(temparr)%len(recids)) or (len(temparr)%npols) or (len(temparr)%3)):
				if verbose:
					print('Length Mismatch !!!!!!! ....... exiting...')
				return 1
			
			nchan	=	int(len(temparr)/(npols*3*len(recids)))
			
			temparr	=	np.reshape(temparr,(npols,3,len(recids),nchan))
			
			for p in range (0,npols):
				for d in range (0,3):						
					#print("Records,channels		"+str(np.shape(temparr)))
				
					for l in range (0,len(recids)):
						#print(data[recids[l]].par(ukey),data[recids[l]].par(vkey),data[recids[l]].par(wkey),bindetails[l,1],bindetails[l,2],bindetails[l,3])
						np.copyto( data.data[ recids[l] ,0,0,0,:,p,d] , temparr[ p, d, l ] )	
			if verbose:
				print('Copied	%d %d	of	%d %d'%(ug,vg,ugrids,vgrids))

	return 0

#	----------------------------------------------------------------------------------------








#	----------------------------------------------------------------------------------------
#
#	Function to convert baselines to binary files and write in scratch direcotory
#
#	----------------------------------------------------------------------------------------

def baselinestobinary(ANTS,data,scratchdir,npols,verbose):
	
	blid	=	[]
	for a in range (1,ANTS):
		for b in range (a+1,ANTS+1):
			blid.append([a,b,256*a+b])
	blid	=	np.array(blid)
	nbase	=	len(blid)
	
	if verbose:
		print('\nWriting total baselines = %d'%nbase)
	
	blstatus	=	[]
	
	for bline in range (0,nbase):
	
		bl		=	blid[bline]		
		bindx	=	np.where(data.par('BASELINE')==bl[2])[0]
		bdata	=	data[bindx]
		
		if(len(bdata)==0):
			if verbose:
				print('Baseline	%d-%d		No data'%(bl[0],bl[1]))
			blstatus.append([bl[0],bl[1],0,0,0])
			continue
		#print(bdata.parnames)
		if verbose:
			print('Baseline	%d-%d		writing...'%(bl[0],bl[1]))
		bindetails	=	[bindx,bdata.par(ukey),bdata.par(vkey),bdata.par(wkey),bdata.par('BASELINE')]
		bindetails	=	np.array(bindetails)
		bindetails	=	np.transpose(bindetails)
			
		np.savetxt(scratchdir+'baselinedetails_%d_%d.txt'%(bl[0],bl[1]),bindetails,fmt='%d	%e	%e	%e	%d')
		
		bldata	=	bdata.data
		blstatus.append([bl[0],bl[1],1,len(bldata),len(bldata[0,0,0,0])])
		
		fname		=	open(scratchdir+'baseline_%d_%d.array'%(bl[0],bl[1]),'wb')
		for p in range (0,npols):
			
			tfdata	=	bldata[:,0,0,0,:,p,:]		#	Time frequency data			
			lent	=	len(tfdata)
			
			for d in range (0,3):
				for r in range (0,lent):
					fname.write(tfdata[r,:,d].astype('float32').tobytes())
		fname.close()	
			
	np.savetxt(scratchdir+'baseline_status.txt', blstatus, fmt='%4d	%4d	%4d	%4d	%4d')	
		
	return 0

#	------------------------------------------------------------------------------------------------------







#	----------------------------------------------------------------------------------------
#
#	Function to convert baselines from binary files in scratch direcotory to UVFITS
#
#	----------------------------------------------------------------------------------------

def baselinesfrombinary(ANTS,data,scratchdir,npols,verbose):
	
	blid	=	[]
	for a in range (1,ANTS):
		for b in range (a+1,ANTS+1):
			blid.append([a,b,256*a+b])
	blid	=	np.array(blid)
	nbase	=	len(blid)

	if verbose:
		print('\nReading total baselines = %d'%nbase)
	
	blstatus	=	np.loadtxt(scratchdir+'baseline_status.txt', dtype=int)
	baselineflag=	np.loadtxt(scratchdir+'badbase.list',dtype=int)
	
	for bline in range (0,nbase):
		
		bl		=	blid[bline]
		
		if(blstatus[bline,2]==0):
			if verbose:
				print('Baseline	%d-%d		No data'%(bl[0],bl[1]))
			continue
		
		bindetails	=	np.loadtxt(scratchdir+'baselinedetails_%d_%d.txt'%(bl[0],bl[1]))
		recids		=	bindetails[:,0].astype('int32')
		baselines	=	bindetails[:,4]
		
		temparr		=	np.fromfile(scratchdir+'baseline_%d_%d_f.array'%(bl[0],bl[1]),dtype='float32',count=-1,sep="")
			
		if((len(temparr)%len(recids)) or (len(temparr)%npols) or (len(temparr)%3)):
			if verbose:
				print('Length Mismatch !!!!!!! ....... exiting...')
			return 1
			
		nchan	=	int(len(temparr)/(npols*3*len(recids)))
		
		temparr	=	np.reshape(temparr,(npols,3,len(recids),nchan))	
		temp00	=	np.zeros(np.shape(temparr), dtype='float32')
				
		for p in range (0,npols):
			for d in range (0,3):
			
				#print("Records,channels		"+str(np.shape(temparr)))
				if (baselineflag[bline,3+p]):
					for l in range (0,len(recids)):
						#print(data[recids[l]].par('UU'),data[recids[l]].par('VV'),data[recids[l]].par('WW'),bindetails[l,1],bindetails[l,2],bindetails[l,3])
						np.copyto( data.data[ recids[l] ,0,0,0,:,p,d] , temparr[ p, d, l ] )
				else:
					for l in range (0,len(recids)):
						#print(data[recids[l]].par('UU'),data[recids[l]].par('VV'),data[recids[l]].par('WW'),bindetails[l,1],bindetails[l,2],bindetails[l,3])
						np.copyto( data.data[ recids[l] ,0,0,0,:,p,d] , temp00[ p, d, l ] )	
		
		del (temp00)
		del (temparr)
		
		if verbose:
			print('Copied	baseline	%d %d'%(bl[0],bl[1]))			
		
	return 0

#	------------------------------------------------------------------------------------------------------






#	-----------------------------------------------------------------------------------------------------
#
#	Function to show time-frequency plot for one single baseline
#
#	-----------------------------------------------------------------------------------------------------

def showbasecomparison(bl,data,data2,dopl,npols,verbose):
	
	flfrac	=	[]
	
	bdata	=	data[data.par('BASELINE')==bl[2]]
	
	if(len(bdata)==0):
		if verbose:
			print('Baseline	%d-%d		No data'%(bl[0],bl[1]))
		return (np.ones(2*npols))
	
	if(dopl):
		fig=plt.figure(figsize=(12,8))
		
	bldata	=	bdata.data
	
	for p in range (0,npols):
			
		tfdata	=	bldata[:,0,0,0,:,p,:]			#	Time frequency data
		tfre	=	tfdata[:,:,0]					#	Real part
		tfim	=	tfdata[:,:,1]					#	Imaginary Part
		tfw		=	tfdata[:,:,2]					#	Weight
		tfamp	=	np.sqrt(tfre*tfre+tfim*tfim)	#	Amplitude
		
		flfrac.append(1.0-float(np.count_nonzero(tfw>0.0))/np.size(tfw))
		
		if(dopl):
			ax=fig.add_subplot(npols,2,npols*p+1)
			plt.imshow(tfamp,interpolation='none',aspect='auto')
			plt.colorbar()
	
	bdata	=	data2[data.par('BASELINE')==bl[2]]
		
	bldata	=	bdata.data
	
	for p in range (0,npols):
			
		tfdata	=	bldata[:,0,0,0,:,p,:]			#	Time frequency data
		tfre	=	tfdata[:,:,0]					#	Real part
		tfim	=	tfdata[:,:,1]					#	Imaginary Part
		tfw		=	tfdata[:,:,2]					#	Weight
		tfamp	=	np.sqrt(tfre*tfre+tfim*tfim)	#	Amplitude
		
		flfrac.append(1.0-float(np.count_nonzero(tfw>0.0))/np.size(tfw))
		
		if(dopl):
			ax=fig.add_subplot(npols,2,npols*p+2)
			plt.imshow(tfamp,interpolation='none',aspect='auto')
			plt.colorbar()
		
	if verbose:	
		print('Baseline	%d-%d	flagged fractions	'%(bl[0],bl[1])),
		for p in range (0,npols):
			print('%.2f %.2f	'%(flfrac[p],flfrac[p+npols])),
		print('\n')
	
	if(dopl):
		plt.show()
	
	return (np.array(flfrac))

#	----------------------------------------------------------------------------------------------------








	
	
	

































