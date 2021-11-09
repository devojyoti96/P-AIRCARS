import sys,os,numpy as np,glob,time,copy
from astropy.io import fits
from casatools import ms as mstools,quanta,msmetadata
from casatasks import *
os.system('rm -rf casa*log')
#Code is written by Leonid Benkevitch and Surajit Mondal
#Code is modified for python3 and CASA6 by Devojyoti Kansabanik, 05 Jan, 2021

def apply_acorr(dat, ants1, ants2, idx_ant, elen, tile, uvw, BW, fout):
    '''
	Amplitude correction for the beamformer-to-receiver cable length differences ('decorrelation')

	Parameters
	----------
	dat : numpy.array
		Data array
	ants1 : numpy.array
		Antenna 1 array
	ants2 : numpy.array
		Antenna 2 array
	idx_ant : numpy.array
		Antenna indices
	elen : numpy.array
		Electronic length of the tiles
	tile : numpy.array
		Tile names
	uvw : numpy.array
		UVW array
	BW : float
		Channel width
	fout : str
		File to wrire correction factors
	Returns
	-------
	numpy.array
		Amplitude corrected data array
    '''
    c_speed_of_light = 299792458. # In m/s
    tauw_sign = +1.
    iau = np.where(ants1 == ants2)[0] # Pointers into ants where a1 == a2
    uant = ants1[iau]  # Unique antennas in the same order as they appear
    nant = len(uant)                  # Number of unique antennas
    antmax = ants1.max()              # Max antenna numeric *name*
    #
    # Fill in the array of antenna indices idx_ant
    #
    iuant = 0
    for ia in uant:
        idx_ant[ia] = iuant
        iuant = iuant + 1

    #
    # Baseline u,v,w coordinates and lengths
    #
    u = uvw[0,:]
    v = uvw[1,:]
    w = uvw[2,:]
    blen = np.sqrt(u**2 + v**2)

    #
    # Amplitude-corrected data block
    #
    cdat = np.zeros_like(dat)

    ndat = len(dat[0,0,:])  # Number of correlations per a pol and a chan

    #
    # Loop over all the baselines
    #

    fout.write('# a1  tl1     elen1     tau1        a2  tl2     ' \
               'elen2     tau2         blen       w      tauw       acor\n') 
    for idat in range(ndat):
        a1 = ants1[idat]                  # First antenna 'name'
        a2 = ants2[idat]                  # Second antenna 'name'
        ia1 = idx_ant[a1]                 # First antenna index
        ia2 = idx_ant[a2]                 # Second antenna index

        #
        # Here we calculate the cable delays for each pair of tiles,
        # tau1 and tau2, from their electrical cable lengths, elen. 
        # Also, from the wave number (in meters), w, we get the
        # geometric delay between the tiles, tauw.
        # The delays are combined in tau:
        #
        tau1 = elen[ia1]/c_speed_of_light
        tau2 = elen[ia2]/c_speed_of_light

        # Let us not take this into account for a time being
        tauw = w[idat]/c_speed_of_light
        
        tau = tau2 - tau1 + tauw_sign*tauw

        acor = np.sinc(BW*tau)

        line = '%4d %4d %9.2f  %12.5e %4d %4d %9.2f  %12.5e ' \
               '%8.2f %8.2f %12.5e %9.6f\n' % \
               (a1, tile[ia1], elen[ia1], tau1,  \
                a2, tile[ia2], elen[ia2], tau2,  \
                blen[idat], w[idat], tauw, acor)

        fout.write(line)
        
        cdat[:,:,idat] = dat[:,:,idat]/acor
     
    return cdat

def decor(msname,metafits,n_tblk,single_time):	
	'''
	Function to correct the MS for amplitude decorrelation. MS convention is also converted from MWA coordinate (X=E->W and Y=N->S) to IAU coordinate (X=S->N, Y=W->E).
	
	Parameters
	----------
	msname : str 
		Name of the measurement set
	metafits : str
		Name of the metafits file of the observation
	n_tblk : int
		Number of time block for one round to avoid memory load in case of large dataset
	single_time : bool 
		Single time block or not
	Returns
	-------
	Decorrelation corrected measurement set in IAU convention
	'''
	code=vishead(vis=msname,mode='get',hdkey='fld_code')[0][0]
	code_list=code.split(',')
	if 'DECORRED' not in code_list:
		c_speed_of_light = 299792458. # m/s
		# Definging CASA6 casatools
		ms=mstools()
		md=msmetadata()
		qa=quanta()
		#
		# The sign of geometric delay term
		#
		tauw_sign = +1.
		#tauw_sign = -1.
			
		#
		# Read the cable delays (electrical cable lengths, elen)
		# from the metafits file.
		#
		ms_dir	=	msname
		md.open(msname)
		num_ant=md.nantennas()
		md.close()
		hl = fits.open(metafits)
		#
		# The second HDU is a fits table with 256 (nPol*nTile) rows.
		# The order of the rows is a little strange, but they can be
		# ordered in the same order as casa as follows:
		# Find the CASA order of the tiles and elen:
		#
		h1 = hl[1]
		tile_order = np.argsort(h1.data['Tile'])
		#
		# The lengths are prefixed by "EL_" which indicates that the
		# lengths given (in metres) are electrical lengths.
		# Cut off the 'EL_' prefixes and  reorder the e-lengths:
		#
		tile = [til for til in h1.data['Tile'][tile_order][::2]]
		tile = np.array(tile, dtype=int)
		elen = [float(Len[3:]) for Len in h1.data['Length'][tile_order][::2]]
		elen = np.array(elen, dtype=np.float64)
		hl.close()


		#raise  SystemExit

		sqrt_2 = np.sqrt(2.)
		r2d = 180./np.pi
		before = time.time()
		fout = open('decor_'+msname.split('/')[-1]+'.txt', 'w')

		ms.open(ms_dir, nomodify=False)

		#
		# Find exact number of individual times in MS
			#
		md = ms.metadata()

		if single_time:
			n_times = 1
		else:
			tims = md.timesforfield(0)
			n_times = len(tims)

		# Convert MS times from metadata to yyyymmdd hh:mm:ss 
		ymds = []
		for tim in tims:
			q = qa.quantity(tim, 's')
			ymd = qa.time(q, form='ymd', prec=8) # 8 means 12:34:56.75 format
			ymds.append(ymd)


		cwids = md.chanwidths(0)  
		cwid = cwids[0]
		IntgTime = md.exposuretime(scan=1,  spwid=0,  polid=0)
		tintg = IntgTime['value']     # Integration time. IntgTime['unit'] = 's' 
		an_samp = cwid*tintg       # = 20000, cplx samples in 0.5s of 40kHz channel
		n_samp = int(an_samp)
		ms.selectinit(datadescid=0)  # Untranslatable CASA dirty swearword


		if single_time:
			# Select the whole MS dataset
			rec =  ms.getdata(['data'])
			dat =  np.copy(rec['data'])
			ants1 = ms.getdata(['antenna1'])['antenna1']
			ants2 = ms.getdata(['antenna2'])['antenna2']
			tims =  ms.getdata(['time'])['time']
			uvw =   ms.getdata(['uvw'])['uvw']
			idx_ant = np.arange(0,num_ant,1) # Fill with 999999 where empty
			
			cdat = apply_acorr(dat, ants1, ants2, idx_ant, elen, tile, \
					       uvw, cwid, fout)

			rec['data'][:,:,:] = cdat
			print ('Time = %.2f' % (tims[0]))  # Only one single time
			code=vishead(vis=msname,mode='get',hdkey='fld_code')[0][0]
			code_list=code.split(',')	
			if 'IAU' not in code_list:
				rec_copy=copy.deepcopy(rec)
				# Converting from MWA coordinate ( X=EW and Y=NS ) to IAU coordinate (X=S->N, Y=W -> E)	.
				# Thus X=>-Y, and Y=>-X. So, XX=>YY,YY=>XX,XY=>YX,YX=>XY
				rec[0,:,:]=rec_copy[3,:,:]
				rec[1,:,:]=rec_copy[2,:,:]
				rec[2,:,:]=rec_copy[1,:,:]
				rec[3,:,:]=rec_copy[0,:,:]
				print ('Measurement set :'+msname+' is being converted in IAU convention')
				if len(code_list)==1 and code_list[0]=='':
					code+='IAU'
				else:
					code+=',IAU'
				vishead(vis=msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
				del rec_copy
			else:
				print ('Measurement set :'+msname+' is already in IAU convention')
			#
			# Put the corrected data back to the MS database
			#
			ms.putdata(rec)

			
		else:
			n_io = n_times//n_tblk
			if n_times % n_tblk != 0: n_io = n_io + 1
			tblk = np.arange(n_io + 1)*n_tblk
			tblk[-1] = n_times
			code=vishead(vis=msname,mode='get',hdkey='fld_code')[0][0]
			code_list=code.split(',')
			if 'IAU' not in code_list:
				do_iau_con=True
				print ('Measurement set :'+msname+' is being converted in IAU convention')
			else:
				print ('Measurement set :'+msname+' is already in IAU convention')
				do_iau_con=False
			for iblk in range(n_io):
				it0 = tblk[iblk]      # Block start (points at the first)
				it1 = tblk[iblk+1]    # Block end (points at the next after the last)
				tim0 = tims[it0]
				tim1 = tims[it1-1]
				ms.selectinit(reset=True)          # Reset the cumulative selection!!!
				ms.select({'time':[tim0, tim1]})
				rec =  ms.getdata(['data'])
				blkdat =  np.copy(rec['data'])
				blkants1 = ms.getdata(['antenna1'])['antenna1']
				blkants2 = ms.getdata(['antenna2'])['antenna2']
				blktims =  ms.getdata(['time'])['time']
				blkuvw =   ms.getdata(['uvw'])['uvw']
				blkcdat = np.zeros_like(blkdat)   # Corrected data to be written to MS

				#
				# Start/end indices of same times into blktims[]
				#
				ibt = np.where(np.diff(blktims) != 0)[0] + 1
				ibt = np.hstack((0, ibt, len(blkants1)))

				idx_ant = np.arange(0,num_ant,1) # Fill with 999999 where empty

				#raise  SystemExit
				#sys.exit(0)

				#
				# Apply amplitude corrections to the (large) piece of correlation data 
				# loaded to RAM. For each single time 'tim' and full set of 
				# baselines (ants1, ants2), a corresponding data block 'dat' is
				# extracted 
				# from 'blkdata', corrected, and stored in 'vdat'.
				#

				for itim in range(it0,it1):
					tim = tims[itim]
					ibt0 = ibt[itim-it0]
					ibt1 = ibt[itim-it0+1]
					dat = blkdat[:,:,ibt0:ibt1]
					ants1 = blkants1[ibt0:ibt1]
					ants2 = blkants2[ibt0:ibt1]
					uvw = blkuvw[:,ibt0:ibt1]

					cdat = apply_acorr(dat, ants1, ants2, idx_ant, elen, tile, \
							           uvw, cwid, fout)

					blkcdat[:,:,ibt0:ibt1] = cdat
					
					tymd = qa.time(qa.quantity(tim, 's'), form='ymd', prec=8)[0]

					print ('Time = %s' % (tymd))

				print ('Block %d-%d done.' % (it0,it1-1))
				if do_iau_con==True:
					blkcdat_copy=copy.deepcopy(blkcdat)

					# Converting from MWA coordinate ( X=EW and Y=NS ) to IAU coordinate (X=S->N, Y=W -> E)	.
					# Thus X=>-Y, and Y=>-X. So, XX=>YY,YY=>XX,XY=>YX,YX=>XY
					blkcdat[0,:,:]=blkcdat_copy[3,:,:]
					blkcdat[1,:,:]=blkcdat_copy[2,:,:]
					blkcdat[2,:,:]=blkcdat_copy[1,:,:]
					blkcdat[3,:,:]=blkcdat_copy[0,:,:]

					del blkcdat_copy
					#
					# Store the corrected correlations in the measurement set
					#
				rec['data'][:,:,:] = blkcdat
				
				
				# Put the corrected data back to the MS database
				#
				ms.putdata(rec)	
		if do_iau_con==True:
			if len(code_list)==1 and code_list[0]=='':
				code+='IAU'
			else:
				code+=',IAU'
			vishead(vis=msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
		ms.close()
		fout.close()
		if len(code_list)==1 and code_list[0]=='':
			code+='DECORRED'
		else:
			code+=',DECORRED'
		vishead(vis=msname,mode='put',hdkey='fld_code',hdvalue=np.array([code]))
	else:
		print ('Amplitude decorrelation is already corrected\n')
	return
