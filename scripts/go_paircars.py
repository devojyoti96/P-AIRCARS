from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from multiprocessing import Process
from paircars.basic_func import *
from paircars.access_ms import *
import os,copy,numpy as np,webbrowser,pickle,paircars,glob,time,getpass,tkinter as tk,subprocess,psutil
from PIL import Image,ImageTk
imagedir=os.path.abspath(os.path.dirname(paircars.__file__))
os.system('rm -rf casa*log')

class PAIRCARS_inputs:
	def __init__(self,root):
		# Master widget		
		self.root=root
		self.root.title('P-AIRCARS')
		self.root.geometry('1600x950')
		self.root.attributes('-alpha',0)
		self.image = Image.open(imagedir+"/mwasun.jpeg")
		self.img_copy= self.image.copy()
		self.background_image = ImageTk.PhotoImage(self.image)
		self.background = Label(self.root,image=self.background_image)
		self.background.pack(fill=BOTH, expand=YES)
		self.background.bind('<Configure>', self._resize_image)
		self.root.resizable(False, False) 
		# NCRA link
		self.image_left1=Image.open(imagedir+'/NCRA_logo.png')
		self.image_left1=self.image_left1.resize((120,120),Image.ANTIALIAS)
		self.left1=ImageTk.PhotoImage(self.image_left1)
		self.ncra_link=Button(self.root,command=self.open_ncra,image=self.left1)
		self.ncra_link.place(x=800,y=630)
		# MWA link
		self.image_left2=Image.open(imagedir+'/MWA_logo.jpeg')
		self.image_left2=self.image_left2.resize((120,120),Image.ANTIALIAS)
		self.left2=ImageTk.PhotoImage(self.image_left2)
		self.mwa_link=Button(self.root,command=self.open_mwa,image=self.left2)
		self.mwa_link.place(x=650,y=630)

		# Frame 1
		self.frame1=Frame(self.root,bg='white',highlightthickness=3)
		self.frame1.place(x=30,y=40,width=600,height=670)
		title1=Label(self.frame1,text='INPUTS',bg='white',fg='black',font=('times new roman',25)).place(relx=0.5,y=20,anchor=CENTER)

		# Data directory
		def delete_entry(event):
			if self.msdir_entry.get()=='Name of the directory of data.....':
				self.msdir_entry.delete(0, "end")
				self.msdir=''
				self.msdir_entry.config(fg='black')
		def restore_entry(event):
			if self.msdir_entry.get()=='':
				self.msdir_entry.delete(0, "end")
				self.msdir='Name of the directory of data.....'
				self.msdir_entry.insert(0,self.msdir)
				self.msdir_entry.config(fg='gray45')
		self.msdir='Name of the directory of data.....'
		msname=Label(self.frame1,text='Data Directory *',bg='white',fg='Black',font=('times new roman',15))
		msname.place(x=10,y=50)
		self.button=ttk.Button(self.frame1,text='browse',command=self.diropen1)
		self.button.place(x=500,y=50)
		self.msdir_entry=Entry(self.frame1,bg='lightgray',textvariable=self.msdir)
		self.msdir_entry.place(x=180,y=55,width=300)
		self.msdir_entry.insert(0,self.msdir)
		self.msdir_entry.config(fg='gray45')
		self.msdir_entry.bind("<FocusIn>",delete_entry)
		self.msdir_entry.bind("<FocusOut>",restore_entry)

		# Basedir
		def delete_entry(event):
			if self.basedir_entry.get()=='Name of the base directory .....':
				self.basedir_entry.delete(0, "end")
				self.basedir=''
				self.basedir_entry.config(fg='black')
		def restore_entry(event):
			if self.basedir_entry.get()=='':
				self.basedir_entry.delete(0, "end")
				self.basedir='Name of the base directory .....'
				self.basedir_entry.insert(0,self.msdir)
				self.basedir_entry.config(fg='gray45')
		self.basedir='Name of the base directory .....'
		basedir=Label(self.frame1,text='Base Directory *',bg='white',fg='Black',font=('times new roman',15))
		basedir.place(x=10,y=80)
		self.button=ttk.Button(self.frame1,text='browse',command=self.diropen2)
		self.button.place(x=500,y=80)
		self.basedir_entry=Entry(self.frame1,bg='lightgray',textvariable=self.basedir)
		self.basedir_entry.place(x=180,y=85,width=300)
		self.basedir_entry.insert(0,self.basedir)
		self.basedir_entry.config(fg='gray45')
		self.basedir_entry.bind("<FocusIn>",delete_entry)
		self.basedir_entry.bind("<FocusOut>",restore_entry)

		# Finalimagedir	
		def delete_entry(event):
			if self.fimagedir_entry.get()=='Name of the directory to store final images .....':
				self.fimagedir_entry.delete(0, "end")
				self.fimagedir=''
				self.fimagedir_entry.config(fg='black')
		def restore_entry(event):
			if self.fimagedir_entry.get()=='':
				self.fimagedir_entry.delete(0, "end")
				self.fimagedir='Name of the directory to store final images .....'
				self.fimagedir_entry.insert(0,self.msdir)
				self.fimagedir_entry.config(fg='gray45')
		self.fimagedir='Name of the directory to store final images .....'
		fimagedir=Label(self.frame1,text='Image Directory',bg='white',fg='Black',font=('times new roman',15))
		fimagedir.place(x=10,y=110)
		self.button=ttk.Button(self.frame1,text='browse',command=self.diropen3)
		self.button.place(x=500,y=110)
		self.fimagedir_entry=Entry(self.frame1,bg='lightgray',textvariable=self.fimagedir)
		self.fimagedir_entry.place(x=180,y=115,width=300)
		self.fimagedir_entry.insert(0,self.fimagedir)
		self.fimagedir_entry.config(fg='gray45')
		self.fimagedir_entry.bind("<FocusIn>",delete_entry)
		self.fimagedir_entry.bind("<FocusOut>",restore_entry)
		
		# Timerange
		def delete_entry(event):
			if self.tmrange_entry.get()=='yy0/mm0/dd0/hh0:mm0:ss0.ff0~yy1/mm1/dd1/hh1:mm1:ss1.ff1,yy2/mm2/dd2/hh2:mm2:ss2.ff2~....':
				self.tmrange_entry.delete(0, "end")
				self.timerange=''
				self.tmrange_entry.config(fg='black')
		def restore_entry(event):
			if self.tmrange_entry.get()=='':
				self.tmrange_entry.delete(0, "end")
				self.timerange='yy0/mm0/dd0/hh0:mm0:ss0.ff0~yy1/mm1/dd1/hh1:mm1:ss1.ff1,yy2/mm2/dd2/hh2:mm2:ss2.ff2~....'
				self.tmrange_entry.insert(0,self.timerange)
				self.tmrange_entry.config(fg='gray45')
		self.timerange='yy0/mm0/dd0/hh0:mm0:ss0.ff0~yy1/mm1/dd1/hh1:mm1:ss1.ff1,yy2/mm2/dd2/hh2:mm2:ss2.ff2~....'
		timerange=Label(self.frame1,text='Time range',bg='white',fg='Black',font=('times new roman',15))
		timerange.place(x=10,y=140)
		self.tmrange_entry=Entry(self.frame1,bg='lightgray',textvariable=self.timerange)
		self.tmrange_entry.place(x=180,y=145,width=400)
		self.tmrange_entry.insert(0,self.timerange)
		self.tmrange_entry.config(fg='gray45')
		self.tmrange_entry.bind("<FocusIn>",delete_entry)
		self.tmrange_entry.bind("<FocusOut>",restore_entry)

		# Chanrange
		def delete_entry(event):
			if self.chanrange_entry.get()=='ch0~ch1,ch2~ch3,....':
				self.chanrange_entry.delete(0, "end")
				self.chanrange=''
				self.chanrange_entry.config(fg='black')
		def restore_entry(event):
			if self.chanrange_entry.get()=='':
				self.chanrange_entry.delete(0, "end")
				self.chanrange='ch0~ch1,ch2~ch3,....'
				self.chanrange_entry.insert(0,self.chanrange)
				self.chanrange_entry.config(fg='gray45')
		self.chanrange='ch0~ch1,ch2~ch3,....'
		chanrange=Label(self.frame1,text='Channel range',bg='white',fg='Black',font=('times new roman',15))
		chanrange.place(x=10,y=170)
		self.chanrange_entry=Entry(self.frame1,bg='lightgray',textvariable=self.chanrange)
		self.chanrange_entry.place(x=180,y=175,width=400)
		self.chanrange_entry.insert(0,self.chanrange)
		self.chanrange_entry.config(fg='gray45')
		self.chanrange_entry.bind("<FocusIn>",delete_entry)
		self.chanrange_entry.bind("<FocusOut>",restore_entry)

		# Caltable
		def delete_entry(event):
			if self.cal_entry.get()=='/path/to/caltable0,/path/to/caltable1,....':
				self.cal_entry.delete(0, "end")
				self.caltable=''
				self.cal_entry.config(fg='black')
		def restore_entry(event):
			if self.cal_entry.get()=='':
				self.cal_entry.delete(0, "end")
				self.caltable='/path/to/caltable0,/path/to/caltable1,....'
				self.cal_entry.insert(0,self.caltable)
				self.cal_entry.config(fg='gray45')
		self.caltable='/path/to/caltable0,/path/to/caltable1,....'
		cal=Label(self.frame1,text='Caltable Names',bg='white',fg='Black',font=('times new roman',15))
		cal.place(x=10,y=200)
		self.cal_entry=Entry(self.frame1,bg='lightgray',textvariable=self.caltable)
		self.cal_entry.place(x=180,y=205,width=400)
		self.cal_entry.insert(0,self.caltable)
		self.cal_entry.config(fg='gray45')
		self.cal_entry.bind("<FocusIn>",delete_entry)
		self.cal_entry.bind("<FocusOut>",restore_entry)

		# Robustness factor
		self.safety=IntVar()
		self.safety.set(1)
		safety=Label(self.frame1,text='Robustness',bg='white',fg='Black',font=('times new roman',15))
		safety.place(x=10,y=230)
		Radiobutton(self.frame1,variable=self.safety,value=0,text='0',bg='white',fg='gray30',highlightbackground = "white").place(x=110,y=235)
		Radiobutton(self.frame1,variable=self.safety,value=1,text='1',bg='white',fg='gray30',highlightbackground = "white").place(x=150,y=235)
		Radiobutton(self.frame1,variable=self.safety,value=2,text='2',bg='white',fg='gray30',highlightbackground = "white").place(x=190,y=235)

		# Quality factor
		self.quality=IntVar()
		self.quality.set(1)
		quality=Label(self.frame1,text='Quality',bg='white',fg='Black',font=('times new roman',15))
		quality.place(x=245,y=230)
		Radiobutton(self.frame1,variable=self.quality,value=0,text='0',bg='white',fg='gray30',highlightbackground = "white").place(x=315,y=235)
		Radiobutton(self.frame1,variable=self.quality,value=1,text='1',bg='white',fg='gray30',highlightbackground = "white").place(x=355,y=235)
		Radiobutton(self.frame1,variable=self.quality,value=2,text='2',bg='white',fg='gray30',highlightbackground = "white").place(x=395,y=235)

		# Verbose
		verbose=Label(self.frame1,text='Verbose',bg='white',fg='Black',font=('times new roman',15))
		verbose.place(x=465,y=230)
		self.verbose=BooleanVar()
		self.verbose.set(False)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.verbose,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=545,y=230)

		# Keep logger
		logger=Label(self.frame1,text='Save Logs',bg='white',fg='Black',font=('times new roman',15))
		logger.place(x=10,y=260)
		self.logger=BooleanVar()
		self.logger.set(False)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.logger,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=105,y=260)
		
		# Interactive
		interactive=Label(self.frame1,text='Interactive',bg='white',fg='Black',font=('times new roman',15))
		interactive.place(x=165,y=260)
		self.interactive=BooleanVar()
		self.interactive.set(False)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.interactive,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=265,y=260)

		# Decor
		dodecor=Label(self.frame1,text='Decorreleation Correction',bg='white',fg='Black',font=('times new roman',15))
		dodecor.place(x=315,y=260)
		self.dodecor=BooleanVar()
		self.dodecor.set(True)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.dodecor,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=545,y=260)

		# Refant
		refants=[str(i) for i in range(0,30)]
		self.refant=StringVar()
		self.refant.set(refants[1])
		refant=Label(self.frame1,text='Reference Antenna',bg='white',fg='Black',font=('times new roman',15))
		refant.place(x=10,y=290)
		self.refant_options=ttk.Combobox(self.frame1,textvariable=self.refant,values=refants,width=3,state='readonly')
		self.refant_options.place(x=185,y=295)

		# Calc image and selfcal parameters
		autocal=Label(self.frame1,text='Auto-calculate Parameters',bg='white',fg='Black',font=('times new roman',15))
		autocal.place(x=315,y=290)
		self.autocal=BooleanVar()
		self.autocal.set(True)
		def onauto():
			if self.autocal.get()==True:
				for child in self.frame2.winfo_children():
					child.configure(state='disable')	
			else:
				for child in self.frame2.winfo_children():
					child.configure(state='normal')
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.autocal,onvalue=True,offvalue=False,command=onauto,\
				highlightbackground = "white")
		c.place(x=545,y=290)

		self.use_wsclean=True
		# Maskfile
		def delete_entry0(event):
			if self.mask_entry.get()=='CASA mask path.....':
				self.mask_entry.delete(0, "end")
				self.maskfile=''
				self.mask_entry.config(fg='black')
		def restore_entry0(event):
			if self.mask_entry.get()=='':
				self.mask_entry.delete(0, "end")
				self.maskfile='CASA mask path.....'
				self.mask_entry.insert(0,self.maskfile)
				self.mask_entry.config(fg='gray45')
		self.maskfile=StringVar()
		self.maskfile='CASA mask path.....'
		maskfile=Label(self.frame1,text='CASA Mask',bg='white',fg='Black',font=('times new roman',15))
		maskfile.place(x=10,y=350)
		self.button1=ttk.Button(self.frame1,text='browse',command=self.diropen4)
		self.button1.place(x=500,y=350)
		self.mask_entry=Entry(self.frame1,bg='lightgray',textvariable=self.maskfile)
		self.mask_entry.place(x=125,y=355,width=360)
		self.mask_entry.insert(0,self.maskfile)
		self.mask_entry.config(fg='gray45')
		self.mask_entry.bind("<FocusIn>",delete_entry0)
		self.mask_entry.bind("<FocusOut>",restore_entry0)
		self.mask_entry.config(state=DISABLED)
		self.button1.config(state=DISABLED)
		self.mask_entry.config(state=DISABLED)

		# Maskfile
		def delete_entry1(event):
			if self.maskstr_entry.get()=='Mask string in CASA format.....':
				self.maskstr_entry.delete(0, "end")
				self.maskstr=''
				self.maskstr_entry.config(fg='black')
		def restore_entry1(event):
			if self.maskstr_entry.get()=='':
				self.maskstr_entry.delete(0, "end")
				self.maskstr='Mask string in CASA format.....'
				self.maskstr_entry.insert(0,self.maskstr)
				self.maskstr_entry.config(fg='gray45')
		self.maskstr=StringVar()
		self.maskstr='Mask string in CASA format.....'
		maskstr=Label(self.frame1,text='CASA Mask String',bg='white',fg='Black',font=('times new roman',15))
		maskstr.place(x=10,y=380)
		self.maskstr_entry=Entry(self.frame1,bg='lightgray',textvariable=self.maskstr)
		self.maskstr_entry.place(x=190,y=385,width=395)
		self.maskstr_entry.insert(0,self.maskstr)
		self.maskstr_entry.config(fg='gray45')
		self.maskstr_entry.bind("<FocusIn>",delete_entry1)
		self.maskstr_entry.bind("<FocusOut>",restore_entry1)
		self.maskstr_entry.config(state=DISABLED)
		self.button1.config(state=DISABLED)
		self.maskstr_entry.config(state=DISABLED)

		# Automask
		automask=Label(self.frame1,text='CASA Auto-masking',bg='white',fg='Black',font=('times new roman',15))
		automask.place(x=190,y=320)
		self.automask=BooleanVar()
		self.automask.set(False)
		self.c1=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.automask,onvalue=True,offvalue=False,highlightbackground = "white")
		self.c1.place(x=375,y=320)
		self.c1.config(state=DISABLED)
		
		# Use WSClean
		def onClicked():
			self.use_wsclean=self.wsclean_input.get()	
			if self.use_wsclean==True:
				self.maskstr_entry.delete(0, "end")
				self.maskstr='Mask string in CASA format.....'
				self.maskstr_entry.insert(0,self.maskstr)
				self.maskstr_entry.config(fg='gray45')
				self.mask_entry.delete(0, "end")
				self.maskfile='CASA mask path.....'
				self.mask_entry.insert(0,self.maskfile)
				self.mask_entry.config(fg='gray45')
				self.c1.deselect()
				self.mask_entry.config(state=DISABLED)
				self.button1.config(state=DISABLED)
				self.maskstr_entry.config(state=DISABLED)
				self.c1.config(state=DISABLED)
			else:
				self.mask_entry.config(state=NORMAL)
				self.button1.config(state=NORMAL)
				self.maskstr_entry.config(state=NORMAL)	
				self.c1.config(state=NORMAL)
		wsclean=Label(self.frame1,text='Use WSClean',bg='white',fg='Black',font=('times new roman',15))
		wsclean.place(x=10,y=320)
		self.wsclean_input=BooleanVar()
		self.wsclean_input.set(True)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.wsclean_input,onvalue=True,offvalue=False,command=onClicked,\
				highlightbackground = "white")
		c.place(x=145,y=320)

		# email
		def delete_entry1(event):
			if self.email_entry.get()=='Enter email address to send notifications.....':
				self.email_entry.delete(0, "end")
				self.email=''
				self.email_entry.config(fg='black')
		def restore_entry1(event):
			if self.email_entry.get()=='':
				self.email_entry.delete(0, "end")
				self.email='Enter email address to send notifications.....'
				self.email_entry.insert(0,self.email)
				self.email_entry.config(fg='gray45')
		self.email=StringVar()
		self.email='Enter email address to send notifications.....'
		email=Label(self.frame1,text='e-mail',bg='white',fg='Black',font=('times new roman',15))
		email.place(x=10,y=410)
		self.email_entry=Entry(self.frame1,bg='lightgray',textvariable=self.email)
		self.email_entry.place(x=80,y=415,width=275)
		self.email_entry.insert(0,self.email)
		self.email_entry.config(fg='gray45')
		self.email_entry.bind("<FocusIn>",delete_entry1)
		self.email_entry.bind("<FocusOut>",restore_entry1)

		# Notification
		self.send_notification=True
		def onClicked():
			self.send_notification=self.notification_input.get()	
			if self.send_notification==False:
				self.email_entry.delete(0, "end")
				self.email='Enter email address to send notifications.....'
				self.email_entry.insert(0,self.email)
				self.email_entry.config(fg='gray45')
				self.email_entry.config(state=DISABLED)
			else:
				self.email_entry.config(state=NORMAL)
		notification=Label(self.frame1,text='Notification',bg='white',fg='Black',font=('times new roman',15))
		notification.place(x=430,y=320)
		self.notification_input=BooleanVar()
		self.notification_input.set(True)
		self.c2=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.notification_input,onvalue=True,offvalue=False,command=onClicked\
			,highlightbackground = "white")
		self.c2.place(x=545,y=320)

		# do bandpass
		bandpass=Label(self.frame1,text='Bandpass',bg='white',fg='Black',font=('times new roman',15))
		bandpass.place(x=360,y=410)
		self.bandpass=BooleanVar()
		self.bandpass.set(True)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.bandpass,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=445,y=410)

		# do polcal
		polcal=Label(self.frame1,text='Polcal',bg='white',fg='Black',font=('times new roman',15))
		polcal.place(x=485,y=410)
		self.polcal=BooleanVar()
		self.polcal.set(True)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.polcal,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=545,y=410)

		# Image frequency interval
		freqintervals=np.linspace(10,30000,int(30000/10))
		freqintervals=["{:0.1f}".format(i) for i in freqintervals.tolist()]
		self.freqint=DoubleVar()
		self.freqint.set(freqintervals[15])
		freqint=Label(self.frame1,text='Frequency Interval (kHz)',bg='white',fg='Black',font=('times new roman',15))
		freqint.place(x=10,y=440)
		self.freqint_options=ttk.Combobox(self.frame1,textvariable=self.freqint,values=freqintervals,width=6)
		self.freqint_options.place(x=230,y=445)

		# Image time interval
		timeintervals=np.linspace(0.5,240,480)
		timeintervals=timeintervals.tolist()
		self.timeint=DoubleVar()
		self.timeint.set(timeintervals[0])
		timeint=Label(self.frame1,text='Time Interval (s)',bg='white',fg='Black',font=('times new roman',15))
		timeint.place(x=360,y=440)
		self.timeint_options=ttk.Combobox(self.frame1,textvariable=self.timeint,values=timeintervals,width=6)
		self.timeint_options.place(x=515,y=445)

		# Image frequency width
		freqwids=np.linspace(10,30000,int(30000/10))
		freqwids=["{:0.1f}".format(i) for i in freqwids.tolist()]
		self.freqwid=DoubleVar()
		self.freqwid.set(freqwids[15])
		freqwid=Label(self.frame1,text='Frequency Width (kHz)',bg='white',fg='Black',font=('times new roman',15))
		freqwid.place(x=10,y=470)
		self.freqwid_options=ttk.Combobox(self.frame1,textvariable=self.freqwid,values=freqwids,width=6)
		self.freqwid_options.place(x=230,y=475)

		# Image time width
		timewids=np.linspace(0.5,10,20)
		timewids=timewids.tolist()
		self.timewid=DoubleVar()
		self.timewid.set(timewids[0])
		timeint=Label(self.frame1,text='Temporal Width (s)',bg='white',fg='Black',font=('times new roman',15))
		timeint.place(x=335,y=470)
		self.timewid_options=ttk.Combobox(self.frame1,textvariable=self.timewid,values=timewids,width=6)
		self.timewid_options.place(x=515,y=475)

		# cpu_frac
		cpufrac=Label(self.frame1,text='CPU (%)',bg='white',fg='Black',font=('times new roman',15))
		cpufrac.place(x=10,y=510)
		self.cpufrac_options=Scale(self.frame1,from_=0.0,to=100.0,orient=HORIZONTAL,resolution=0.5,length=200,bg='white',highlightbackground = "white")
		self.cpufrac_options.place(x=100,y=500)
		self.cpufrac_options.set(50)

		# clear screen
		clearscreen=Label(self.frame1,text='Clear Virtual Screens',bg='white',fg='Black',font=('times new roman',15))
		clearscreen.place(x=355,y=510)
		self.clearscreen=BooleanVar()
		self.clearscreen.set(True)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.clearscreen,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=545,y=510)

		# savemodel
		savemodel=Label(self.frame1,text='Save models',bg='white',fg='Black',font=('times new roman',15))
		savemodel.place(x=10,y=540)
		self.savemodel=BooleanVar()
		self.savemodel.set(False)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.savemodel,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=115,y=540)

		# saveresiduals
		saveresiduals=Label(self.frame1,text='Save residuals',bg='white',fg='Black',font=('times new roman',15))
		saveresiduals.place(x=185,y=540)
		self.saveresiduals=BooleanVar()
		self.saveresiduals.set(False)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.saveresiduals,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=305,y=540)

		# x cutout
		def delete_entry(event):
			if self.xcut_entry.get()==3:
				self.xcut_entry.delete(0, "end")
				self.xcut=''
				self.xcut_entry.config(fg='black')
		def restore_entry(event):
			if self.xcut_entry.get()=='':
				self.xcut_entry.delete(0, "end")
				self.xcut=3
				self.xcut_entry.insert(0,self.xcut)
				self.xcut_entry.config(fg='gray45')
		self.xcut=DoubleVar()
		self.xcut=3
		xcut=Label(self.frame1,text='XY Cutout (degree)',bg='white',fg='Black',font=('times new roman',15))
		xcut.place(x=375,y=540)
		self.xcut_entry=Entry(self.frame1,bg='lightgray',textvariable=self.xcut)
		self.xcut_entry.place(x=545,y=545,width=30)
		self.xcut_entry.insert(0,self.xcut)
		self.xcut_entry.config(fg='gray45')
		self.xcut_entry.bind("<FocusIn>",delete_entry)
		self.xcut_entry.bind("<FocusOut>",restore_entry)

		# Flag
		flag=Label(self.frame1,text='Perform Flagging',bg='white',fg='Black',font=('times new roman',15))
		flag.place(x=10,y=570)
		self.flag=BooleanVar()
		self.flag.set(True)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.flag,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=160,y=570)

		# aNKFlag
		ankflag=Label(self.frame1,text='Use aNKflag',bg='white',fg='Black',font=('times new roman',15))
		ankflag.place(x=220,y=570)
		self.ankflag=BooleanVar()
		self.ankflag.set(True)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.ankflag,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=330,y=570)
	
		# HPC
		ishpc=Label(self.frame1,text='HPC environment',bg='white',fg='Black',font=('times new roman',15))
		ishpc.place(x=380,y=570)
		self.ishpc=BooleanVar()
		self.ishpc.set(False)	
		def onhpc():
			if self.ishpc.get()==True:
				for child in self.frame3.winfo_children():
					child.configure(state='normal')	
			elif self.ishpc.get()==False:
				for child in self.frame3.winfo_children():
					child.configure(state='disabled')
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.ishpc,command=onhpc,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=540,y=570)

		# Load input file
		def delete_entry(event):
			if self.loadinput_entry.get()=='Load P-AIRCARS input file.....':
				self.loadinput_entry.delete(0, "end")
				self.loadinput=''
				self.loadinput_entry.config(fg='black')
		def restore_entry(event):
			if self.loadinput_entry.get()=='':
				self.loadinput_entry.delete(0, "end")
				self.loadinput='Load P-AIRCARS input file.....'
				self.loadinput_entry.insert(0,self.loadinput)
				self.loadinput_entry.config(fg='gray45')
		self.loadinput='Load P-AIRCARS input file.....'
		loadinput=Label(self.frame1,text='Input file',bg='white',fg='Black',font=('times new roman',15))
		loadinput.place(x=10,y=600)
		self.button=ttk.Button(self.frame1,text='load',command=self.fileopen1)
		self.button.place(x=490,y=600)
		self.loadinput_entry=Entry(self.frame1,bg='lightgray',textvariable=self.loadinput)
		self.loadinput_entry.place(x=110,y=605,width=360)
		self.loadinput_entry.insert(0,self.loadinput)
		self.loadinput_entry.config(fg='gray45')
		self.loadinput_entry.bind("<FocusIn>",delete_entry)
		self.loadinput_entry.bind("<FocusOut>",restore_entry)


		self.fresh=BooleanVar()
		self.fresh.set(True)	
		self.restart=BooleanVar()
		self.restart.set(False)
		# Restart
		def onrestart():
			if self.restart.get()==True:
				self.fresh.set(False)
		def onfresh():
			if self.fresh.get()==True:
				self.restart.set(False)
		restart=Label(self.frame1,text='Restart P-AIRCARS',bg='white',fg='Black',font=('times new roman',15))
		restart.place(x=10,y=630)	
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.restart,command=onrestart,onvalue=True,offvalue=False,\
				highlightbackground = "white")
		c.place(x=180,y=630)

		# fresh start
		fresh=Label(self.frame1,text='Fresh Start P-AIRCARS',bg='white',fg='Black',font=('times new roman',15))
		fresh.place(x=330,y=630)
		c=Checkbutton(self.frame1,text='',bg='white',fg='Black',font=('times new roman',15),variable=self.fresh,command=onfresh,onvalue=True,offvalue=False,\
				highlightbackground = "white")
		c.place(x=540,y=630)

		# Open log
		self.log_button=Button(self.root,text='Open Log',command=self.logview,font=('times new roman',15))
		self.log_button.place(x=660,y=770)

		# Download data
		self.download_button=Button(self.root,text='Download Data',command=self.download_data,font=('times new roman',15))
		self.download_button.place(x=785,y=770)

		# Frame 2
		self.frame2=Frame(self.root,bg='white',highlightthickness=3)
		self.frame2.place(x=970,y=40,width=600,height=450)
		title2=Label(self.frame2,text='ADVANCED INPUTS',bg='white',fg='black',font=('times new roman',25)).place(relx=0.5,y=20,anchor=CENTER)
		
		# Cell size
		def delete_entry(event):
			if self.cellsize_entry.get()=='xxarcsec or xxarcmin':
				self.cellsize_entry.delete(0, "end")
				self.cellsize=''
				self.cellsize_entry.config(fg='black')
		def restore_entry(event):
			if self.cellsize_entry.get()=='':
				self.cellsize_entry.delete(0, "end")
				self.cellsize='xxarcsec or xxarcmin'
				self.cellsize_entry.insert(0,self.cellsize)
				self.cellsize_entry.config(fg='gray45')
		self.cellsize=StringVar()
		self.cellsize='xxarcsec or xxarcmin'
		cellsize=Label(self.frame2,text='Pixelsize',bg='white',fg='Black',font=('times new roman',15))
		cellsize.place(x=10,y=50)
		self.cellsize_entry=Entry(self.frame2,bg='lightgray',textvariable=self.cellsize)
		self.cellsize_entry.place(x=95,y=55,width=150)
		self.cellsize_entry.insert(0,self.cellsize)
		self.cellsize_entry.config(fg='gray45')
		self.cellsize_entry.bind("<FocusIn>",delete_entry)
		self.cellsize_entry.bind("<FocusOut>",restore_entry)

		# Image size
		def delete_entry(event):
			if self.imsize_entry.get()=='Number of pixels':
				self.imsize_entry.delete(0, "end")
				self.imsize=''
				self.imsize_entry.config(fg='black')
		def restore_entry(event):
			if self.imsize_entry.get()=='':
				self.imsize_entry.delete(0, "end")
				self.imsize='Number of pixels'
				self.imsize_entry.insert(0,self.imsize)
				self.imsize_entry.config(fg='gray45')
		self.imsize=StringVar()
		self.imsize='Number of pixels'
		imsize=Label(self.frame2,text='Number of pixels',bg='white',fg='Black',font=('times new roman',15))
		imsize.place(x=270,y=50)
		self.imsize_entry=Entry(self.frame2,bg='lightgray',textvariable=self.imsize)
		self.imsize_entry.place(x=420,y=53,width=150)
		self.imsize_entry.insert(0,self.imsize)
		self.imsize_entry.config(fg='gray45')
		self.imsize_entry.bind("<FocusIn>",delete_entry)
		self.imsize_entry.bind("<FocusOut>",restore_entry)
			
		# Multiscale scales
		def delete_entry(event):
			if self.scales_entry.get()=='0,3,9,6':
				self.scales_entry.delete(0, "end")
				self.scales=''
				self.scales_entry.config(fg='black')
		def restore_entry(event):
			if self.scales_entry.get()=='':
				self.scales_entry.delete(0, "end")
				self.scales='0,3,9,6'
				self.scales_entry.insert(0,self.scales)
				self.scales_entry.config(fg='gray45')
		self.scales=StringVar()
		self.scales='0,3,9,6'
		scales=Label(self.frame2,text='Multiscale scales',bg='white',fg='Black',font=('times new roman',15))
		scales.place(x=10,y=90)
		self.scales_entry=Entry(self.frame2,bg='lightgray',textvariable=self.scales)
		self.scales_entry.place(x=160,y=93,width=150)
		self.scales_entry.insert(0,self.scales)
		self.scales_entry.config(fg='gray45')
		self.scales_entry.bind("<FocusIn>",delete_entry)
		self.scales_entry.bind("<FocusOut>",restore_entry)

		# UVtaper
		def delete_entry(event):
			if self.uvtaper_entry.get()=='xxlambda or xxklambda':
				self.uvtaper_entry.delete(0, "end")
				self.uvtaper=''
				self.uvtaper_entry.config(fg='black')
		def restore_entry(event):
			if self.uvtaper_entry.get()=='':
				self.uvtaper_entry.delete(0, "end")
				self.uvtaper='xxlambda or xxklambda'
				self.uvtaper_entry.insert(0,self.uvtaper)
				self.uvtaper_entry.config(fg='gray45')
		self.uvtaper=StringVar()
		self.uvtaper='xxlambda or xxklambda'
		uvtaper=Label(self.frame2,text='UV-taper',bg='white',fg='Black',font=('times new roman',15))
		uvtaper.place(x=320,y=90)
		self.uvtaper_entry=Entry(self.frame2,bg='lightgray',textvariable=self.uvtaper)
		self.uvtaper_entry.place(x=405,y=93,width=165)
		self.uvtaper_entry.insert(0,self.uvtaper)
		self.uvtaper_entry.config(fg='gray45')
		self.uvtaper_entry.bind("<FocusIn>",delete_entry)
		self.uvtaper_entry.bind("<FocusOut>",restore_entry)

		# Start sigma
		sigmas=np.arange(7,12.25,0.25).tolist()
		self.sigma=DoubleVar()
		self.sigma.set(sigmas[12])
		sigma=Label(self.frame2,text='Start Sigma',bg='white',fg='Black',font=('times new roman',15))
		sigma.place(x=10,y=130)
		self.sigma_options=ttk.Combobox(self.frame2,textvariable=self.sigma,values=sigmas,width=5,state='readonly')
		self.sigma_options.place(x=125,y=135)

		# Sigma step
		steps=["{:0.1f}".format(i) for i in np.arange(0.1,2.1,0.1).tolist()]
		self.step=DoubleVar()
		self.step.set(steps[4])
		step=Label(self.frame2,text='Sigma Step',bg='white',fg='Black',font=('times new roman',15))
		step.place(x=200,y=130)
		self.step_options=ttk.Combobox(self.frame2,textvariable=self.step,values=steps,width=3,state='readonly')
		self.step_options.place(x=310,y=135)


		# Minimum sigma
		minsigmas=["{:0.2f}".format(i) for i in np.arange(5,11,0.25).tolist()]
		self.minsigma=DoubleVar()
		self.minsigma.set(minsigmas[12])
		minsigma=Label(self.frame2,text='Minimum Sigma',bg='white',fg='Black',font=('times new roman',15))
		minsigma.place(x=365,y=130)
		self.minsigma_options=ttk.Combobox(self.frame2,textvariable=self.minsigma,values=minsigmas,width=5,state='readonly')
		self.minsigma_options.place(x=515,y=135)


		# Res flux
		resfracs=["{:0.2f}".format(i) for i in np.arange(0.05,0.55,0.05).tolist()]
		self.resfrac=DoubleVar()
		self.resfrac.set(resfracs[2])
		resfrac=Label(self.frame2,text='Residual Flux Fraction',bg='white',fg='Black',font=('times new roman',15))
		resfrac.place(x=10,y=170)
		self.resfrac_options=ttk.Combobox(self.frame2,textvariable=self.resfrac,values=resfracs,width=4,state='readonly')
		self.resfrac_options.place(x=215,y=175)
	
		# UVrange
		def delete_entry(event):
			if self.uvrange_entry.get()=='CASA uvrange format':
				self.uvrange_entry.delete(0, "end")
				self.uvrange=''
				self.uvrange_entry.config(fg='black')
		def restore_entry(event):
			if self.uvrange_entry.get()=='':
				self.uvrange_entry.delete(0, "end")
				self.uvrange='CASA uvrange format'
				self.uvrange_entry.insert(0,self.uvrange)
				self.uvrange_entry.config(fg='gray45')
		self.uvrange=StringVar()
		self.uvrange='CASA uvrange format'
		uvrange=Label(self.frame2,text='UV range',bg='white',fg='Black',font=('times new roman',15))
		uvrange.place(x=315,y=175)
		self.uvrange_entry=Entry(self.frame2,bg='lightgray',textvariable=self.uvrange)
		self.uvrange_entry.place(x=405,y=175,width=165)
		self.uvrange_entry.insert(0,self.uvrange)
		self.uvrange_entry.config(fg='gray45')
		self.uvrange_entry.bind("<FocusIn>",delete_entry)
		self.uvrange_entry.bind("<FocusOut>",restore_entry)

		# Skip frequency
		freqintervals=np.linspace(10,2000,int(2000/10))
		freqintervals=["{:0.1f}".format(i) for i in freqintervals.tolist()]
		self.skipfreq=DoubleVar()
		self.skipfreq.set(freqintervals[127])
		skipfreq=Label(self.frame2,text='Skip Frequency (kHz)',bg='white',fg='Black',font=('times new roman',15))
		skipfreq.place(x=10,y=215)
		self.skipfreq_options=ttk.Combobox(self.frame2,textvariable=self.skipfreq,values=freqintervals,width=6)
		self.skipfreq_options.place(x=210,y=220)

		# Skip time
		timeintervals=np.linspace(0.5,240,480)
		timeintervals=timeintervals.tolist()
		self.skiptime=DoubleVar()
		self.skiptime.set(timeintervals[239])
		skiptime=Label(self.frame2,text='Skip Time (s)',bg='white',fg='Black',font=('times new roman',15))
		skiptime.place(x=370,y=215)
		self.skiptime_options=ttk.Combobox(self.frame2,textvariable=self.skiptime,values=timeintervals,width=6)
		self.skiptime_options.place(x=505,y=220)

		# Gain minsnr
		minsnrs=["{:0.1f}".format(i) for i in np.arange(3,5.5,0.5).tolist()]
		self.minsnr=DoubleVar()
		self.minsnr.set(minsnrs[2])
		minsnr=Label(self.frame2,text='Gain Minimum SNR',bg='white',fg='Black',font=('times new roman',15))
		minsnr.place(x=10,y=255)
		self.minsnr_options=ttk.Combobox(self.frame2,textvariable=self.minsnr,values=minsnrs,width=4)
		self.minsnr_options.place(x=195,y=260)

		# DRdeltarms
		def delete_entry(event):
			if self.drrms_entry.get()=='DR rms step':
				self.drrms_entry.delete(0, "end")
				self.drrms=''
				self.drrms_entry.config(fg='black')
		def restore_entry(event):
			if self.drrms_entry.get()=='':
				self.drrms_entry.delete(0, "end")
				self.drrms='DR rms step'
				self.drrms_entry.insert(0,self.drrms)
				self.drrms_entry.config(fg='gray45')
		self.drrms=DoubleVar()
		self.drrms='DR rms step'
		drrms=Label(self.frame2,text='DR RMS step',bg='white',fg='Black',font=('times new roman',15))
		drrms.place(x=280,y=255)
		self.drrms_entry=Entry(self.frame2,bg='lightgray',textvariable=self.drrms)
		self.drrms_entry.place(x=405,y=260,width=165)
		self.drrms_entry.insert(0,self.drrms)
		self.drrms_entry.config(fg='gray45')
		self.drrms_entry.bind("<FocusIn>",delete_entry)
		self.drrms_entry.bind("<FocusOut>",restore_entry)

		# DRdeltaneg
		def delete_entry(event):
			if self.drneg_entry.get()=='DR negative step':
				self.drneg_entry.delete(0, "end")
				self.drneg=''
				self.drneg_entry.config(fg='black')
		def restore_entry(event):
			if self.drneg_entry.get()=='':
				self.drneg_entry.delete(0, "end")
				self.drrms='DR negative step'
				self.drneg_entry.insert(0,self.drneg)
				self.drneg_entry.config(fg='gray45')
		self.drneg=DoubleVar()
		self.drneg='DR negative step'
		drneg=Label(self.frame2,text='DR Negative step',bg='white',fg='Black',font=('times new roman',15))
		drneg.place(x=10,y=290)
		self.drneg_entry=Entry(self.frame2,bg='lightgray',textvariable=self.drneg)
		self.drneg_entry.place(x=175,y=295,width=125)
		self.drneg_entry.insert(0,self.drneg)
		self.drneg_entry.config(fg='gray45')
		self.drneg_entry.bind("<FocusIn>",delete_entry)
		self.drneg_entry.bind("<FocusOut>",restore_entry)

		# minDR
		def delete_entry(event):
			if self.mindr_entry.get()=='Minimum DR':
				self.mindr_entry.delete(0, "end")
				self.mindr=''
				self.mindr_entry.config(fg='black')
		def restore_entry(event):
			if self.mindr_entry.get()=='':
				self.mindr_entry.delete(0, "end")
				self.mindr='Minimum DR'
				self.mindr_entry.insert(0,self.mindr)
				self.mindr_entry.config(fg='gray45')
		self.mindr=IntVar()
		self.mindr='Minimum DR'
		mindr=Label(self.frame2,text='Minimum DR',bg='white',fg='Black',font=('times new roman',15))
		mindr.place(x=320,y=290)
		self.mindr_entry=Entry(self.frame2,bg='lightgray',textvariable=self.mindr)
		self.mindr_entry.place(x=455,y=295,width=115)
		self.mindr_entry.insert(0,self.mindr)
		self.mindr_entry.config(fg='gray45')
		self.mindr_entry.bind("<FocusIn>",delete_entry)
		self.mindr_entry.bind("<FocusOut>",restore_entry)

		# maxDR
		def delete_entry(event):
			if self.maxdr_entry.get()=='Maximum DR':
				self.maxdr_entry.delete(0, "end")
				self.maxdr=''
				self.maxdr_entry.config(fg='black')
		def restore_entry(event):
			if self.maxdr_entry.get()=='':
				self.maxdr_entry.delete(0, "end")
				self.maxdr='Maximum DR'
				self.maxdr_entry.insert(0,self.maxdr)
				self.maxdr_entry.config(fg='gray45')
		self.maxdr=IntVar()
		self.maxdr='Maximum DR'
		maxdr=Label(self.frame2,text='Maximum DR',bg='white',fg='Black',font=('times new roman',15))
		maxdr.place(x=10,y=330)
		self.maxdr_entry=Entry(self.frame2,bg='lightgray',textvariable=self.maxdr)
		self.maxdr_entry.place(x=135,y=335,width=115)
		self.maxdr_entry.insert(0,self.maxdr)
		self.maxdr_entry.config(fg='gray45')
		self.maxdr_entry.bind("<FocusIn>",delete_entry)
		self.maxdr_entry.bind("<FocusOut>",restore_entry)

		# Minselfcal
		def delete_entry(event):
			if self.selfsnr_entry.get()=='Selfcal SNR':
				self.selfsnr_entry.delete(0, "end")
				self.selfsnr=''
				self.selfsnr_entry.config(fg='black')
		def restore_entry(event):
			if self.selfsnr_entry.get()=='':
				self.selfsnr_entry.delete(0, "end")
				self.selfsnr='Selfcal SNR'
				self.selfsnr_entry.insert(0,self.selfsnr)
				self.selfsnr_entry.config(fg='gray45')
		self.selfsnr=IntVar()
		self.selfsnr='Selfcal SNR'
		selfsnr=Label(self.frame2,text='Minimum Selfcal SNR',bg='white',fg='Black',font=('times new roman',15))
		selfsnr.place(x=270,y=330)
		self.selfsnr_entry=Entry(self.frame2,bg='lightgray',textvariable=self.selfsnr)
		self.selfsnr_entry.place(x=465,y=335,width=105)
		self.selfsnr_entry.insert(0,self.selfsnr)
		self.selfsnr_entry.config(fg='gray45')
		self.selfsnr_entry.bind("<FocusIn>",delete_entry)
		self.selfsnr_entry.bind("<FocusOut>",restore_entry)


		# Extra time
		def delete_entry(event):
			if self.extra_entry.get()=='Extra time':
				self.extra_entry.delete(0, "end")
				self.extra=''
				self.extra_entry.config(fg='black')
		def restore_entry(event):
			if self.extra_entry.get()=='':
				self.extra_entry.delete(0, "end")
				self.extra='Extra time'
				self.extra_entry.insert(0,self.extra)
				self.extra_entry.config(fg='gray45')
		self.extra=IntVar()
		self.extra='Extra time'
		extra=Label(self.frame2,text='Extra time averaging (s)',bg='white',fg='Black',font=('times new roman',15))
		extra.place(x=10,y=370)
		self.extra_entry=Entry(self.frame2,bg='lightgray',textvariable=self.extra)
		self.extra_entry.place(x=210,y=375,width=95)
		self.extra_entry.insert(0,self.extra)
		self.extra_entry.config(fg='gray45')
		self.extra_entry.bind("<FocusIn>",delete_entry)
		self.extra_entry.bind("<FocusOut>",restore_entry)

		# max avg time
		def delete_entry(event):
			if self.maxtime_entry.get()=='Max time':
				self.maxtime_entry.delete(0, "end")
				self.maxtime=''
				self.maxtime_entry.config(fg='black')
		def restore_entry(event):
			if self.maxtime_entry.get()=='':
				self.maxtime_entry.delete(0, "end")
				self.maxtime='Max time'
				self.maxtime_entry.insert(0,self.maxtime)
				self.maxtime_entry.config(fg='gray45')
		self.maxtime=IntVar()
		self.maxtime='Max time'
		maxtime=Label(self.frame2,text='Maximum average time (s)',bg='white',fg='Black',font=('times new roman',15))
		maxtime.place(x=285,y=370)
		self.maxtime_entry=Entry(self.frame2,bg='lightgray',textvariable=self.maxtime)
		self.maxtime_entry.place(x=509,y=375,width=75)
		self.maxtime_entry.insert(0,self.maxtime)
		self.maxtime_entry.config(fg='gray45')
		self.maxtime_entry.bind("<FocusIn>",delete_entry)
		self.maxtime_entry.bind("<FocusOut>",restore_entry)

		# weight
		def delete_entry(event):
			if self.weight_entry.get()=='uniform/natural/briggs':
				self.weight_entry.delete(0, "end")
				self.weight=''
				self.weight_entry.config(fg='black')
		def restore_entry(event):
			if self.weight_entry.get()=='':
				self.weight_entry.delete(0, "end")
				self.weight='uniform/natural/briggs'
				self.weight_entry.insert(0,self.weight)
				self.weight_entry.config(fg='gray45')
		self.weight=StringVar()
		self.weight='uniform/natural/briggs'
		weight=Label(self.frame2,text='Weighting mode',bg='white',fg='Black',font=('times new roman',15))
		weight.place(x=10,y=405)
		self.weight_entry=Entry(self.frame2,bg='lightgray',textvariable=self.weight)
		self.weight_entry.place(x=160,y=410,width=155)
		self.weight_entry.insert(0,self.weight)
		self.weight_entry.config(fg='gray45')
		self.weight_entry.bind("<FocusIn>",delete_entry)
		self.weight_entry.bind("<FocusOut>",restore_entry)

		# robust
		robusts=[float("{:0.1f}".format(i)) for i in np.arange(-1.0,1.1,0.1).tolist()]
		self.robust=DoubleVar()
		self.robust.set(robusts[-1])
		robust=Label(self.frame2,text='Robust',bg='white',fg='Black',font=('times new roman',15))
		robust.place(x=450,y=405)
		self.robust_options=ttk.Combobox(self.frame2,textvariable=self.robust,values=robusts,width=4,state='readonly')
		self.robust_options.place(x=525,y=410)

		for child in self.frame2.winfo_children():
			child.configure(state='disable')


		# HPC specific inputs

		# Frame 3
		self.frame3=Frame(self.root,bg='white',highlightthickness=3)
		self.frame3.place(x=970,y=500,width=600,height=410)
		title3=Label(self.frame3,text='HPC SETTINGS',bg='white',fg='black',font=('times new roman',25)).place(relx=0.5,y=20,anchor=CENTER)


		for child in self.frame3.winfo_children():
			child.configure(state='disable')
		
		# Validating inputs button
		self.button=Button(self.root,text='Validate Inputs',font=('times new roman',18),command=self.validate_input,height=2,width=15)
		self.button.place(x=30,y=770)

		# Save inputs button
		self.button=Button(self.root,text='Save Inputs',font=('times new roman',18),command=self.save_input,height=2,width=10)
		self.button.place(x=255,y=770)

		# Start button
		self.button=Button(self.root,text='Run P-AIRCARS',font=('times new roman',18),command=self.run_paircars,height=2,width=15)
		self.button.place(x=420,y=770)

	def download_data(self):
		self.popupwin1()
			
	def popupwin1(self):
		self.top_window= Toplevel(self.root)
		self.top_window.overrideredirect(True) # turns off title bar, geometry
		self.top_window.geometry('900x350+30+100') # set new geometry
		self.root.withdraw()
		# make a frame for the title bar
		title_bar = Frame(self.top_window, relief='raised', bd=0)
		label=Label(title_bar,text='Download MWA data',font=('times new roman',25))
		label.pack(side=TOP,pady=5)
		# a canvas for the main area of the window
		global top1
		top1 = Canvas(self.top_window)
		# put a close button on the title bar
		var=tk.IntVar()
		close_button = Button(top1, text='Close this Window', command=lambda:[var.set(1),self.close_win(self.top_window)])
		# pack the widgets
		title_bar.pack(pady=5,side=TOP,fill="x")
		close_button.pack(side=BOTTOM)
		top1.pack(expand=1, fill="both")

		def delete_entry(event):
			if self.apikey_entry.get()=='MWA ASVO API key':
				self.apikey_entry.delete(0, "end")
				self.apikey=''
				self.apikey_entry.config(fg='black')
		def restore_entry(event):
			if self.apikey_entry.get()=='':
				self.apikey_entry.delete(0, "end")
				self.apikey='MWA ASVO API key'
				self.apikey_entry.insert(0,self.apikey)
				self.apikey_entry.config(fg='gray45')
		self.apikey='MWA ASVO API key'
		apikey=Label(top1,text='MWA ASVO API key *',fg='Black',font=('times new roman',15))
		apikey.place(x=10,y=10)
		self.apikey_entry=Entry(top1,bg='lightgray',textvariable=self.apikey)
		self.apikey_entry.place(x=230,y=13,width=650)
		if self.apikey_entry.get()=='':
			self.apikey_entry.insert(0,self.apikey)
		self.apikey_entry.config(fg='gray45')
		self.apikey_entry.bind("<FocusIn>",delete_entry)
		self.apikey_entry.bind("<FocusOut>",restore_entry)
		
		# Data directory
		def delete_entry(event):
			if self.datadir_entry.get()=='Name of the directory of data':
				self.datadir_entry.delete(0, "end")
				self.datadir=''
				self.datadir_entry.config(fg='black')
		def restore_entry(event):
			if self.datadir_entry.get()=='':
				self.datadir_entry.delete(0, "end")
				self.datadir='Name of the directory of data'
				self.datadir_entry.insert(0,self.datadir)
				self.datadir_entry.config(fg='gray45')
		self.datadir='Name of the directory of data'
		datadir=Label(top1,text='Data Directory *',fg='Black',font=('times new roman',15))
		datadir.place(x=10,y=50)
		button2=ttk.Button(top1,text='browse',command=self.diropen5)
		button2.place(x=800,y=50)
		self.datadir_entry=Entry(top1,bg='lightgray',textvariable=self.datadir)
		self.datadir_entry.place(x=180,y=55,width=600)
		if self.datadir_entry.get()=='':
			self.datadir_entry.insert(0,self.datadir)
		self.datadir_entry.config(fg='gray45')
		self.datadir_entry.bind("<FocusIn>",delete_entry)
		self.datadir_entry.bind("<FocusOut>",restore_entry)

		# Timerange
		def delete_entry(event):
			if self.timerange_entry.get()=='yy0/mm0/dd0/hh0:mm0:ss0.ff0~yy1/mm1/dd1/hh1:mm1:ss1.ff1':
				self.timerange_entry.delete(0, "end")
				self.timerange1=''
				self.tmrange_entry.config(fg='black')
		def restore_entry(event):
			if self.timerange_entry.get()=='':
				self.timerange_entry.delete(0, "end")
				self.timerange='yy0/mm0/dd0/hh0:mm0:ss0.ff0~yy1/mm1/dd1/hh1:mm1:ss1.ff1'
				self.timerange_entry.insert(0,self.timerange1)
				self.timerange_entry.config(fg='gray45')
		self.timerange1='yy0/mm0/dd0/hh0:mm0:ss0.ff0~yy1/mm1/dd1/hh1:mm1:ss1.ff1'
		timerange1=Label(top1,text='Time range *',fg='Black',font=('times new roman',15))
		timerange1.place(x=10,y=90)
		self.timerange_entry=Entry(top1,bg='lightgray',textvariable=self.timerange1)
		self.timerange_entry.place(x=180,y=95,width=700)
		if self.timerange_entry.get()=='':
			self.timerange_entry.insert(0,self.timerange1)
		self.timerange_entry.config(fg='gray45')
		self.timerange_entry.bind("<FocusIn>",delete_entry)
		self.timerange_entry.bind("<FocusOut>",restore_entry)

		# Cal download
		caldownload=Label(top1,text='Download calibration data',fg='Black',font=('times new roman',15))
		caldownload.place(x=10,y=130)
		self.caldownload=BooleanVar()
		self.caldownload.set(True)
		c=Checkbutton(top1,text='',fg='Black',font=('times new roman',15),variable=self.caldownload,onvalue=True,offvalue=False,highlightbackground = "white")
		c.place(x=245,y=130)

		# Project ID
		def delete_entry(event):
			if self.projectid_entry.get()=='G0002':
				self.projectid_entry.delete(0, "end")
				self.projectid=''
				self.projectid_entry.config(fg='black')
		def restore_entry(event):
			if self.projectid_entry.get()=='':
				self.projectid_entry.delete(0, "end")
				self.projectid='G0002'
				self.projectid_entry.insert(0,self.projectid)
				self.projectid_entry.config(fg='gray45')
		self.projectid='G0002'
		projectid=Label(top1,text='Project ID',fg='Black',font=('times new roman',15))
		projectid.place(x=300,y=130)
		self.projectid_entry=Entry(top1,bg='lightgray',textvariable=self.projectid)
		self.projectid_entry.place(x=400,y=135,width=100)
		if self.projectid_entry.get()=='':
			self.projectid_entry.insert(0,self.projectid)
		self.projectid_entry.config(fg='gray45')
		self.projectid_entry.bind("<FocusIn>",delete_entry)
		self.projectid_entry.bind("<FocusOut>",restore_entry)


		# Obs ID
		def delete_entry(event):
			if self.obsid_entry.get()=='OBSID1,OBSID2':
				self.obsid_entry.delete(0, "end")
				self.obsid=''
				self.obsid_entry.config(fg='black')
		def restore_entry(event):
			if self.obsid_entry.get()=='':
				self.obsid_entry.delete(0, "end")
				self.obsid='OBSID1,OBSID2'
				self.obsid_entry.insert(0,self.obsid)
				self.obsid_entry.config(fg='gray45')
		self.obsid='OBSID1,OBSID2'
		obsid=Label(top1,text='Observation ID',fg='Black',font=('times new roman',15))
		obsid.place(x=520,y=130)
		self.obsid_entry=Entry(top1,bg='lightgray',textvariable=self.obsid)
		self.obsid_entry.place(x=660,y=135,width=220)
		if self.obsid_entry.get()=='':
			self.obsid_entry.insert(0,self.obsid)
		self.obsid_entry.config(fg='gray45')
		self.obsid_entry.bind("<FocusIn>",delete_entry)
		self.obsid_entry.bind("<FocusOut>",restore_entry)


		var=tk.IntVar()
		button1= Button(top1, text="Start download", command=lambda:[var.set(1),self.start_download()])
		button1.pack(side=BOTTOM,pady=5)
		button1=close_button
		button1.wait_variable(var)
		self.root.deiconify()
		del self.datadir,self.apikey
		return

	def start_download(self):
		data_dir=''
		time_range=''
		api_key=''
		cal_download=True
		project_id='G0002'
		obs_id=''

		if self.apikey_entry.get()!='' and self.apikey_entry.get()!='MWA ASVO API key':
			api_key=self.apikey_entry.get()

		if self.datadir_entry.get()!='' and self.datadir_entry.get()!='Name of the directory of data':
			data_dir=self.datadir_entry.get()

		if self.timerange_entry.get()!='' and self.timerange_entry.get()!='yy0/mm0/dd0/hh0:mm0:ss0.ff0~yy1/mm1/dd1/hh1:mm1:ss1.ff1':
			time_range=self.timerange_entry.get()

		project_id=self.projectid_entry.get()
		cal_download=self.caldownload.get()
		if self.obsid_entry.get()!='' and self.obsid_entry.get()!='OBSID1,OBSID2':
			obs_id=self.obsid_entry.get()
		if api_key=='':
			messagebox.showerror("No API key", "Please provide API key")
		elif data_dir=='':
			messagebox.showerror("No Data Directory", "Please provide data directory")
		elif time_range=='':
			messagebox.showerror("No time range", "Please provide time range")
		else:
			subprocess.Popen(["start_download",api_key,data_dir,time_range,project_id,obs_id,str(cal_download)])
			self.close_win(self.top_window)
			self.root.deiconify()
		return 
			

	def _resize_image(self,event):
		new_width = event.width
		new_height = event.height
		self.image = self.img_copy.resize((new_width, new_height))
		self.background_image = ImageTk.PhotoImage(self.image)
		self.background.configure(image =  self.background_image)

	def diropen1(self):
		self.msdir=filedialog.askdirectory(initialdir='/',title='Choose Data Directory')
		if self.msdir_entry.get()!='':
			self.msdir_entry.delete(0, "end") # delete all the text in the entry
			self.msdir_entry.insert(0, '') #Insert blank for user input
		self.msdir_entry.insert(END,self.msdir)
		self.msdir_entry.config(fg='black')
	
	def diropen2(self):
		self.basedir=filedialog.askdirectory(initialdir='/',title='Choose Base Directory')
		if self.basedir_entry.get()!='':
			self.basedir_entry.delete(0, "end") # delete all the text in the entry
			self.basedir_entry.insert(0, '') #Insert blank for user input
		self.basedir_entry.insert(END,self.basedir)
		self.basedir_entry.config(fg='black')
	
	def diropen3(self):
		self.fimagedir=filedialog.askdirectory(initialdir='/',title='Choose Final Image Directory')
		if self.fimagedir_entry.get()!='':
			self.fimagedir_entry.delete(0, "end") # delete all the text in the entry
			self.fimagedir_entry.insert(0, '') #Insert blank for user input
		self.fimagedir_entry.insert(END,self.fimagedir)
		self.fimagedir_entry.config(fg='black')

	def diropen4(self):
		self.maskfile=filedialog.askdirectory(initialdir='/',title='Choose CASA mask')
		if self.mask_entry.get()!='':
			self.mask_entry.delete(0, "end") # delete all the text in the entry
			self.mask_entry.insert(0, '') #Insert blank for user input
		self.mask_entry.insert(END,self.maskfile)
		self.mask_entry.config(fg='black')

	def diropen5(self):
		self.datadir=filedialog.askdirectory(initialdir='/',title='Choose Data Directory')
		if self.datadir_entry.get()!='':
			self.datadir_entry.delete(0, "end") # delete all the text in the entry
			self.datadir_entry.insert(0, '') #Insert blank for user input
		self.datadir_entry.insert(END,self.datadir)
		self.datadir_entry.config(fg='black')
	

	def fileopen1(self,filename=''):
		if filename=='':
			if os.path.exists(self.loadinput_entry.get())==False:
				try:
					self.loadinput=filedialog.askopenfilename(filetypes =[('P-AIRCARS Input File', '*.paircars')],initialdir=os.getcwd(),title='Choose P-AIRCARS input file')
				except:
					self.loadinput=''
				if self.loadinput_entry.get()!='':
					self.loadinput_entry.delete(0, "end") # delete all the text in the entry
					self.loadinput_entry.insert(0, '') #Insert blank for user input
			else:
				self.loadinput=self.loadinput_entry.get()
		else:
			self.loadinput=filename
		if self.loadinput!=None:
			if self.loadinput_entry.get()=='Load P-AIRCARS input file.....' or self.loadinput_entry.get()!='':
				self.loadinput_entry.delete(0,"end")
			self.loadinput_entry.insert(END,self.loadinput)
			self.loadinput_entry.config(fg='black')
			if os.path.exists(self.loadinput)==False:
				self.loadinput_entry.insert(0,'P-AIRCARS input file does not exists.')
				self.loadinput_entry.config(fg='gray45')	
			else:
				inputdic=pickle.load(open(self.loadinput,'rb'))
				self.msdir_entry.delete(0,"end")
				self.msdir_entry.insert(END,inputdic['msdir'])
				self.msdir_entry.config(fg='black')
				self.basedir_entry.delete(0,"end")
				self.basedir_entry.insert(END,inputdic['basedir'])
				self.basedir_entry.config(fg='black')
				self.tmrange_entry.delete(0,"end")
				self.tmrange_entry.insert(END,inputdic['timerange'])
				self.tmrange_entry.config(fg='black')
				self.chanrange_entry.delete(0,"end")
				self.chanrange_entry.insert(END,inputdic['chanrange'])
				self.chanrange_entry.config(fg='black')
				self.fimagedir_entry.delete(0,"end")
				self.fimagedir_entry.insert(END,inputdic['final_image_dir'])
				self.fimagedir_entry.config(fg='black')
				self.cal_entry.delete(0,"end")
				self.cal_entry.insert(END,','.join(inputdic['calibrator_caltable']))
				self.cal_entry.config(fg='black')
				self.safety.set(inputdic['safety_factor'])
				self.quality.set(inputdic['quality_factor'])
				self.verbose.set(inputdic['verbose'])
				self.logger.set(inputdic['keep_logger'])
				self.interactive.set(inputdic['interactive'])
				self.dodecor.set(inputdic['do_decor_correction'])
				self.refant.set(inputdic['ref_ant'])
				self.mask_entry.delete(0,"end")
				if os.path.isdir(inputdic['maskfile']):
					self.mask_entry.insert(END,inputdic['maskfile'])
				else:
					self.mask_entry.insert(END,'')
				self.mask_entry.config(fg='black')
				self.maskstr_entry.delete(0,"end")
				self.maskstr_entry.insert(END,inputdic['maskstr'])
				self.maskstr_entry.config(fg='black')
				if inputdic['calc_image_parameters']==True and inputdic['calc_selfcalib_params']==True:
					self.autocal.set(True)
				else:
					self.autocal.set(False)
				self.wsclean_input.set(inputdic['use_wsclean'])
				self.automask.set(inputdic['want_auto_masking'])
				self.notification_input.set(inputdic['send_notification'])
				self.email_entry.delete(0,"end")
				self.email_entry.insert(END,inputdic['email'])
				self.email_entry.config(fg='black')
				self.bandpass.set(inputdic['do_bandpass'])
				self.polcal.set(inputdic['do_polcal'])
				self.freqint_options.set(inputdic['image_delta_freq'])
				self.timeint_options.set(inputdic['image_delta_time'])
				self.freqwid_options.set(inputdic['image_freq'])
				self.timewid_options.set(inputdic['image_time'])
				self.cpufrac_options.set(inputdic['cpu_frac'])
				self.clearscreen.set(inputdic['clear_screen'])
				self.xcut_entry.delete(0,"end")
				self.xcut_entry.insert(END,float(inputdic['cutoutbox'].split(',')[0]))
				self.xcut_entry.config(fg='black')
				self.savemodel.set(inputdic['savemodel'])
				self.saveresiduals.set(inputdic['saveresidual'])
				self.flag.set(inputdic['want_uvsub_flag'])
				self.ankflag.set(inputdic['use_ankflagger'])
				self.ishpc.set(inputdic['hpc_environment'])	
				if inputdic['calc_image_parameters']==True and inputdic['calc_selfcalib_params']==True:
					self.cellsize_entry.delete(0,"end")
					self.cellsize_entry.insert(END,inputdic['cellsize'])
					self.cellsize_entry.config(fg='black')
					self.imsize_entry.delete(0,"end")
					self.imsize_entry.insert(END,inputdic['imsize'][0])
					self.imsize_entry.config(fg='black')
					self.scales_entry.delete(0,"end")
					self.scales_entry.insert(END,','.join(inputdic['multiscale_scales']))
					self.scales_entry.config(fg='black')
					self.uvtaper_entry.delete(0,"end")
					self.uvtaper_entry.insert(END,inputdic['uvtaper'])
					self.uvtaper_entry.config(fg='black')
					self.sigma_options.set(inputdic['start_sigma'])
					self.step_options.set(inputdic['sigma_step'])
					self.minsigma_options.set(inputdic['min_sigma'])
					self.resfrac_options.set(inputdic['residual_frac'])
					self.uvrange_entry.delete(0,"end")
					self.uvrange_entry.insert(END,inputdic['uvrange_to_cal'])
					self.uvrange_entry.config(fg='black')
					self.skipfreq_options.set(inputdic['skip_freq'])
					self.skiptime_options.set(inputdic['skip_time'])
					self.minsnr_options.set(inputdic['gain_minsnr'])
					self.drrms_entry.delete(0,"end")
					self.drrms_entry.insert(END,inputdic['DR_delta_rms'])
					self.drrms_entry.config(fg='black')
					self.drneg_entry.delete(0,"end")
					self.drneg_entry.insert(END,inputdic['DR_delta_neg'])
					self.drneg_entry.config(fg='black')
					self.mindr_entry.delete(0,"end")
					self.mindr_entry.insert(END,inputdic['min_DR'])
					self.mindr_entry.config(fg='black')
					self.maxdr_entry.delete(0,"end")
					self.maxdr_entry.insert(END,inputdic['max_DR'])
					self.maxdr_entry.config(fg='black')
					self.selfsnr_entry.delete(0,"end")
					self.selfsnr_entry.insert(END,inputdic['min_selfcal_snr'])
					self.selfsnr_entry.config(fg='black')
					self.extra_entry.delete(0,"end")
					self.extra_entry.insert(END,inputdic['extra_time'])
					self.extra_entry.config(fg='black')
					self.maxtime_entry.delete(0,"end")
					self.maxtime_entry.insert(END,inputdic['max_time_avg'])
					self.maxtime_entry.config(fg='black')
					self.weight_entry.delete(0,"end")
					self.weight_entry.insert(END,inputdic['weight'])
					self.weight_entry.config(fg='black')
					self.robust_options.set(inputdic['robust'])
		else:
			self.loadinput_entry.insert(0,'Load P-AIRCARS input file.....')
			self.loadinput_entry.config(fg='gray45')
	
	def getdata(self):
		if os.path.exists('inputs.py'):
			os.system('rm -rf inputs.py')
		fil=open('inputs.py','w')
		if self.msdir_entry.get()!='Name of the directory of data.....':
			msdir=self.msdir_entry.get()
		else:
			msdir=''
		self.msdir_input=msdir
		fil.write('msdir=\''+msdir+'\'\n')
		if self.basedir_entry.get()!='Name of the base directory .....':
			basedir=self.basedir_entry.get()
		else:
			basedir=''
		self.basedir_input=basedir
		fil.write('basedir=\''+basedir+'\'\n')
		fil.write('paircars_dir=\''+basedir+'\'\n')
		if self.tmrange_entry.get()!='yy0/mm0/dd0/hh0:mm0:ss0.ff0~yy1/mm1/dd1/hh1:mm1:ss1.ff1,yy2/mm2/dd2/hh2:mm2:ss2.ff2~....':
			timerange=self.tmrange_entry.get()
		else:
			timerange=''
		if timerange!='':
			try:
				time_list=[]
				for i in timerange.split(','):
					l=i.split('~')
					for j in l:
						time_list.append(j)
				timerange_list_mjdsecs=sorted([float("{:.2f}".format(timestamp_to_mjdsec(i,format=0))) for i in time_list])
			except:
				timerange=''
		fil.write('timerange=\''+timerange+'\'\n')
		if self.chanrange_entry.get()!='ch0~ch1,ch2~ch3,....':
			chanrange=self.chanrange_entry.get()
		else:
			chanrange=''
		if chanrange!='':
			try:
				chanrange_list=chanrange.split(',')
				for chan in chanrange_list:
					s_chan=int(chan.split('~')[0])
					e_chan=int(chan.split('~')[-1])
			except:
				chanrange=''
		fil.write('chanrange=\''+chanrange+'\'\n')
		if self.fimagedir_entry.get()!='Name of the directory to store final images .....':
			final_image_dir=self.fimagedir_entry.get()
			savedir=final_image_dir
		else:
			final_image_dir=''
			savedir=''
		self.fimagedir_input=final_image_dir
		fil.write('final_image_dir=\''+final_image_dir+'\'\n')
		fil.write('savedir=\''+savedir+'\'\n')
		if self.caltable!='/path/to/caltable0,/path/to/caltable1,....':
			calibrator_caltable=self.cal_entry.get().split(',')
		else:
			calibrator_caltable=[]	
		fil.write('calibrator_caltable='+str(calibrator_caltable)+'\n')	
		safety_factor=self.safety.get()
		fil.write('safety_factor='+str(safety_factor)+'\n')
		quality_factor=self.quality.get()
		fil.write('quality_factor='+str(quality_factor)+'\n')
		verbose=self.verbose.get()
		fil.write('verbose='+str(verbose)+'\n')
		keep_logger=self.logger.get()
		fil.write('keep_logger='+str(keep_logger)+'\n')
		interactive=self.interactive.get()
		fil.write('interactive='+str(interactive)+'\n')
		do_decor_correction=self.dodecor.get()
		fil.write('do_decor_correction='+str(do_decor_correction)+'\n')
		ref_ant=self.refant.get()
		fil.write('ref_ant='+str(ref_ant)+'\n')
		if self.maskfile!='CASA mask path.....':
			maskfile=self.maskfile
		else:
			maskfile=''
		fil.write('maskfile=\''+str(maskfile)+'\'\n')
		if self.maskstr!='Mask string in CASA format.....':
			maskstr=self.maskstr
		else:
			maskstr=''
		fil.write('maskstr=\''+str(maskstr)+'\'\n')
		if self.autocal.get()==True:
			calc_image_parameters=True
			calc_selfcalib_params=True
		else:
			calc_image_parameters=False
			calc_selfcalib_params=False
		fil.write('calc_image_parameters='+str(calc_image_parameters)+'\n')
		fil.write('calc_selfcalib_params='+str(calc_selfcalib_params)+'\n')
		use_wsclean=self.wsclean_input.get()
		fil.write('use_wsclean='+str(use_wsclean)+'\n')
		want_auto_masking=self.automask.get()
		fil.write('want_auto_masking='+str(want_auto_masking)+'\n')
		send_notification=self.notification_input.get()
		fil.write('send_notification='+str(send_notification)+'\n')
		email=self.email_entry.get()
		if email=='Enter email address to send notifications.....':
			email=''
		fil.write('email=\''+str(email)+'\'\n')
		do_bandpass=self.bandpass.get()
		fil.write('do_bandpass='+str(do_bandpass)+'\n')
		do_polcal=self.polcal.get()
		fil.write('do_polcal='+str(do_polcal)+'\n')
		image_delta_freq=self.freqint_options.get()
		fil.write('image_delta_freq='+str(image_delta_freq)+'\n')
		image_delta_time=self.timeint_options.get()
		fil.write('image_delta_time='+str(image_delta_time)+'\n')
		image_freq=self.freqwid_options.get()
		fil.write('image_freq='+str(image_freq)+'\n')
		image_time=self.timewid_options.get()
		fil.write('image_time='+str(image_time)+'\n')
		cpu_frac=self.cpufrac_options.get()/100.0
		fil.write('cpu_frac='+str(cpu_frac)+'\n')
		fil.write('instance='+str(int((psutil.cpu_count()*cpu_frac)/1.5))+'\n')
		clear_screen=self.clearscreen.get()
		fil.write('clear_screen='+str(clear_screen)+'\n')
		if self.xcut_entry==3 or self.xcut_entry.get()=='':
			cutoutbox='3,3'
		else:
			cutoutbox=str(self.xcut_entry.get())+','+str(self.xcut_entry.get())
		fil.write('cutoutbox=\''+str(cutoutbox)+'\'\n')
		savemodel=self.savemodel.get()
		fil.write('savemodel='+str(savemodel)+'\n')
		saveresidual=self.saveresiduals.get()
		fil.write('saveresidual='+str(saveresidual)+'\n')
		want_uvsub_flag=self.flag.get()
		fil.write('want_uvsub_flag='+str(want_uvsub_flag)+'\n')
		use_ankflagger=self.ankflag.get()
		fil.write('use_ankflagger='+str(use_ankflagger)+'\n')
		hpc_environment=self.ishpc.get()
		fil.write('hpc_environment='+str(hpc_environment)+'\n')
		fil.seek(0)
		fil.close()
		
	def get_advanced_data(self):
		fil=open('inputs.py','a')
		if self.cellsize_entry.get()=='xxarcsec or xxarcmin':
			cellsize='20arcsec'
		else:
			cellsize=self.cellsize_entry.get()
		fil.write('cellsize=\''+str(cellsize)+'\'\n')
		if self.imsize_entry.get()=='Number of pixels':
			imsize=[1280]
		else:
			try:
				imsize=[int(self.imsize_entry.get())]
			except:
				imsize=[1280]
		fil.write('imsize='+str(imsize)+'\n')
		if self.scales_entry.get()=='0,3,6,9' or self.scales_entry.get()=='':
			multiscale_scales=[0,3,6,9]
		else:
			multiscale_scales=[int(i) for i in self.scales_entry.get().split(',')]
		fil.write('multiscale_scales='+str(multiscale_scales)+'\n')
		if self.uvtaper_entry.get()=='xxlambda or xxklambda':
			uvtaper=''
		else:
			uvtaper=self.uvtaper_entry.get()
		fil.write('uvtaper=\''+str(uvtaper)+'\'\n')
		start_sigma=float(self.sigma_options.get())
		fil.write('start_sigma='+str(start_sigma)+'\n')
		sigma_step=float(self.step_options.get())
		fil.write('sigma_step='+str(sigma_step)+'\n')
		min_sigma=float(self.minsigma_options.get())
		fil.write('min_sigma='+str(min_sigma)+'\n')
		residual_frac=float(self.resfrac_options.get())
		fil.write('residual_frac='+str(residual_frac)+'\n')
		if self.uvrange_entry.get()=='CASA uvrange format':
			uvrange_to_cal=''
		else:
			uvrange_to_cal=self.uvrange_entry.get()
		fil.write('uvrange_to_cal=\''+str(uvrange_to_cal)+'\'\n')
		skip_freq=float(self.skipfreq_options.get())
		fil.write('skip_freq='+str(skip_freq)+'\n')
		skip_time=float(self.skiptime_options.get())
		fil.write('skip_time='+str(skip_time)+'\n')
		gain_minsnr=float(self.minsnr_options.get())
		fil.write('gain_minsnr='+str(gain_minsnr)+'\n')
		if self.drrms_entry.get()=='DR rms step':
			DR_delta_rms=30.0
		else:
			DR_delta_rms=float(self.drrms_entry.get())
		fil.write('DR_delta_rms='+str(DR_delta_rms)+'\n')
		if self.drneg_entry.get()=='DR negative step':
			DR_delta_neg=30.0
		else:
			DR_delta_neg=float(self.drneg_entry.get())
		fil.write('DR_delta_neg='+str(DR_delta_neg)+'\n')
		if self.mindr_entry.get()=='Minimum DR':
			min_DR=35
		else:
			min_DR=int(self.mindr_entry.get())
		fil.write('min_DR='+str(min_DR)+'\n')
		if self.maxdr_entry.get()=='Maximum DR':
			max_DR=3500
		else:
			max_DR=int(self.maxdr_entry.get())
		fil.write('max_DR='+str(max_DR)+'\n')
		if self.selfsnr_entry.get()=='Selfcal SNR':
			min_selfcal_snr=4.0
		else:
			min_selfcal_snr=float(self.selfsnr_entry.get())
		fil.write('min_selfcal_snr='+str(min_selfcal_snr)+'\n')
		if self.extra_entry.get()=='Extra time':
			extra_time=5.0
		else:
			extra_time=float(self.extra_entry.get())
		fil.write('extra_time='+str(extra_time)+'\n')
		if self.weight_entry.get()=='uniform/natural/briggs':
			weight='briggs'
		else:
			weight=float(self.weight_entry.get())
		fil.write('weight=\''+str(weight)+'\'\n')
		robust=float(self.robust_options.get())
		fil.write('robust='+str(robust)+'\n')
		if self.maxtime_entry.get()=='Max time':
			maxtime_entry=5.0
		else:
			maxtime_entry=float(self.maxtime_entry.get())
		fil.write('max_time_avg='+str(maxtime_entry)+'\n')
		fil.seek(0)
		fil.close()

	def open_ncra(self):
		self.ncralink='http://www.ncra.tifr.res.in/ncra/ncra/research/research-at-ncra-tifr/research-areas/the-sun-and-the-heliosphere/SolarPhysics'
		webbrowser.open(self.ncralink, new=2)
		return

	def open_mwa(self):
		self.mwalink='https://www.mwatelescope.org'
		webbrowser.open(self.mwalink, new=2)
		return

	def validate_input(self,show_message=True):
		self.getdata()
		self.get_advanced_data()
		if self.msdir_input=='' or self.basedir_input=='' or self.msdir_input=='Name of the directory of data.....' or self.basedir_input=='Name of the base directory .....':
			messagebox.showerror("Validate inputs", "Mandatory inputs are missing")	
			return 1
		else:
			if show_message:	
				messagebox.showinfo("Validate inputs", "Inputs are correct")
			return 0

	def run(self):	
		cwd=os.getcwd()
		validate=self.validate_input(show_message=False)
		if validate==0:
			cwd=os.getcwd()
			if os.path.isdir(self.basedir_input)==False:
				os.makedirs(self.basedir_input)
		res='go'
		if len(glob.glob(self.basedir_input+'/.*paircars_running'))>0:
			res=messagebox.askquestion('P-AIRCARS Running', 'P-AIRCARS already running in the base directory. Do you really want to over run?')
		return res

	def run_paircars(self):
		keyres=self.run()
		cwd=os.getcwd()
		if keyres=='no':
			messagebox.showinfo("Run P-AIRCARS", "P-AIRCARS is not started.")
			return
		else:
			if self.clearscreen.get()==True:
				running_jobids=[i.split('/')[-1].split('_paircars_running')[0][1:] for i in glob.glob(self.basedir_input+'/.*paircars_running')]
				screen_list=[os.path.basename(i) for i in glob.glob('/var/run/screen/S-'+str(getpass.getuser())+'/*')]
				paircars_homedir=os.path.expanduser('~')+'/.paircars'
				job_id_file=os.path.expanduser('~')+'/paircars_job_id.p'
				if os.path.isfile(job_id_file)==False:
					job_ids=[]
				else:
					job_ids=pickle.load(open(job_id_file,'rb'))
				for jobs in running_jobids:
					if os.path.exists(self.basedir_input+'/'+str(jobs)+'_pids.log'):
						pids=np.loadtxt(self.basedir_input+'/'+str(jobs)+'_pids.log',unpack=True).astype('int')
						if pids.shape==():
							pids=np.array([int(pids)])
						if len(pids)>0:
							for pid in pids:
								a=os.system('kill -0 '+str(int(pid))+' >/dev/null 2>&1')
								if a==0:
									os.system('kill -0 '+str(int(pid))+' >/dev/null 2>&1')
								else:
									pass
					for screen_name in screen_list:
						if jobs in screen_name:
							a=os.system('screen -S '+screen_name+' -X quit >/dev/null 2>&1')
							if a!=0:
								a=os.system('screen -wipe '+screen_name+' >/dev/null 2>&1')
					if int(jobs) in job_ids:
						job_ids.remove(int(jobs))
						if os.path.isfile(paircars_homedir+'/'+str(jobs)+'_inputs.paircars'):
							os.system('rm -rf '+paircars_homedir+'/'+str(jobs)+'_inputs.paircars')
				pickle.dump(job_ids,open(job_id_file,'wb'))							
			os.system('rm -rf '+self.basedir_input+'/.Finished* '+self.basedir_input+'/.*_paircars_running '+self.basedir_input\
						+'/*P-AIRCARS_mainlog_* '+self.basedir_input+'/*.batch '+self.basedir_input+'/*.log')
			if os.path.isfile(self.basedir_input+'/selfcal_inputs.py'):
				os.system('rm -rf '+self.basedir_input+'/selfcal_inputs.py')
			a=os.system('cp -r inputs.py '+self.basedir_input+'/selfcal_inputs.py')
			if os.WEXITSTATUS(a)!=0:
				messagebox.showerror("Copy error", "Can not copy the input file in base directory. Check the write permission of the base directory.")	
				os.system('rm -rf inputs.py')
				return 
			else:
				os.chdir(self.basedir_input)
				sys.path.append(os.getcwd())
				a=os.system('validating_paircars_input')
				if os.WEXITSTATUS(a)!=0:
					messagebox.showerror("Input error", "Error in inputs. Check the terminal for more information.")
					os.system('rm -rf inputs.py')
					return
				elif os.path.isdir(self.msdir_input)==False:
					messagebox.showerror('Data error','Data directory does not exist')
					os.system('rm -rf inputs.py')
					return
				else:	
					msfiles=glob.glob(self.msdir_input+'/*.ms')
					if len(msfiles)==0:
						messagebox.showerror('No data','No measurement set present in data directory')
						os.system('rm -rf inputs.py')
						return
					job_id_file=os.path.expanduser('~')+'/paircars_job_id.p'
					job_ids=[]
					if os.path.isfile(job_id_file)==False:
						job_id=np.random.randint(10,1000)
						job_ids.append(job_id)
						pickle.dump(job_ids,open(job_id_file,'wb'))
					else:
						job_ids=pickle.load(open(job_id_file,'rb'))
						while True:
							job_id=np.random.randint(10,1000)
							if job_id in job_ids:
								continue
							else:
								job_ids.append(job_id)
								pickle.dump(job_ids,open(job_id_file,'wb'))
								break
					fil=open(self.basedir_input+'/selfcal_inputs.py','a')
					fil.write('job_id='+str(job_id)+'\n')
					fil.seek(0)
					fil.close()
					save_input_file=self.save_input(savefile=self.basedir_input+'/inputs.paircars')
					paircars_homedir=os.path.expanduser('~')+'/.paircars'
					if os.path.isdir(paircars_homedir)==False:
						os.makedirs(paircars_homedir)
					os.system('cp -r '+save_input_file+' '+paircars_homedir+'/'+str(job_id)+'_inputs.paircars')
					os.chdir(self.basedir_input)
					screen_name='P-AIRCARS_mainlog_'+str(job_id)
					cmd='start_paircars --fresh '+str(self.fresh.get())+' --restart '+str(self.restart.get())
					os.system('echo "'+cmd+'" > '+screen_name+'.batch')
					screen_cmd='sh '+screen_name+'.batch'
					os.system('screen -S '+screen_name+' -X quit')	
					time.sleep(0.5)
					os.system('screen -mdS '+screen_name)
					time.sleep(0.5)
					os.system('screen -S '+screen_name+' -X stuff \"'+screen_cmd+'\n"')	
					os.chdir(cwd)
					os.system('rm -rf inputs.py')
					messagebox.showinfo('P-AIRCARS','P-AIRCARS has started. \nJob ID : '+str(job_id))
					return

	def save_input(self,savefile=''):
		dic={}
		if self.msdir_input!='Name of the directory of data.....':
			dic['msdir']=self.msdir_input
		else:
			dic['msdir']=''
		if self.basedir_input!='Name of the base directory .....':
			dic['basedir']=self.basedir_entry.get()
		else:
			dic['basedir']=''
		if self.tmrange_entry.get()!='yy0/mm0/dd0/hh0:mm0:ss0.ff0~yy1/mm1/dd1/hh1:mm1:ss1.ff1,yy2/mm2/dd2/hh2:mm2:ss2.ff2~....':
			dic['timerange']=self.tmrange_entry.get()
		else:
			dic['timerange']=''
		if self.chanrange_entry.get()!='ch0~ch1,ch2~ch3,....':
			dic['chanrange']=self.chanrange_entry.get()
		else:
			dic['chanrange']=''
		if self.fimagedir_entry.get()!='Name of the directory to store final images .....':
			dic['final_image_dir']=self.fimagedir_entry.get()
			dic['savedir']=dic['final_image_dir']
		else:
			dic['final_image_dir']=''
			dic['savedir']=''
		if self.caltable!='/path/to/caltable0,/path/to/caltable1,....':
			dic['calibrator_caltable']=self.cal_entry.get().split(',')
		else:
			dic['calibrator_caltable']=[]	
		dic['safety_factor']=int(self.safety.get())
		dic['quality_factor']=int(self.quality.get())
		dic['verbose']=self.verbose.get()
		dic['keep_logger']=self.logger.get()
		dic['interactive']=self.interactive.get()
		dic['do_decor_correction']=self.dodecor.get()
		dic['ref_ant']=self.refant.get()
		if self.maskfile!='CASA mask path.....':
			dic['maskfile']=self.maskfile
		else:
			dic['maskfile']=''
		if self.maskstr!='Mask string in CASA format.....':
			dic['maskstr']=self.maskstr
		else:
			dic['maskstr']=''
		if self.autocal.get()==True:
			dic['calc_image_parameters']=True
			dic['calc_selfcalib_params']=True
		else:
			dic['calc_image_parameters']=False
			dic['calc_selfcalib_params']=False
		dic['use_wsclean']=self.use_wsclean
		dic['want_auto_masking']=self.automask.get()
		dic['send_notification']=self.send_notification
		email=self.email_entry.get()
		if email=='Enter email address to send notifications.....':
			email=''
		dic['email']=email
		dic['do_bandpass']=self.bandpass.get()
		dic['do_polcal']=self.polcal.get()
		dic['image_delta_freq']=float(self.freqint_options.get())
		dic['image_delta_time']=float(self.timeint_options.get())
		dic['image_freq']=float(self.freqwid_options.get())
		dic['image_time']=float(self.timewid_options.get())
		dic['cpu_frac']=float(self.cpufrac_options.get())
		dic['clear_screen']=self.clearscreen.get()
		if self.xcut_entry==3 or self.xcut_entry.get()=='':
			dic['cutoutbox']='3,3'
		else:
			dic['cutoutbox']=str(self.xcut_entry.get())+','+str(self.xcut_entry.get())
		dic['savemodel']=self.savemodel.get()
		dic['saveresidual']=self.saveresiduals.get()
		dic['want_uvsub_flag']=self.flag.get()
		dic['use_ankflagger']=self.ankflag.get()
		dic['hpc_environment']=self.ishpc.get()
		if dic['calc_image_parameters']==True and dic['calc_selfcalib_params']==True:
			if self.cellsize_entry.get()=='xxarcsec or xxarcmin':
				dic['cellsize']='20arcsec'
			else:
				dic['cellsize']=self.cellsize_entry.get()
			if self.imsize_entry.get()=='Number of pixels':
				dic['imsize']=[1280]
			else:
				try:
					dic['imsize']=[int(self.imsize_entry.get())]
				except:
					dic['imsize']=[1280]
			if self.scales_entry.get()=='0,3,6,9' or self.scales_entry.get()=='':
				dic['multiscale_scales']=[0,3,6,9]
			else:
				dic['multiscale_scales']=self.scales_entry.get().split(',')
			if self.uvtaper_entry.get()=='xxlambda or xxklambda':
				dic['uvtaper']=''
			else:
				dic['uvtaper']=self.uvtaper_entry.get()
			dic['start_sigma']=float(self.sigma_options.get())
			dic['sigma_step']=float(self.step_options.get())
			dic['min_sigma']=float(self.minsigma_options.get())
			dic['residual_frac']=float(self.resfrac_options.get())
			if self.uvrange_entry.get()=='CASA uvrange format':
				dic['uvrange_to_cal']=''
			else:
				dic['uvrange_to_cal']=self.uvrange_entry.get()
			dic['skip_freq']=float(self.skipfreq_options.get())
			dic['skip_time']=float(self.skiptime_options.get())
			dic['gain_minsnr']=float(self.minsnr_options.get())
			if self.drrms_entry.get()=='DR rms step':
				dic['DR_delta_rms']=30.0
			else:
				dic['DR_delta_rms']=float(self.drrms_entry.get())
			if self.drneg_entry.get()=='DR negative step':
				dic['DR_delta_neg']=30.0
			else:
				dic['DR_delta_neg']=float(self.drneg_entry.get())
			if self.mindr_entry.get()=='Minimum DR':
				dic['min_DR']=35
			else:
				dic['min_DR']=int(self.mindr_entry.get())
			if self.maxdr_entry.get()=='Maximum DR':
				dic['max_DR']=3500
			else:
				dic['max_DR']=int(self.maxdr_entry.get())
			if self.selfsnr_entry.get()=='Selfcal SNR':
				dic['min_selfcal_snr']=4.0
			else:
				dic['min_selfcal_snr']=float(self.selfsnr_entry.get())
			if self.extra_entry.get()=='Extra time':
				dic['extra_time']=5.0
			else:
				dic['extra_time']=float(self.extra_entry.get())
			if self.maxtime_entry.get()=='Max time':
				dic['max_time_avg']=10.0
			else:
				dic['max_time_avg']=float(self.maxtime_entry.get())
			if self.weight_entry.get()=='uniform/natural/briggs':
				dic['weight']='briggs'
			else:
				dic['weight']=float(self.weight_entry.get())
			dic['robust']=float(self.robust_options.get())
		if savefile=='':
			paircars_input_file=self.save_file()
		else:
			paircars_input_file=savefile
		try:
			pickle.dump(dic,open(paircars_input_file,'wb'))
			if savefile=='':
				messagebox.showinfo("Input file save", "Input file saved as : "+paircars_input_file)
		except:
			if savefile=='':
				messagebox.showerror("Input file save error", "Could not save input file : "+paircars_input_file)
		return paircars_input_file
		
	def save_file(self):
		if self.basedir=='Name of the base directory .....':
			save_dir=os.getcwd()
		else:
			save_dir=self.basedir
		self.input_file=filedialog.asksaveasfilename(initialfile=save_dir+'/inputs.paircars',defaultextension=".paircars",filetypes=[("P-AIRCARS input","*.paircars")])
		return self.input_file

	def open_logger(self):
		subprocess.call(["log_viewer",self.basedir_entry.get()+'/Logs_and_Errors'])
		return

	def close_win(self,top):
		top.destroy()
	
	def insert_val(self,e):
		e.insert(0,0)

	def popupwin(self):
		top= Toplevel(self.root)
		top.geometry("350x250")
		self.root.eval(f'tk::PlaceWindow {str(top)} center')
		label=Label(top,text='P-AIRCARS Job ID')
		label.pack(pady=20)
		self.popentry=IntVar()
		self.pop_entry= Entry(top,width= 15,textvariable=self.popentry)
		self.pop_entry.pack()
		var=tk.IntVar()
		button= Button(top, text="Ok", command=lambda:[var.set(1),self.close_win(top)])
		button.pack(pady=35, side= TOP)
		button.wait_variable(var)
		return self.popentry.get()
		
	def logview(self):
		if self.basedir_entry.get()=='Name of the base directory .....' or self.basedir_entry.get()=='':
			paircars_homedir=os.path.expanduser('~')+'/.paircars'
			job_id=self.popupwin()
			if os.path.isfile(paircars_homedir+'/'+str(job_id)+'_inputs.paircars'):
				self.fileopen1(paircars_homedir+'/'+str(job_id)+'_inputs.paircars')
				if os.path.isdir(self.basedir_entry.get()+'/Logs_and_Errors')==False:
					messagebox.showerror('No Logs','No logs available for P-AIRCARS Job with Job ID : '+str(job_id)+'\n')
					return
				else:
					p=Process(target=self.open_logger)
					p.start()
			else:
				messagebox.showerror('No Jobs','No P-AIRCARS Job running for Job ID : '+str(job_id)+'\n')
				return
		else:
			if os.path.isdir(self.basedir_entry.get()+'/Logs_and_Errors')==False:
				os.makedirs(self.basedir_entry.get()+'/Logs_and_Errors')
			p=Process(target=self.open_logger)
			p.start()
		return
		
		
if __name__=='__main__':
	paircars_homedir=os.path.expanduser('~')+'/.paircars'
	if os.path.isdir(paircars_homedir)==False:
		os.makedirs(paircars_homedir)
	root=Tk()
	obs=PAIRCARS_inputs(root)
	root.mainloop()
