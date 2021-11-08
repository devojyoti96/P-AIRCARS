import sys,os
from itertools import islice
from subprocess import Popen, PIPE
from textwrap import dedent
from threading import Thread

import tkinter as tk # Python 3
from tkinter import *
from tkinter import ttk
from queue import Queue, Empty # Python 3

def iter_except(function, exception):
    """Works like builtin 2-argument `iter()`, but stops on `exception`."""
    try:
        while True:
            yield function()
    except exception:
        return

class DisplaySubprocessOutputDemo:
	def __init__(self, root, apikey, destdir, timerange, projectid, obsid, cal_download=True):
		self.root = root
		self.root.title('MWA data download')
		self.root.geometry('1000x800')
		# start dummy subprocess to generate some output
		process_list=['download_mwa_data','--API_key='+str(apikey),'--dest_dir='+str(destdir),\
		 '--timerange='+str(timerange),'--cal_download='+str(cal_download)]
		if obsid!='':
			process_list.append('--cal_obsids='+str(obsid))
		if projectid!='':
			process_list.append('--project_ID='+str(projectid))
		self.process = Popen(process_list, stdout=PIPE)

		# launch thread to read the subprocess output
		#   (put the subprocess output into the queue in a background thread,
		#    get output from the queue in the GUI thread.
		#    Output chain: process.readline -> queue -> label)
		q = Queue(maxsize=1024)  # limit output buffering (may stall subprocess)
		t = Thread(target=self.reader_thread, args=[q])
		t.daemon = True # close pipe if GUI process exits
		t.start()

		self.main_frame=Frame(self.root)
		self.main_frame.pack(fill=BOTH,expand=1)
	
		self.my_canvas=Canvas(self.main_frame)
		self.my_canvas.pack(side=LEFT,fill=BOTH,expand=1)

		self.my_scrollbar=ttk.Scrollbar(self.main_frame,orient=VERTICAL,command=self.my_canvas.yview)
		self.my_scrollbar.pack(side=RIGHT,fill=Y)

		self.my_canvas.configure(yscrollcommand=self.my_scrollbar.set)
		self.my_canvas.bind('<Configure>',self.on_configure)

		self.second_frame=Frame(self.my_canvas)
		self.my_canvas.create_window((0,0),window=self.second_frame,anchor="nw")

		# show subprocess' stdout in GUI
		self.label=Label(self.second_frame, text="  ", font=("times new roman", 15),justify='left')
		self.label.pack(ipadx=4, padx=4, ipady=4, pady=4, fill='both')
		self.update(q) # start update loop

	def reader_thread(self, q):
		"""Read subprocess output and put it into the queue."""
		try:
			with self.process.stdout as pipe:
				for line in iter(pipe.readline, b''):
					q.put(line)
		finally:
			q.put(None)

	def update(self, q):
		"""Update GUI with items from the queue."""
		for line in iter_except(q.get_nowait, Empty): # display all content
			if line is None:
				self.quit()
				return
			else:
				self.label['text'] = str(self.label['text'])+line.decode() # update GUI
				break # display no more than one line per 40 milliseconds
		self.root.after(40, self.update, q) # schedule next update

	def quit(self):
		self.process.kill() # exit subprocess if GUI is closed (zombie!)
		self.root.destroy()
	
	def on_configure(self,event):
		# update scrollregion after starting 'mainloop'
		# when all widgets are in canvas
		self.my_canvas.configure(scrollregion=self.my_canvas.bbox('all'))


root = tk.Tk()
app = DisplaySubprocessOutputDemo(root,sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],cal_download=eval(sys.argv[6]))
root.protocol("WM_DELETE_WINDOW", app.quit)
# center window
root.eval('tk::PlaceWindow %s center' % root.winfo_pathname(root.winfo_id()))
root.mainloop()

