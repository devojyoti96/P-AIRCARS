from tkinter import *	
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
import os,copy,numpy as np,webbrowser,pickle,paircars,glob,time,getpass,tkinter as tk,sys
from PIL import Image,ImageTk
imagedir=os.path.abspath(os.path.dirname(paircars.__file__))


class LogView(object):
	def __init__(self, master, path):
		self.nodes = dict()
		self.root = root
		self.file_path=path
		self.path=path
		self.root.geometry('1500x800')
		self.root.attributes('-alpha',0)
		self.root.title('P-AIRCARS Logger')
		self.tree =ttk.Treeview(self.root, selectmode ='extended')
		self.T = Text(self.root,height=900,width=900)
		self.l = Label(self.root, text = "Logger")
		self.tree.pack(side ='left',fill='y')
		self.verscrlbar = ttk.Scrollbar(self.root,orient ="vertical", command = self.tree.yview)
		self.verscrlbar.pack(side ='left', fill ='y')
		self.tree.configure(yscrollcommand = self.verscrlbar.set)
		self.tree.column("#0", width = 590, anchor ='w')
		self.tree.heading('#0', text='Log Files', anchor='c')		
		abspath = os.path.abspath(self.path)
		self.insert_node('', os.path.basename(abspath), abspath)
		self.refresh_button=Button(self.root,text='Refresh Logger',command=self.refresh_logger)
		self.refresh_button.pack(pady=10)
		self.refresh_tree_button=Button(self.root,text='Refresh Log files',command=self.refresh_tree)
		self.refresh_tree_button.pack(expand=True)
		self.tree.bind('<<TreeviewOpen>>', self.open_node)
		self.tree.bind('<Double-Button>', self.get_data)
		self.verscrlbar1 = ttk.Scrollbar(self.root,orient ="vertical", command = self.T.yview)
		self.verscrlbar1.pack(side ='right', fill ='y')
		self.T.configure(yscrollcommand = self.verscrlbar1.set)
		self.l.pack(expand=True)
		self.T.pack(expand=True)

	def insert_node(self, parent, text, abspath):
		node = self.tree.insert(parent, 'end', text=text, values=(abspath,),open=False)
		if os.path.isdir(abspath):
			self.nodes[node] = abspath
			self.tree.insert(node, 'end')

	def open_node(self, event):
		node = self.tree.focus()
		abspath = self.nodes.pop(node, None)
		if abspath:
			self.tree.delete(self.tree.get_children(node))
			for p in sorted(os.listdir(abspath)):
				self.insert_node(node, p, os.path.join(abspath, p))
	
	def get_data(self,event):
		val=(self.tree.focus())
		if self.tree.item(val)['values'][0]==self.path:
			self.file_path=''
		else:
			self.file_path=self.tree.item(val)['values'][0]
			if os.path.isdir(self.file_path)==False:
				fil=open(self.file_path,'r')
				text=fil.readlines()
				self.text_widget(text,os.path.basename(self.file_path))
				fil.close()
		return self.file_path
	
	def delete_entry(self):
		self.T.configure(state='normal')
		self.T.delete('1.0',END)		

	def refresh_logger(self):
		if os.path.isfile(self.file_path):
			try:
				fil=open(self.file_path,'r')
				text=fil.readlines()
				self.text_widget(text,os.path.basename(self.file_path))
				fil.close()
			except:
				messagebox.showerror('Open error', 'Error in opening log file : '+self.file_path)
		else:
			messagebox.showerror('No log file', 'Log file : '+self.file_path+' does not exist')

	def text_widget(self,file_text,file_name):
		self.delete_entry()
		for i in file_text:
			self.T.insert(END,i)
		self.T.see(END)
		self.T.configure(state='disabled')
		self.l.config(text=file_name)

	def refresh_tree(self):
		x=self.tree.get_children()
		for i in x:
			self.tree.delete(i)
		abspath = os.path.abspath(self.path)
		self.insert_node('', os.path.basename(abspath), abspath)
		self.tree.bind('<<TreeviewOpen>>', self.open_node)
		self.tree.bind('<Double-Button>', self.get_data)
		
if __name__ == '__main__':
	root = tk.Tk()
	app = LogView(root, path=sys.argv[1])
	root.mainloop()
