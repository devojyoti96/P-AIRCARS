from setuptools import setup
import os,sys
python_version=float('.'.join(sys.version.split(' ')[0].split('.')[:-1]))
print ('Python version :'+str(sys.version.split(' ')[0]),python_version)
if python_version!=3.6 and python_version!=3.7:
	print ('Python version is less than 3.6 or grater than 3.7. aNKflag can only run with python 3.6 and 3.7\n')	
	os._exit(0)
try:
	import numpy as np
except:
	os.system('python3.6 -m pip install numpy')

cwd=os.getcwd()
LD_LIBRARY_PATH=cwd+'/gsl/lib'
INCLUDE_PATH=cwd+'/gsl/include/'

os.chdir(cwd+'/aNKflag')
makefil=open('Makefile','r')
lines=makefil.readlines()
for i in range(len(lines)):
	if 'GSL_INCLUDE_DIR=' in lines[i]:
		lines[i]='GSL_INCLUDE_DIR='+INCLUDE_PATH+'\n'

for i in range(len(lines)):
	if 'GSL_LIBRARIES=' in lines[i]:
		lines[i]='GSL_LIBRARIES=-L'+LD_LIBRARY_PATH+' -Wl,\"-R '+LD_LIBRARY_PATH+'\"\n'

makefil.close()

lines=lines[:18]
with open("Makefile", "w") as output:
    for line in lines:
        output.write(line)
output.close()


if os.path.isfile('ankflag')==True:
	os.system('make clean')
os.system('make')
np.save('LDPATH',LD_LIBRARY_PATH)
os.chdir(cwd)
try:
	import casatools
	print ('casatools is already installed\n')
except:
	os.system('python3.6 -m pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools --user')

try:
	import casatasks
	print ('casatasks is already installed\n')
except:
	os.system('python3.6 -m pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatasks --user')

os.system('rm -rf casa*log')

setup(
    name='aNKflag',
    version='1.0',
    packages=['aNKflag'],
	package_data={'aNKflag':['*.c', '*.npy', 'ankflag', '*.h','*.dat']},
    author='Apurba Bera, Python wrapper by Devojyoti Kansabanik',
    description='Flagger',
    install_requires=["numpy", "astropy", "matplotlib", "scipy>=0.15.1"],
    )

