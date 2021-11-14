from setuptools import setup,find_packages
import os,sys,shutil,subprocess,glob,pip
from distutils.sysconfig import get_python_lib

os.environ['PATH']='/usr/local/bin:/usr/local/sbin:/bin:/usr/bin:/sbin:/usr/sbin'

cwd=os.getcwd()
if os.path.exists('mwa_pb/data/mwa_full_embedded_element_pattern.h5')==False:
	print ('Download MWA primary beam data.......\n')
	os.chdir('mwa_pb/data')
	os.system('wget -q -c http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5 --no-check-certificate')
	os.chdir(cwd)
def install(package):
    if hasattr(pip, 'main'):
        pip.main(['install', package])
    else:
        pip._internal.main(['install', package])

cwd=os.getcwd()
# Installing Libraries locally
install_ini_dir=os.path.dirname(get_python_lib())
if os.path.isdir(install_ini_dir+'/paircars_libraries')==False:
	os.makedirs(install_ini_dir+'/paircars_libraries')
try:
	os.system('rm -rf '+install_ini_dir+'/paircars_libraries/*')
	print ('Installing libraries.....\n')
	os.system('cp -r libraries/local '+install_ini_dir+'/paircars_libraries/local')
except:
	pass
os.chdir(install_ini_dir+'/paircars_libraries')
pwd=os.getcwd()
install_dir=pwd+'/local'
if os.path.isdir(install_dir)==False:
	os.makedirs(install_dir)


a=os.system(install_dir+'/bin/fftw-wisdom > tmp')
if a!=0:
	print ('Installing FFTW....\n')
	if os.path.isdir('fftw-3.3.8')==False:
		os.system('wget -q -c http://www.fftw.org/fftw-3.3.8.tar.gz --no-check-certificate >> tmp')
		shutil.unpack_archive('fftw-3.3.8.tar.gz')
		os.system('rm -rf fftw-3.3.8.tar.gz')
	os.chdir('fftw-3.3.8')
	os.system('make clean >> tmp')
	os.system('./configure --prefix='+install_dir+' --enable-threads --enable-openmp  --enable-shared >> tmp')
	os.system('make >> tmp')
	os.system('make install >> tmp')
	os.system('rm -rf tmp')
	os.chdir(pwd)

if os.path.exists(install_dir+'/lib/liblapack.so')==False:
	print ('Installing BLAS and LAPACK......\n')
	if os.path.isdir('lapack-3.10.0')==False:
		os.system('wget -q -c https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz --no-check-certificate >> tmp')
		shutil.unpack_archive('v3.10.0.tar.gz')
		os.system('rm -rf v3.10.0.tar.gz')
	os.chdir('lapack-3.10.0')
	if os.path.isdir('build')==False:
		os.mkdir('build')
	os.chdir('build')
	os.system('rm -rf *')
	os.system('cmake -DCMAKE_INSTALL_LIBDIR='+install_dir+'/lib -DBUILD_SHARED_LIBS=ON .. >> tmp')
	os.system('make >> tmp')
	os.system('make install >> tmp')
	os.system('rm -rf tmp')
	os.chdir(pwd)

if os.path.exists(install_dir+'/lib/libwcs.so')==False:
	print ('Installing wcslib....\n')
	if os.path.isdir('wcslib-4.24')==False:
		os.system('wget -q -c http://ftp.debian.org/debian/pool/main/w/wcslib/wcslib_4.24.orig.tar.bz2 --no-check-certificate')
		shutil.unpack_archive('wcslib_4.24.orig.tar.bz2')
		os.system('rm -rf wcslib_4.24.orig.tar.bz2')
	os.chdir('wcslib-4.24')
	os.system('./configure --without-pgplot --prefix='+install_dir+' >> tmp')
	os.system('make >> tmp')
	os.system('make install >> tmp')
	os.system('rm -rf tmp')
	os.chdir(pwd)

if os.path.exists(install_dir+'/bin/bison')==False:
	print ('Installing bison....\n')
	if os.path.isdir('bison-2.3')==False:
		os.system('wget -q -c http://ftp.gnu.org/gnu/bison/bison-2.3.tar.gz --no-check-certificate')
		shutil.unpack_archive('bison-2.3.tar.gz')
		os.system('rm -rf bison-2.3.tar.gz')
	os.chdir('bison-2.3')
	os.system('./configure --prefix='+install_dir+' >> tmp')
	os.system('make >> tmp')
	os.system('make install >> tmp')
	os.system('rm -rf tmp')
	os.chdir(pwd)

if os.path.exists(install_dir+'/bin/flex')==False:
	print ('Installing flex......\n')
	if os.path.isdir('flex-2.6.0'):
		os.system('wget -q -c https://liquidtelecom.dl.sourceforge.net/project/flex/flex-2.6.0.tar.bz2 --no-check-certificate')
		shutil.unpack_archive('flex-2.6.0.tar.bz2')
		os.system('rm -rf flex-2.6.0.tar.bz2')
	os.chdir('flex-2.6.0')
	os.system('./configure --prefix='+install_dir+' >> tmp')
	os.system('make >> tmp')
	os.system('make install >> tmp')
	os.system('rm -rf tmp')
	os.chdir(pwd)

if os.path.exists(install_dir+'/lib/libcasa_casa.so')==False:
	if os.path.isdir('casacore-3.1.1')==False:
		print ('Installing CASAcore....\n')
		os.system('wget -q -c https://github.com/casacore/casacore/archive/v3.1.1.tar.gz --no-check-certificate >> tmp')
		shutil.unpack_archive('v3.1.1.tar.gz')
		os.system('rm -rf v3.1.1.tar.gz')
	os.chdir('casacore-3.1.1')
	if os.path.isdir('build')==False:
		os.mkdir('build')
	os.chdir('build')
	os.system('rm -rf *')
	os.system('make clean >> tmp')
	os.system('export LD_LIBRARY_PATH='+install_dir+'$LD_LIBRARY_PATH')
	os.system('cmake ../ -DCMAKE_PREFIX_PATH='+install_dir+' -DCMAKE_INSTALL_PREFIX='+install_dir+' -DBUILD_PYTHON3=OFF -DBUILD_PYTHON=OFF >> tmp')
	os.system('make >> tmp')
	os.system('make install >> tmp')
	os.system('rm -rf tmp')
	os.chdir(pwd)

if os.path.exists(install_dir+'/lib/libcfitsio.so')==False:
	print ('Installing CFITSIO.....\n')
	if os.path.isdir('cfitsio-4.0.0')==False:
		os.system('wget -q -c http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.0.0.tar.gz --no-check-certificate >> tmp')
		shutil.unpack_archive('cfitsio-4.0.0.tar.gz')
		os.system('rm -rf cfitsio-4.0.0.tar.gz')
	os.chdir('cfitsio-4.0.0')
	os.system('./configure --prefix='+install_dir+' >> tmp')
	os.system('make shared >> tmp')
	os.system('make install >> tmp')
	os.system('rm -rf tmp')
	os.chdir(pwd)

if len(glob.glob(install_dir+'/lib/libboost*.so'))==0:
	print ('Installing BOOST....\n')
	if os.path.isdir('boost_1_77_0')==False:
		os.system('wget -q -c https://boostorg.jfrog.io/artifactory/main/release/1.77.0/source/boost_1_77_0.tar.gz --no-check-certificate >> tmp')
		shutil.unpack_archive('boost_1_77_0.tar.gz')
		os.system('rm -rf boost_1_77_0.tar.gz')
	os.chdir('boost_1_77_0')
	os.system('./bootstrap.sh --prefix='+install_dir+' --includedir='+install_dir+'/include --libdir='+install_dir+'/lib >> tmp')
	os.system('./b2 --prefix='+install_dir+' --link=shared --threading=multi install >> tmp')
	os.system('rm -rf tmp')
	os.chdir(pwd)

if os.path.exists(install_dir+'/lib/libgsl.so')==False:
	print ('Installing GSL....\n')
	if len(glob.glob('gsl-*'))==0:
		os.system('wget -q -c https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz --no-check-certificate >> tmp')
		shutil.unpack_archive('gsl-latest.tar.gz')
		os.system('rm -rf gsl-latest.tar.gz')
	gsl_dir=glob.glob('gsl-*')[0]
	os.chdir(gsl_dir)
	os.system('./configure --prefix='+install_dir+' >> tmp')
	os.system('make >> tmp')
	os.system('make install >> tmp')
	os.system('rm -rf tmp')
	os.chdir(pwd)

a=os.system('wsclean > wsclean_test')
os.system('rm -rf wsclean_test')
if a!=0:
	print ('Installing WSClean....\n')
	if os.path.isdir('wsclean-2.7')==False:
		os.system('wget -q -c https://sourceforge.net/projects/wsclean/files/wsclean-2.7/wsclean-2.7.tar.bz2 --no-check-certificate >> tmp')
		shutil.unpack_archive('wsclean-2.7.tar.bz2')
		os.system('rm -rf wsclean-2.7.tar.bz2')
	os.chdir('wsclean-2.7')
	if os.path.isdir('build')==False:
		os.makedirs('build')
	os.chdir('build')
	os.system('cmake ../ -DCXX11=ON -DCMAKE_PREFIX_PATH='+install_dir+' -DCMAKE_INSTALL_PREFIX='+install_dir+' >> tmp')
	os.system('make >> tmp')
	os.system('make install >> tmp')
	os.system('rm -rf tmp')
	os.chdir(pwd)
os.system('rm -rf tmp')

try:
	import numpy as np
except:
	install('numpy>=1.19.0')
	import numpy as np

a=os.system('which screen > tmp 2> tmp.error')
if a!=0:
	if os.path.isdir('screen-4.8.0')==False:
		os.system('wget "https://ftp.gnu.org/gnu/screen/screen-4.8.0.tar.gz" --no-check-certificate')
		shutil.unpack_archive('screen-4.8.0.tar.gz')
		os.system('rm -rf screen-4.8.0.tar.gz')
	os.chdir('screen-4.8.0')
	os.system('./configure --prefix='+install_dir)
	os.system('make install && install -m 644 screen-4.8.0/etc/etcscreenrc '+install_dir+'/etc/screenrc')
	screen_path=install_dir+'/bin'
	np.save(cwd+'/paircars_client/screen_path',screen_path)
else:
	fil=open('tmp','r')
	screen_path=os.path.dirname(fil.read())
	fil.seek(0)
	fil.close()
	np.save(cwd+'/paircars_client/screen_path',screen_path)
os.system('rm -rf tmp tmp.error')

os.chdir(cwd)
cwd=os.getcwd()
LD_LIBRARY_PATH=install_dir+'/lib'
INCLUDE_PATH=install_dir+'/include/'

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
os.chdir(cwd)

np.save(cwd+'/aNKflag/LDPATH',LD_LIBRARY_PATH)
os.system('cp -r scripts/run_intensity_selfcal.py scripts/run_intensity_selfcal')
os.system('cp -r scripts/run_bandpass_selfcal.py scripts/run_bandpass_selfcal')
os.system('cp -r scripts/run_pol_selfcal.py scripts/run_pol_selfcal')
os.system('cp -r scripts/control_paircars.py scripts/control_paircars')
os.system('cp -r scripts/validating_paircars_input.py scripts/validating_paircars_input')
os.system('cp -r scripts/manage_database.py scripts/manage_database')
os.system('cp -r scripts/parallel_ms_split.py scripts/parallel_ms_split')
os.system('cp -r scripts/final_imaging.py scripts/final_imaging')
os.system('cp -r scripts/compress_caltables.py scripts/compress_caltables')
os.system('cp -r scripts/run_paircars.py scripts/run_paircars')
os.system('cp -r scripts/start_paircars.py scripts/start_paircars')
os.system('cp -r scripts/go_paircars.py scripts/go-paircars')
os.system('cp -r scripts/log_viewer.py scripts/log_viewer')
os.system('cp -r scripts/track_final_imaging.py scripts/track_final_imaging')
os.system('cp -r scripts/download_mwa_data.py scripts/download_mwa_data')
os.system('cp -r scripts/start_download.py scripts/start_download')
setup(
    name='paircars',
    version='1.0.0',
    packages=['paircars','aNKflag','CALIBRATE','mwa_pb','paircars_casatasks','mantaray'],
    package_data={'paircars':['libpaircars.so','MWA_OBSids.npy','flux_scale_polyfit.npy','*.png','*.jpeg'],'aNKflag':['*.c', '*.npy', 'ankflag', '*.h','*.dat'],\
					'CALIBRATE':['calibrate_tools/*'],'mwa_pb':['data/*.fits', 'data/*.txt', 'data/*.h5', 'data/*.fab', 'data/*.dat']},
    author='Devojyoti Kansabanik',
    author_email='dkansabanik@ncra.tifr.res.in',
    description='PAIRCARS',
    install_requires=["extension-helpers","bokeh==2.4.0","pyparsing==2.4.7","numpy>=1.19.0", "astropy==4.3", "skyfield", "matplotlib", "scipy>=0.15.1", "h5py","julian","psutil",\
					"casatools","casatasks","casadata","cmake",'requests>=2.18.3','websocket_client','colorama','shadems'],
    scripts=['scripts/run_intensity_selfcal','scripts/run_bandpass_selfcal','scripts/run_pol_selfcal','scripts/control_paircars','scripts/validating_paircars_input',\
		'scripts/manage_database','scripts/parallel_ms_split','scripts/final_imaging','scripts/compress_caltables','scripts/run_paircars','scripts/start_paircars',\
		'scripts/go-paircars','scripts/log_viewer','scripts/track_final_imaging','scripts/start_download','scripts/download_mwa_data','scripts/beam_correct_image_CASA_mwa.py',\
		'scripts/beam_correct_image_CASA.py','scripts/beam_correct_image_IAU.py','scripts/beam_correct_image.py','scripts/beamtest.py','scripts/beam_ra_dec.py',\
		'scripts/calc_jones.py','scripts/make_beam_test.py','scripts/mwa_sensitivity.py','scripts/plot_skymap.py','scripts/primarybeammap_tant_test.py','scripts/track_and_suppress.py'],
	 extras_require={'skymap':["ephem", "Pillow"]},python_requires='>=3.6.1',
	)

os.system('rm -rf scripts/parallel_ms_split scripts/final_imaging scripts/run_intensity_selfcal scripts/run_bandpass_selfcal scripts/run_pol_selfcal scripts/validating_paircars_input scripts/control_paircars scripts/manage_database scripts/compress_caltables scripts/run_paircars scripts/go-paircars scripts/start_paircars scripts/log_viewer scripts/track_final_imaging scripts/download_mwa_data scripts/start_download')

setup(name='mantaray-client',
      version='1.0.0',
      packages=['mantaray','mantaray.api','mantaray.scripts'],
      install_requires=['requests>=2.18.3',
                        'websocket_client',
                        'colorama'],
 	entry_points={'console_scripts': ['mwa_client = mantaray.scripts.mwa_client:main']},
	python_requires='>=3.6.1',
	)

setup(
    name='jprq',
    version='2.0.1',
    author='Azimjon Pulatov',
    author_email='azimjohn@yahoo.com',
    maintainer='Azimjon Pulatov',
    maintainer_email='azimjohn@yahoo.com',
    url='https://github.com/azimjohn/jprq-py',
    license='MIT',
    packages=['jprq'],
    entry_points={
        'console_scripts': [
            'jprq = jprq.main:main',
        ]
    },
    install_requires=['certifi==2019.9.11','websockets==9.1','aiohttp==3.7.4','bson~=0.5.10','click==8.0.3'],
    python_requires='>=3.6.1',
)

