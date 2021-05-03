from setuptools import setup
import os

os.system('cp -r scripts/run_intensity_selfcal.py scripts/run_intensity_selfcal')
os.system('cp -r scripts/run_bandpass_selfcal.py scripts/run_bandpass_selfcal')
os.system('cp -r scripts/run_pol_selfcal.py scripts/run_pol_selfcal')
os.system('cp -r scripts/control_paircars.py scripts/control_paircars')
os.system('cp -r scripts/validating_paircars_input.py scripts/validating_paircars_input')
os.system('cp -r scripts/manage_database.py scripts/manage_database')
setup(
    name='paircars',
    version='1.1.0',
    packages=['paircars'],
	package_data={'paircars':['libpaircars.so']},
    author='Devojyoti Kansabanik',
    author_email='Andrew.Williams@curtin.edu.au',
    description='PAIRCARS',
    scripts=['scripts/run_intensity_selfcal','scripts/run_bandpass_selfcal','scripts/run_pol_selfcal','scripts/control_paircars','scripts/validating_paircars_input','scripts/manage_database'],
    install_requires=["numpy", "astropy", "skyfield", "matplotlib", "scipy>=0.15.1", "h5py"],
    extras_require={'skymap':["ephem", "Pillow"]}   # Needed only to generate sky maps in mwa_pb/skymap.py
)

os.system('rm -rf scripts/run_intensity_selfcal scripts/run_bandpass_selfcal scripts/run_pol_selfcal scripts/validating_paircars_input scripts/control_paircars scripts/manage_database')

setup(
    name='paircars_casatasks',
    version='1.1.0',
    packages=['paircars_casatasks'],
    author='Devojyoti Kansabanik',
    author_email='Andrew.Williams@curtin.edu.au',
    description='PAIRCARS',
    install_requires=["numpy", "astropy", "skyfield", "matplotlib", "scipy>=0.15.1", "h5py"],
    extras_require={'skymap':["ephem", "Pillow"]}   # Needed only to generate sky maps in mwa_pb/skymap.py
)

setup(
    name='CALIBRATE',
    version='1.1.0',
    packages=['CALIBRATE'],
    package_data={'CALIBRATE':['calibrate_tools/*']},
    author='Devojyoti Kansabanik',
    author_email='Andrew.Williams@curtin.edu.au',
    description='PAIRCARS',
    install_requires=["numpy", "astropy", "skyfield", "matplotlib", "scipy>=0.15.1", "h5py"],
    extras_require={'skymap':["ephem", "Pillow"]}   # Needed only to generate sky maps in mwa_pb/skymap.py
)

