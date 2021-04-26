from setuptools import setup
import os

os.system('cp -r scripts/run_intensity_selfcal.py scripts/run_intensity_selfcal')
os.system('cp -r scripts/run_bandpass_selfcal.py scripts/run_bandpass_selfcal')
os.system('cp -r scripts/run_pol_selfcal.py scripts/run_pol_selfcal')

setup(
    name='paircars',
    version='1.1.0',
    packages=['paircars'],
    author='Devojyoti Kansabanik',
    author_email='Andrew.Williams@curtin.edu.au',
    description='PAIRCARS',
    scripts=['scripts/run_intensity_selfcal','scripts/run_bandpass_selfcal','scripts/run_pol_selfcal'],
    install_requires=["numpy", "astropy", "skyfield", "matplotlib", "scipy>=0.15.1", "h5py"],
    extras_require={'skymap':["ephem", "Pillow"]}   # Needed only to generate sky maps in mwa_pb/skymap.py
)

os.system('rm -rf scripts/run_intensity_selfcal scripts/run_bandpass_selfcal scripts/run_pol_selfcal')

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

