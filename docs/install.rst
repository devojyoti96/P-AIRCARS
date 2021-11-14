Installation
============
P-AIRCARS is a self-contained python package. It is build on python3. Python version greater than 3.6 is required for P-AIRACRS. It is recommended to use python virtual environment to install P-AIRCARS, but it is not a necessary requirement. Installation steps are as follows.

Creating virtual environment
****************************
Install python virtual environment:

>> python3 -m pip install --upgrade pip
 
>> python3 -m pip install virtualenv

Activate virtual environment:

>> python3 -m venv /path/to/virtualenv/paircars

>> source /path/to/virtualenv/paircars

Check the virtual environment path:

>> which python3

It should show "/path/to/virtualenv/paircars/bin/python3"

Obtaining P-AIRCARS
*******************
P-AIRCARS source code can be downloaded from https://github.com/devojyoti96/P-AIRCARS.git. 

P-AIRCARS is not public now. It will be public very soon. If you want to use P-AIRCARS before it becomes public, please reach us at dkansabanik@ncra.tifr.res.in, paircarsnotification@gmail.com

Installing P-AIRCARS
********************
Move to the P-AIRCARS directory.

Run *setup.py*

>> python3 setup.py install

That's all. It will automatically install all required packages and libraries. This installation process is in very early stage. If you find any issues during installation, please reach us at dkansabanik@ncra.tifr.res.in, paircarsnotification@gmail.com. 

.. note::
   If you installed P-AIRCARS in a virtual environment, please make sure you are inside the virtual environment. P-AIRCARS can not be accessed outside the virtual environment in this case.

   >> which python3

   It should show "/path/to/virtualenv/paircars/bin/python3". If you are not in virtual environment, please activate it. 

