DATA ANALYSIS SOFTWARE


Installation Instructions for Windows
=====================================

1. Download the latest version of Python 3 from http://www.python.org (Python 3.10.1 as of 12/20/2021) and install it. Make sure to select the box to "Add Python to PATH" during the installation process to make it easier for your computer to locate the installed program. Once the installation finishes, there should be a "Setup was successful" confirmation window. You can ignore the option to "Disable path length limit."


2. Download the plethysmography project files from this github page (https://github.com/bjrenquist/plethysmography) by clicking the green "Code" button and selecting "Download ZIP". Right click and "Extract All" to unzip the downloaded folder and move it to the location of your choice (i.e. Desktop, Documents, etc.). Open this folder and navigate to the subfolder "plethysmography-main\data_analysis\respiratory_analysis_v1_5", which contains four Python (.py) files. Create a shortcut to the "respiration_analysis_v1-5.py" file and place the shortcut on your Desktop for quick access.


3. Install the Python modules that are required to run the respiration_analysis software. This step can be intimidating, but if you are careful and double check spelling at each step, it should go smoothly.
   

	a. Click on the Windows Start Menu icon and start typing "Command Prompt". Right click the option with the same name that Windows suggests, and select "Run as administrator". Allow changes to this system when the warning dialog appears.


	b. Check if Python is setup correctly. Type the following into the prompt and hit Enter:
		
		python --version

	  The current version of Python should be displayed if all is well (i.e. Python 3.8.3).


	c. Update Python's built-in module installer, pip. Enter the following at the command prompt:

		python -m pip install -U pip

	  This will update pip to the latest version available.


	d. Install the following modules by entering the following lines one at a time:

		pip install openpyxl

		pip install lxml

		pip install numpy

		pip install scipy

		pip install matplotlib

	  You should see a confirmation that each module was successfully installed.


4. Now that all the required modules are installed, you should be able to run the respiration_analysis
   program using the shortcut that you created in step 2. When you run it, two windows should open: a
   black Windows terminal and the Respiration Analysis window with Voltage and Pressure graphs. The
   Windows terminal can be largely ignored. If you run into a software bug during analysis, an error
   should be displayed in this window (useful for troubleshooting).
