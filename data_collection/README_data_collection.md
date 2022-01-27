<b>DATA COLLECTION SOFTWARE</b>

This directory contains software used for collecting data with the plethysmography system. The software was written using LabVIEW and is able to collect and display realtime chamber pressure data from two mice simultaneously. This software was designed for use with the USB-1208FS-Plus data aqcuisition unit, made by Measurement Computing (https://www.mccdaq.com/). Using a different data aqcuisition unit is possible, but may require alterations to the included LabVIEW .vi files.
<br>
<br>
<b>Installation</b>
<br>
You will first need to install LabVIEW 2017 Version 17 (32-bit) from National Instruments (https://www.ni.com/). Other versions of LabVIEW may or may not be compatible with this software.

If using the USB-1208FS-Plus data acquisition unit, you will also need to download the MCC DAQ Software bundle from Measurement Computing (https://www.mccdaq.com/). This bundle includes InstaCal, which is necessary for initial setup of the data acquisition unit, and the ULx for NI LabVIEW driver, which allows LabVIEW to control the data acquisition unit.

After setting up the data acquisition unit and connecting it to the computer which will collect the data, you must run the InstaCal program. It should automatically detect the USB-1208FS-Plus board and when prompted you should add it to known boards for this computer. This step must be performed before running the data collection program with LabVIEW, or you will get an error. You may also wish to run some diagnostic tests using InstaCal to verify your hardware is working correctly. 

Assuming you have downloaded the plethysmography repository from the main page, navigate to the "plethysmography-main > data_collection > mouse_breathing_monitor" directory. It is important that all files within this directory remain unaltered. The best way to avoid changing anything is to create a shortcut to the "Mouse Breathing Monitor V6.vi" file and move the shortcut to your desktop. The data collection software should launch when you open the shortcut file.
<br>
<br>
<b>Operation</b>
<br>
After opening the "Mouse Breathing Monitor V6" file, a window showing the data collection interface will load. This interface is split into two section, one for each input channel. The two channels operate independently, data collection can be started or stopped for one channel without impacting collection from the other.

For each channel, a real-time graph of the pressure transducer's voltage output is displayed on the right. Scaling of the time axis is done by changing the "Seconds to Display" setting at the top of the interface window. An indicator light to the left of the voltage axis will illuminate if the corresponding channel is collecting data. 

The Mouse ID input field next to each graph allows you to add a unique identifier for the mouse being measured. This ID will be added to the filename when the collected data is saved. 

The Calibration Volume (mL) field is used to add an event marker of the specified volume in milliliters to the datafile when the LOG button is pressed. By default the range of allowable volumes is 0 - 0.5 mL.

The Cuff field indicates the voltage output of the pressure transducer monitoring the balloon cuff. By default the indicator light will illuminate if the voltage increases past the value 3.0, suggesting the balloon is overinflated and may be affecting the mouse's ability to breath normally.

The START button will initiate data collection for the corresponding channel when pressed. To stop data collection, press the same button again, which is now labeled STOP.

When the data collection software is opened, LabVIEW will automatically start running the program and pressing either START button will initiate data collection. To stop running the program, press the STOP PROGRAM button at the top right of the interface. At this point, none of the buttons will function until the program is restarted. To make any coding changes to the program, you must stop the program first. To restart the program, press the white arrow labeled "Run" from the toolbar at the top of the window.

Collected data is stored to the desktop in a directory labeled "Mouse Breathing Data" and sorted by the data and time when collection was started. The data is stored in a RCALV file, which has a ".rcalv" file extension. The filename includes the data, time, mouse ID, and channel number used for collection. The RCALV file is structured specifically so it can be loaded by the data analysis software. You can open it with a simple text editor to examine the raw data, but be careful not to change the contents, or the analysis software may fail to load it properly.
<br>
<br>
<b>Troubleshooting</b>
<br>
"Dependency loaded from new path" warning: This likely means that one of the LabVIEW .vi files has been moved or deleted from the "mouse_breathing_monitor" directory. Either redownload the repository to fix the directory or tell LabVIEW where to find the moved files when prompted.

"Error 100000 occurred at ULx Stop Task.vi": Check that your data acquisition unit is plugged into the computer and has gone through the one-time setup with InstaCal BEFORE opening LabVIEW or the data collection software.
