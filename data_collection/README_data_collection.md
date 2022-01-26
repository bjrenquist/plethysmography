<b>DATA COLLECTION</b>

This directory contains software used for collecting data with the plethysmography system. The software was written using LabVIEW and is able to collect and display realtime chamber pressure data from two mice simultaneously. This software was designed for use with the USB-1208FS-Plus data aqcuisition unit, made by Measurement Computing (https://www.mccdaq.com/). Using a different data aqcuisition unit is possible, but may require alterations to the included LabVIEW .vi files.

<b>Installation</b><br>
You will first need to install LabVIEW 2017 Version 17 (32-bit) from National Instruments (https://www.ni.com/). Other versions of LabVIEW may or may not be compatible with this software.

If using the USB-1208FS-Plus data acquisition unit, you will also need to download the MCC DAQ Software bundle from Measurement Computing (https://www.mccdaq.com/). This bundle includes InstaCal, which is necessary for initial setup of the data acquisition unit, and the ULx for NI LabVIEW driver, which allows LabVIEW to control the data acquisition unit.

After setting up the data acquisition unit and connecting it to the computer which will collect the data, you must run the InstaCal program. It should automatically detect the USB-1208FS-Plus board and when prompted you should add it to known boards for this computer. This step must be performed before running the data collection program with LabVIEW, or you will get an error. You may also wish to run some diagnostic tests using InstaCal to verify your hardware is working correctly. 

Assuming you have downloaded the plethysmography repository from the main page, navigate to the "plethysmography-main > data_collection > mouse_breathing_monitor" directory. It is important that all files within this directory remain unaltered. The best way to avoid changing anything is to create a shortcut to the "Mouse Breathing Monitor V6.vi" file and move the shortcut to your desktop. The data collection software should launch when you open the shortcut file.

<b>Operation</b><br>


<b>Troubleshooting</b><br>
