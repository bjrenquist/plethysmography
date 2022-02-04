Plethysmography System Repository
=================================

The purpose of this repository is to collect all software and 3D printer files used in the creation and use of the Renquist lab's leak-free head out plethysmography system for mice. To get started, make sure you are on the main page of the repository (https://github.com/bjrenquist/plethysmography), then click on the green "Code" button and select "Download ZIP". The entire repository will then download to your computer as a .zip directory. Extract this directory using your tool of choice and move the "plethysmography-main" directory to a safe, convenient location (Documents, Desktop, etc.).

The repository is split into three sections:

Printed Parts
-------------
This directory contains a set of Blender and STL files used to make the 3D printed parts of the plethysmography system. See the included README_printed_parts file for more information.

Data Collection
---------------
This directory contains software used for collecting data with the plethysmography system. The software was written using LabVIEW and is able to collect and display realtime chamber pressure data from two mice simultaneously. See the included README_data_collection file for more information.

Data Analysis
-------------
This directory contains software used for analyzing collected data. The software was written using Python and is able to automatically analyze the chamber pressure signal obtained with the Data Collection software. It identifies breathing landmarks for each breath and computes 1-minute window averages of common breath characteristics. The analysis is automatically exported to an Excel file. See the included README_data_analysis file for more information.
