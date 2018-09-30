# 384well_qPCR_analysis
Taqman qPCR analysis from 384 well plates

This python program automates qPCR analysis, made for 384-well plate analysis.  
Program will match technical replicates by names, and gets the dCt.  
Two output sheets are generated: 
  1. ALL RESULTS - results of dCts without deleting any outliers or high control (hex/GAPDH) values
  2. FINAL RESULTS - results of dCts after deleting any outliers or high control (hex/GAPDH) values
  
  
  
Rules for qPCR analysis input:
 
1. name of FAM and HEX samples must be identical except that one ends with - FAM.txt, and the other - HEX.txt
	ex: pPCRname - FAM.txt 
		pPCRname - HEX.txt

2. All FAM and Hex raw data text files must be in a folder titled: raw_qPCR. Anf the python code must be in a folder above them. 
	example file paths: 
		Main folder > python code (qPCR analysis.py)
		Main folder > raw_qPCR > HEX.txt and FAM.txt 

3. name of naming sheet must be the same as the samples except begins with 'NAMING ' and ends with .csv
	ex: NAMING pPCRname.csv

4. Naming sheet columns must be:
	A. Position (position in 384 well)
	B. Sample_name (name from excel naming sheet)

5. naming sheet must be in raw_qPCR folder
	example file paths: 
		Main folder > python code (qPCR analysis.py)
		Main folder > raw_qPCR > NAMING sheet







