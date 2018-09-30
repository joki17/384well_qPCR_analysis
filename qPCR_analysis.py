import csv 
import os
import pandas as pd
import numpy as np

def hex_filename (fam_filename):
	x = fam_filename.find('FAM.txt') #finds string index before FAM.txt
	hex_file = fam_filename[:x] + 'HEX.txt'
	return hex_file

def get_df_fam (fam_file):
	#gets df, index column = Position, leaves out first row, header = first row
	df_fam = pd.read_table(fam_file, sep="\t", index_col=2, skiprows=1, header=0)
	return df_fam

def get_df_hex (hex_file):
	#gets df, index column = Position, leaves out first row, header = first row
	df_hex = pd.read_table(hex_file, sep="\t", index_col=2, skiprows=1, header=0)
	return df_hex

def get_mean_stdev (df_calcs):
	mean1 = df_calcs.groupby('names', as_index=False)['dCt'].mean()
	stdev1 = df_calcs.groupby('names', as_index=False)['dCt'].apply(np.std) #note: stdev1 = series

	results1 = pd.concat([mean1, stdev1], axis=1, sort=False)
	results1 = results1.rename(columns={0 : 'stdev_mean'})
	
	results1 = pd.merge(df_calcs, results1, left_on = 'names', right_on = 'names')
	results1 = results1.rename(columns={'dCt_y' : 'dCt_mean', 'dCt_x':'dCt'})
	return results1

def one_sample_df(all_samples_df):
	#get df with only column 0, ext gene name
	df1 = all_samples_df.iloc[:,0]

	array = []
	count_array = []

	#count starts at -1, to match df numbering
	count = -1
	for item in df1:
		count += 1

		if item not in array:
			array.append(item)
			count_array.append(count)	

	return (count_array)


#input file should be in same path as cwd, or else manually change fhand file path. 
dir_path = os.getcwd()
print ('cwd: ' + dir_path)

for dirname, dirs, files in os.walk(dir_path + '/raw_qPCR'):

	fam_files = list(files) 
	for fam_file in fam_files:
		if fam_file.endswith('FAM.txt'):
			print (fam_file)
			hex_file = hex_filename(fam_file)
			print (hex_file)


			with open('raw_qPCR/' + fam_file) as f, open('raw_qPCR/' + hex_file) as h:
				
				df_fam = get_df_fam(f)
				df_hex = get_df_hex(h)

				df_names=pd.read_csv('raw_qPCR/NAMING ' + fam_file[:-10] + '.csv', encoding = "utf-8", index_col=0)
				
				df_calcs = pd.DataFrame()
				df_calcs['names'] = df_names.Sample_name
				#input Cps into new df 
				df_calcs['target'] = df_fam.Cp
				df_calcs['ctrl'] = df_hex.Cp
				df_calcs['dCt'] = df_fam.Cp - df_hex.Cp
				
				#sort the new database from names column
				df_calcs = df_calcs.sort_values('names')


				#get mean and SD of values before deleting any outliers
				results1 = get_mean_stdev(df_calcs)
				print('results1 df: \n', results1)


				#DELETE high GAPDH control values
				#high gapdh dataframe if want to view deleted values
				df_hiGAPDH = results1[results1['ctrl'] > 25]

				#drop rows with GAPDH value < 25, keeps index same as results1 just missing vals
				df_goodGAPDH = results1[results1['ctrl'] < 25]
				df_goodGAPDH = df_goodGAPDH.drop(columns="dCt_mean")
				df_goodGAPDH = df_goodGAPDH.drop(columns="stdev_mean")
				#print ('good GAPDH: \n', df_goodGAPDH)


				df_calcs2 = pd.DataFrame()

				for i, g in df_goodGAPDH.groupby('names'):	# i = 'names' value, g = group of same name values 

					#https://stackoverflow.com/questions/23199796/detect-and-exclude-outliers-in-pandas-dataframe
					#only keeps values within 1.25 standard deviations of dCt_x column
					new = g[np.abs(g.dCt-g.dCt.mean()) <= (1.25*g.dCt.std())]

					df_calcs2 = df_calcs2.append(new)


				results2 = get_mean_stdev(df_calcs2)
				print ('results2 df: \n', results2)

				
				all_result_file = 'ALL RESULTS_' + fam_file[:-10] + '.csv'
				with open(all_result_file,'w') as csvfile:
					results2.to_csv(all_result_file, sep=',', index=True, encoding='utf-8', quoting=None)
				

				one_sample_array = one_sample_df(results2)

				finaldf_results3 = results2.loc[one_sample_array]

				finaldf_results3 = finaldf_results3.drop(columns="target")
				finaldf_results3 = finaldf_results3.drop(columns="ctrl")
				finaldf_results3 = finaldf_results3.drop(columns="dCt")

				print('results3 df: \n', finaldf_results3)


				final_result_file = 'FINAL RESULTS_' + fam_file[:-10] + '.csv'
				with open(final_result_file,'w') as csvfile:
					finaldf_results3.to_csv(final_result_file, sep=',', index=True, encoding='utf-8', quoting=None)








