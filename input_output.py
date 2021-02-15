#requirements : install openpyxl, pandas, numpy
import pandas as pd
import numpy as np
#Reading each sheet separately
df = pd.read_excel("Option_3_DataSet.xlsx", sheet_name='Vol and Swaps Rates usd')

annuity_df = pd.read_excel("Option_3_DataSet.xlsx", sheet_name='Option_3_DataSet')
annuity_df = annuity_df.drop('Currency', axis = 1)
#Editing the columns names beofre setting them
df.iloc[0] = df.iloc[0].str.replace('USSW', '')
df.iloc[0] = df.iloc[0].str.replace(' Curncy', '')
df.iloc[0] = df.iloc[0].str.replace('USSN', '')
df.iloc[0] = df.iloc[0].str.replace(' CMPN', '')

#Setting the right names for columns
df.columns = df.iloc[0]
df = df.rename(columns={'Time\\ Underlying or Vol': 'Time'})
#Setting right names for rows
df['Time'] = pd.to_datetime(df['Time'])
df = df.set_index('Time')
#Deleting useless columns
df = df.drop(df.columns[[221,220,219,218,217,216,215,214,213,212]], axis=1)
df = df.drop([0])

#Creating a dataset for swaprates
swaprates = df.iloc[:, :15]
#Creating a dataset for volatilities
dfvols = df.iloc[:, 15:211]
#Journalizing the values
dfvols = dfvols.applymap(lambda x: np.sqrt(x/252))
