# to write files similar to IOS output files:

"""
author: Lu Guan
date: Oct. 06, 2020
about: This script is for processing RBR CTD data and producing .ctd files in IOS Header format.

Modified July 2021 - September 2021 by Samantha Huntington


"""
#globals().clear()

import sys
import os
import pyrsktools
import itertools
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import pyproj
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from copy import copy, deepcopy
from scipy import signal
import gsw
import xarray as xr
from matplotlib import pyplot as plt
import glob
from datetime import datetime
from datetime import timedelta
from decimal import Decimal
import random
import cartopy.crs as ccrs
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter, LatitudeLocator)
import openpyxl
from openpyxl import load_workbook






#----------------------- step 1: .rsk file in EPdesktop structure - Export profile data to .csv files from .rsk files--------------------------------

#run function to export files
#sh - call function at bottom of function

#EXPORT_FILES(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2019-107-example/', file = '066024_20190823_0915_CTD_Data.rsk', year = '2019', cruise_number='107')

def EXPORT_FILES(dest_dir, file, year, cruise_number, event_start):  # event_start is taken from the Event_List or LOG
    """
    Read in a rsk file and output in csv format
    Inputs:
        - folder, file, year, cruise_number: rsk-format file containing raw RBR data
    Outputs:
        - csv file: csv files containing the profile data
    """



    # function to export .csv files from .rsk file.
    filename = str(dest_dir) + str(file) #full path and name of .rsk file
    rsk = pyrsktools.open(filename) # load up an RSK

    #check the number of profiles
    n_profiles = len(list(rsk.profiles())) # get the number of profiles recorded


    ctd_data = pd.DataFrame() # set an empty pandas dataframe to store data

    #print(list(rsk.channels.items()))
    #export .csv file for each profile
    for i in range(0, n_profiles, 1):

        downcast = list(itertools.islice(rsk.casts(pyrsktools.Region.CAST_DOWN), i, i+1))[0].npsamples() #separate samples for each downcast file
        downcast_dataframe = pd.DataFrame(data=downcast, columns=downcast.dtype.names) # convert data into pandas data frame

        upcast = list(itertools.islice(rsk.casts(pyrsktools.Region.CAST_UP), i, i+1))[0].npsamples()
        upcast_dataframe = pd.DataFrame(data=upcast, columns=upcast.dtype.names)
        #print(downcast[0])

        column_names = list(downcast.dtype.names) #read the column names
        #print(list(rsk.channels))
        #print(column_names)
        #print(len(column_names))
        col_range = len(column_names)                                                           # added this to get them all
        column_names[0] = 'Time(yyyy-mm-dd HH:MM:ss.FFF)' #update time
        for j in range(1, col_range, 1):       # update the column names
            column_names[j]= column_names[j][0: -3] + "(" + list(rsk.channels.items())[j-1][1][4] + ")"

        downcast_dataframe.columns = column_names # update column names in downcast data frame
        downcast_dataframe["Cast_direction"] = "d"  # add a column for cast direction
        downcast_dataframe["Event"] = i+ event_start # add a column for event number - get event_start from Logs
        upcast_dataframe.columns = column_names
        upcast_dataframe["Cast_direction"] = "u"
        upcast_dataframe["Event"] = i+event_start
        downcast_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(i+event_start).zfill(4) + "_DOWNCAST.csv" #downcast file name
        upcast_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(i+event_start).zfill(4) + "_UPCAST.csv" #upcast file name
        profile_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(i+event_start).zfill(4) + ".csv" # profile name
        output_filename = year + "-" + cruise_number + "_CTD_DATA.csv" # all data file name

        profile_data = pd.concat([downcast_dataframe, upcast_dataframe]) #combine downcast and upcast into one profile
        ctd_data = ctd_data.append(profile_data, ignore_index=True) #combine profiles into one file


        #downcast_dataframe.to_csv(folder + downcast_name)
        #upcast_dataframe.to_csv(folder + upcast_name)
        profile_data.to_csv(dest_dir + profile_name) #export each profile

    ctd_data.to_csv(dest_dir + output_filename) #export all data in one .csv file

#EXPORT_FILES(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2021-037/', file = '204848_20210203_2358.rsk', year = '2021', cruise_number='037', event_start=3)

### ----------------------------  Exploring if we can input mulitple .rsk files----------------------

def EXPORT_MULTIFILES(dest_dir, year, cruise_number, all_last, event_start):  # ALL or LAST, depends on event log, see which profiles to extract
    """
    Read in a directory of rsk files and output in csv format
    Inputs:
        - folder, file, year, cruise_number: rsk-format file containing raw RBR data
    Outputs:
        - csv file: csv files containing the profile data
    """
    files = os.listdir(dest_dir)  # list all the files in dest_dir
    files = list(filter(lambda f: f.endswith('.rsk'), files))  # keep the rsk files only
    n_files = len(files) #get the number of files

    current_profile = event_start

    # function to export .csv files from .rsk files .
    for k in range(0, n_files, 1):
        filename = str(dest_dir) + str(files[k]) #full path and name of .rsk file
        rsk = pyrsktools.open(filename) # load up an RSK

        #check the number of profiles
        n_profiles = len(list(rsk.profiles())) # get the number of profiles recorded
        #last_profile = n_profiles[:-1] # get the last profile only
        last_profile = list(rsk.profiles())[:-1]




        ctd_data = pd.DataFrame() # set an empty pandas dataframe to store data

        if all_last == 'ALL':  # get all profiles within each .rsk file

        #export .csv file for each profile
            for i in range(0, n_profiles, 1):


                downcast = list(itertools.islice(rsk.casts(pyrsktools.Region.CAST_DOWN), i, i+1))[0].npsamples() #separate samples for each downcast file
                downcast_dataframe = pd.DataFrame(data=downcast, columns=downcast.dtype.names) # convert data into pandas data frame

                upcast = list(itertools.islice(rsk.casts(pyrsktools.Region.CAST_UP), i, i+1))[0].npsamples()
                upcast_dataframe = pd.DataFrame(data=upcast, columns=upcast.dtype.names)

                column_names = list(downcast.dtype.names) #read the column names

                col_range = len(column_names)                                                           # added this to get them all
                column_names[0] = 'Time(yyyy-mm-dd HH:MM:ss.FFF)' #update time

                for j in range(1, col_range, 1):       # update the column names
                    column_names[j]= column_names[j][0: -3] + "(" + list(rsk.channels.items())[j-1][1][4] + ")"
                downcast_dataframe.columns = column_names  # update column names in downcast data frame
                upcast_dataframe.columns = column_names

                downcast_dataframe["Cast_direction"] = "d"  # add a column for cast direction
                downcast_dataframe["Event"] = current_profile # add a column for event number - count the profiles

                upcast_dataframe["Cast_direction"] = "u"
                upcast_dataframe["Event"] = current_profile


                downcast_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(current_profile).zfill(4) + "_DOWNCAST.csv" #downcast file name
                upcast_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(current_profile).zfill(4) + "_UPCAST.csv" #upcast file name
                profile_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(current_profile).zfill(4) + ".csv" # profile name

                output_filename = year + "-" + cruise_number + "_CTD_DATA.csv" # all data file name

                current_profile = current_profile + 1  # sequential profiles through the files
                profile_data = pd.concat([downcast_dataframe, upcast_dataframe]) #combine downcast and upcast into one profile
                #profile_data['Event'] = e_num + 1

                #ctd_data = ctd_data.append(profile_data, ignore_index=True) #combine profiles into one file # not needed for multifile


                #downcast_dataframe.to_csv(folder + downcast_name)
                #upcast_dataframe.to_csv(folder + upcast_name)
                profile_data.to_csv(dest_dir + profile_name) #export each profile

            #ctd_data.to_csv(dest_dir + output_filename) #export all data in one .csv file # not needed for multifile

        elif all_last == 'LAST':  # get the last profile in each .rsk file  This is sometimes needed, depends on Logs
                # export .csv file for each profile

            i = n_profiles-1 # get the last profile only

            downcast = list(itertools.islice(rsk.casts(pyrsktools.Region.CAST_DOWN), i, i +1))[0].npsamples()  # separate samples for each downcast file
            downcast_dataframe = pd.DataFrame(data=downcast, columns=downcast.dtype.names)  # convert data into pandas data frame

            upcast = list(itertools.islice(rsk.casts(pyrsktools.Region.CAST_UP), i, i + 1))[0].npsamples()
            upcast_dataframe = pd.DataFrame(data=upcast, columns=upcast.dtype.names)

            column_names = list(downcast.dtype.names)  # read the column names
            col_range = len(column_names)  # added this to get them all
            column_names[0] = 'Time(yyyy-mm-dd HH:MM:ss.FFF)'  # update time

            for j in range(1, col_range, 1):  # update the column names
                column_names[j] = column_names[j][0: -3] + "(" + list(rsk.channels.items())[j - 1][1][4] + ")"


            downcast_dataframe.columns = column_names  # update column names in downcast data frame
            downcast_dataframe["Cast_direction"] = "d"  # add a column for cast direction
            downcast_dataframe["Event"] = current_profile  # add a column for event number
            upcast_dataframe.columns = column_names
            upcast_dataframe["Cast_direction"] = "u"
            upcast_dataframe["Event"] = current_profile
            downcast_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(current_profile).zfill(4) + "_DOWNCAST.csv"  # downcast file name
            upcast_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(current_profile).zfill(4) + "_UPCAST.csv"  # upcast file name
            profile_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(current_profile).zfill(4) + ".csv"  # profile name
            profile_number = len(list(profile_name)) / 2
            output_filename = year + "-" + cruise_number + "_CTD_DATA.csv"  # all data file name

            # event number counter - sequential with need to manually change to match event log if required
            current_profile = current_profile + 1
            profile_data = pd.concat([downcast_dataframe, upcast_dataframe])  # combine downcast and upcast into one profile
            #ctd_data = ctd_data.append(profile_data, ignore_index=True) #combine profiles into one file # not needed for multifile

            # downcast_dataframe.to_csv(folder + downcast_name)
            # upcast_dataframe.to_csv(folder + upcast_name)
            #ctd_data.to_csv((dest_dir + profile_name))
            profile_data.to_csv(dest_dir + profile_name)  # export each profile

            #ctd_data.to_csv(dest_dir + output_filename) #export all data in one .csv file # not needed for multifile

#EXPORT_MULTIFILES(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2021-039_Python_Run/', year = '2021', cruise_number='039', all_last='LAST', event_start=1)

#------------------------- step 1a - alternative if rsk file cant be used and an excel file is used instead...........

def READ_EXCELrsk(dest_dir, year, cruise_number, event_start, all_last):
   """
    #function to read in an excel (.xlxs) file exported from RUSKIN software using the rbr .rsk file.
   """

   files = os.listdir(dest_dir)  # list all the files in dest_dir
   files = list(filter(lambda f: f.endswith('.xlsx'), files))  # keep the rsk xlsx files only (make sure no other xlsx)
   n_files = len(files)  # get the number of files

   current_profile = event_start

   for k in range(0, n_files, 1):
       filename = str(dest_dir) + str(files[k])  # full path and name of .rsk file

       # extract a dataframe from the excel sheets
       df1 = pd.read_excel(filename, sheet_name='Data', skiprows=[0])  #engine=openpyxl
       #print('printing df1')
       #print(df1)


       #print(df1.sheetnames)

       df2 = pd.read_excel(filename, sheet_name='Profile_annotation', skiprows=[0])
       df3 = pd.read_excel(filename, sheet_name='Metadata', skiprows=[0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13], usecols=[2])
       #print(df1.keys)
       df1['Time'] = pd.to_datetime(df1['Time'])
       df2['Time 1'] = pd.to_datetime(df2['Time 1'])
       df2['Time 2'] = pd.to_datetime(df2['Time 2'])

       down_times = pd.DataFrame()
       up_times = pd.DataFrame()

       # find the start and end times for each profile
       down_times['Start'] = df2['Time 1'][1::3]
       down_times['End'] = df2['Time 2'][1::3]
       down_times.index = range(len(down_times))
       up_times['Start'] = df2['Time 1'][2::3]
       up_times['End'] = df2['Time 2'][2::3]
       up_times.index = range(len(up_times))

       n_times = len(list(down_times['Start']))

       ctd_data = pd.DataFrame()

       for i in range(0, n_times, 1):
           down_start_time = down_times['Start'][i]
           down_end_time = down_times['End'][i]
           up_start_time = up_times['Start'][i]
           up_end_time = up_times['End'][i]

           # extract data for each profile - using start and end times
           downcast = df1[(df1['Time'] > down_start_time) & (df1['Time'] <=down_end_time)]
           downcast['Cast_direction'] = 'd'
           # Add event numbers - need an event start number (i + 1 only works if one file and the events start at 1.
           # alter code to read in multiple files and create an Event_Start_List to loop through?
           # or if events aren't sequential a list of all event numbers to add at the end?
           downcast['Event'] = i+event_start     #get this from the log
           upcast = df1[(df1['Time'] > up_start_time) & (df1['Time'] <= up_end_time)]
           upcast['Cast_direction'] = 'u'
           upcast['Event'] = i+event_start       # get this from the log
           profile_data = pd.concat([downcast, upcast])  # combine downcast and upcast into one profile

           if all_last == 'LAST':
               profile_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(i + current_profile).zfill(
                   4) + ".csv"  # profile name #i + from log
           elif all_last == 'ALL':
               profile_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(i + event_start).zfill(
                   4) + ".csv"  # profile name #i + from log
               output_filename = year + "-" + cruise_number + "_CTD_DATA.csv"  # all data file name

               ctd_data = ctd_data.append(profile_data, ignore_index=True)  # combine profiles into one file
               profile_data.rename({'Time': 'Time(yyyy-mm-dd HH:MM:ss.FFF)'}, axis=1, inplace=True)
               profile_data.to_csv(dest_dir + profile_name)  # export each profile
               current_profile = current_profile + 1


           output_filename = year + "-" + cruise_number + "_CTD_DATA.csv"  # all data file name



           #output_filename = year + "-" + cruise_number + "_CTD_DATA.csv"  # all data file name#ctd_data.to_csv(dest_dir + output_filename)
       #ctd_data.to_csv(dest_dir + output_filename)

#READ_EXCELrsk(dest_dir='C:/Projects/RBR_CTD_Pycharm/2021-038/', year='2021', cruise_number='038', event_start = 341, all_last = 'ALL')

#----------------------- step 1b: .rsk file in full structure - combine profiles into one .csv file --------------------------------

#MERGE_FILES(dest_dir = '/home/guanl/Desktop/Projects/RBR/Processing/2020-085/CTD data/', year = '2020', cruise_number= '085')

def MERGE_FILES(dest_dir, year, cruise_number):
    """
    Read in multiple csv file and output one csv format
    Inputs:
        - folder, year, cruise_number
    Outputs:
        - csv file: csv files containing the profile data
   """
    files = os.listdir(dest_dir) # list all the files in dest_dir
    files = list(filter(lambda fname: 'profile' in fname, files)) # keep the profile csv files only
    event_csv = pd.read_csv(dest_dir + year + '-'
                                              + cruise_number + '_header-merge.csv')

    files.sort(key=lambda x: int(x[-8:-4])) #reorder the files according to profile number
    n_profiles = len(files)
    ctd_data = pd.DataFrame()  # set an empty pandas dataframe to store data
    for i in range(0, n_profiles, 1):
        input_filename = str(dest_dir) + str(files[i])
        #data = pd.read_csv(input_filename, sep=',', skiprows=range(0, 20), encoding= 'unicode_escape') #original by Lu
        data = pd.read_csv(input_filename, sep=',', encoding='unicode_escape')
        data['Event'] = int(files[i][-8:-4])
        #data['Event'] = event_csv['LOC:Event Number'][i] # Use this if you want non-sequential events, but won't work auto processing (cast i+1 etc...)
        #data['Time'] = pd.to_datetime(data['Time'])        #might not be needed
        ctd_data = pd.concat([ctd_data, data], ignore_index=True)
    if ctd_data.columns.to_list()[0] == '//Time(yyyy-mm-dd HH:MM:ss.FFF)':
        ctd_data.rename({'//Time(yyyy-mm-dd HH:MM:ss.FFF)': 'Time(yyyy-mm-dd HH:MM:ss.FFF)'}, axis=1, inplace=True)


    output_filename = year + '-' + cruise_number + '_CTD_DATA.csv'

    ctd_data.to_csv(dest_dir + output_filename, index = False)

#MERGE_FILES(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2021-038/', year = '2021', cruise_number='038')

#----------------------------   Step 2. Create Metadata dictionay    ---------------------------------------------------
def CREATE_META_DICT(dest_dir, file, year, cruise_number):
    """
     Read in a csv file and output a metadata dictionary
     Inputs:
         - folder, file, year, cruise_number: rsk-format file containing raw RBR data & csv file containing metadata
     Outputs:
         - metadata dictionary
     """
    meta_dict = {}
    # function to export .csv files from .rsk file.
    rsk_filename = str(dest_dir) + str(file) #full path and name of a .rsk file
    rsk = pyrsktools.open(rsk_filename) # load up an RSK (only one even if many - IF settings and metadata are teh

    header_input_name = str(year) + '-' + str(cruise_number) + '_header-merge.csv'
    header_input_filename = dest_dir + header_input_name
    header = pd.read_csv(header_input_filename, header=0)

    # get the time interval for the IOS Header (sampling period)
    time_input_name = str(year) + '-' + str(cruise_number) + '_CTD_DATA.csv'
    time_input_filename = dest_dir + time_input_name
    time_input = pd.read_csv(time_input_filename)
    time_interval = pd.to_datetime(time_input['Time(yyyy-mm-dd HH:MM:ss.FFF)'][2]) - pd.to_datetime(
        time_input['Time(yyyy-mm-dd HH:MM:ss.FFF)'][1])
    time_interval = str(time_interval)
    time_interval = time_interval[-8:-3]

    csv_input_name = str(year) + '-' + str(cruise_number) + '_METADATA.csv'
    csv_input_filename = dest_dir + csv_input_name

    meta_csv = pd.read_csv(csv_input_filename)

    #meta_dict['number_of_profiles'] = len(list(rsk.profiles()))
    meta_dict['number_of_profiles'] = meta_csv['Value'][meta_csv['Name'] == 'number_of_profiles'].values[0]         # CHOOSE !!! if multiple rsk
    meta_dict['Processing_Start_time'] = datetime.now()
    meta_dict['Instrument_information'] = rsk.instrument
    meta_dict['Sampling_Interval'] = time_interval
    #meta_dict['RSK_filename'] = rsk.name
    meta_dict['RSK_filename'] = meta_csv['Value'][meta_csv['Name'] == 'RSK_filename'].values[0:] # if more than one (list??)
    meta_dict['Channels'] = list(rsk.channels.keys())
    meta_dict['Channel_details'] = list(rsk.channels.items())
    meta_dict['Data_description'] = meta_csv['Value'][meta_csv['Name'] == 'Data_description'].values[0]
    meta_dict['Final_file_type'] = meta_csv['Value'][meta_csv['Name'] == 'Final_file_type'].values[0]
    meta_dict['Number_of_channels'] = meta_csv['Value'][meta_csv['Name'] == 'Number_of_channels'].values[0]
    #meta_dict['Number_of_channels'] = len(list(rsk.channels.items())
    meta_dict['Mission'] = meta_csv['Value'][meta_csv['Name'] == 'Mission'].values[0]
    meta_dict['Agency'] = meta_csv['Value'][meta_csv['Name'] == 'Agency'].values[0]
    meta_dict['Country'] = meta_csv['Value'][meta_csv['Name'] == 'Country'].values[0]
    meta_dict['Project'] = meta_csv['Value'][meta_csv['Name'] == 'Project'].values[0]
    meta_dict['Scientist'] = meta_csv['Value'][meta_csv['Name'] == 'Scientist'].values[0]
    meta_dict['Platform'] = meta_csv['Value'][meta_csv['Name'] == 'Platform'].values[0]
    meta_dict['Instrument_Model'] = meta_csv['Value'][meta_csv['Name'] == 'Instrument_Model'].values[0]
    meta_dict['Serial_number'] = meta_csv['Value'][meta_csv['Name'] == 'Serial_number'].values[0]
    meta_dict['Instrument_type'] = meta_csv['Value'][meta_csv['Name'] == 'Instrument_type'].values[0]
    meta_dict['Location'] = header
    return meta_dict

#metadata = CREATE_META_DICT(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2021-070/', file = '204694_20211018_1808.rsk', year = '2021', cruise_number= '070')


#-------------------------------  step 3. Add 6 line headers to CTD_DATA.csv file--------------------------------------------------
# Prepare data file with six line header for further applications in IOS Shell
#ADD_6LINEHEADER(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2019-107-example/', year = '2019', cruise_number = '107')

#### NOT USING THIS FUNCTION ################
def ADD_6LINEHEADER(dest_dir, year, cruise_number):  # Original by Lu, adapted below to account for differing cols
    """
     Read in a csv file and output in csv format for IOSShell
     Inputs:
         - folder, file, year, cruise: csv-format file containing raw RBR CTD data exported from rsk file
     Outputs:
         - csv file: csv files containing 6-header line for IOSShell
     """
    # Add six-line header to the .csv file.
    # This file could be used for data processing via IOSShell
    input_name = str(year) + '-' + str(cruise_number) + '_CTD_DATA.csv'
    output_name = str(year) + "-" + str(cruise_number) + '_CTD_DATA-6linehdr.csv'
    input_filename = dest_dir + input_name
    ctd_data = pd.read_csv(input_filename, header=0)


    ctd_data['Time(yyyy-mm-dd HH:MM:ss.FFF)'] = ctd_data['Time(yyyy-mm-dd HH:MM:ss.FFF)'].str[:19]
    ctd_data['Date'] = pd.to_datetime(ctd_data['Time(yyyy-mm-dd HH:MM:ss.FFF)']) #add new column of Date
    ctd_data['Date'] = [d.date() for d in ctd_data['Date']]
    ctd_data['TIME:UTC'] = pd.to_datetime(ctd_data['Time(yyyy-mm-dd HH:MM:ss.FFF)']) #add new column of time
    ctd_data['TIME:UTC'] = [d.time() for d in ctd_data['TIME:UTC']]
    ctd_data['Date'] = pd.to_datetime(ctd_data['Date'], format='%Y-%m-%d').dt.strftime('%d/%m/%Y')



    ctd_data = ctd_data.drop(ctd_data.columns[[0,1, 6, 13, 14]], 1) #  This Changes!!!
    #cols = list(ctd_data.columns)
    #cols = cols[0:4] + cols[5:]
    #ctd_data=ctd_data[cols]
    columns = ctd_data.columns.tolist()


    column_names = dict.fromkeys(ctd_data.columns, '') #set empty column names
    ctd_data = ctd_data.rename(columns=column_names)   #remove column names

    #write header information into a dataframe
    channel = ['Y', 'Y', 'N', 'Y', 'Y', 'Y', 'Y', 'Y', 'N', 'N', 'N', 'Y', 'Y','Y']
    index = ['Conductivity', 'Temperature', 'Pressure_Air', 'Fluorescence', 'Oxygen:Dissolved:Saturation', 'Pressure', 'Depth', 'Salinity:CTD', ' ', ' ', 'Cast_direction', 'Event_number', 'Date', 'TIME:UTC']
    unit = ['mS/cm', 'deg C (ITS90)', 'decibar', 'mg/m^3', '%', 'decibar', 'Metres', 'PSS-78', ' ', ' ', 'n/a', 'n/a', 'n/a', 'n/a']
    input_format = ['R4', 'R4', 'R4', 'R4', 'R4', 'R4', 'R4', 'R4', ' ', ' ', ' ', 'I4', 'D:dd/mm/YYYY', 'T:HH:MM:SS']
    output_format = ['R4:F11.4', 'R4:F9.4', 'R4:F7.1', 'R4:F8.3', 'R4:F11.4', 'R4:F7.1', 'R4:F7.1', 'R4:F9.4', ' ', ' ', ' ', 'I:I4', 'D:YYYY/mm/dd', 'T:HH:MM:SS']
    na_value = ['-99', '-99', '-99', '-99', '-99', '-99', '-99', '-99', '-99', '', '', '', '', '']
    header = pd.DataFrame([channel, index, unit, input_format, output_format, na_value])
    column_names_header = dict.fromkeys(header.columns, '') #set empty column names
    header = header.rename(columns=column_names_header)

    ctd_data_header = header.append(ctd_data)
    ctd_data_header.to_csv(dest_dir + output_name, index = False, header = False)


#ADD_6LINEHEADER(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2020-085/', year = '2020', cruise_number = '085')


#____________________________________ trying to loop through column names to create lists of channels for header

# Prepare data file with six line header for further applications in IOS Shell
#ADD_6LINEHEADER(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2019-107-example/', year = '2019', cruise_number = '107')

#### USING THIS FUNCTION ###########
def ADD_6LINEHEADER_2(dest_dir, year, cruise_number):
    """
     Read in a csv file and output in csv format for IOSShell
     Filter through the csv and remove un-needed columns.
     Inputs:
         - folder, file, year, cruise: csv-format file containing raw RBR CTD data exported from rsk file
     Outputs:
         - csv file: csv files containing 6-header line for IOSShell
     """
    # Add six-line header to the .csv file.
    # This file could be used for data processing via IOSShell
    input_name = str(year) + '-' + str(cruise_number) + '_CTD_DATA.csv'
    output_name = str(year) + "-" + str(cruise_number) + '_CTD_DATA-6linehdr.csv'
    input_filename = dest_dir + input_name
    ctd_data = pd.read_csv(input_filename, header=0)


    ctd_data['Time(yyyy-mm-dd HH:MM:ss.FFF)'] = ctd_data['Time(yyyy-mm-dd HH:MM:ss.FFF)'].str[:19]
    ctd_data['Date'] = pd.to_datetime(ctd_data['Time(yyyy-mm-dd HH:MM:ss.FFF)']) #add new column of Date
    ctd_data['Date'] = [d.date() for d in ctd_data['Date']]
    ctd_data['TIME:UTC'] = pd.to_datetime(ctd_data['Time(yyyy-mm-dd HH:MM:ss.FFF)']) #add new column of time
    ctd_data['TIME:UTC'] = [d.time() for d in ctd_data['TIME:UTC']]
    ctd_data['Date'] = pd.to_datetime(ctd_data['Date'], format='%Y-%m-%d').dt.strftime('%d/%m/%Y')

    # make a list of columns to drop
    drop_list = ['Temperature .1(Â°C)', 'speedofsound(m/s)', 'Temperature .1', 'Temperature.1', 'Speed of sound ', 'Speed of sound',
                'specificconductivity(ÂµS/cm)', 'Density anomaly', 'Density anomaly ', 'Dissolved Oâ,,Concentration',
                 'Dissolved OÃ¢Â‚Â‚ concentration ', 'Dissolved OÃ¢Â‚Â‚ concentration', 'Specific conductivity ', 'Dissolved O₂ concentration',
                 'Specific conductivity', 'Dissolved Oâ\x82\x82 concentration ', 'Dissolved Oâ\x82\x82 concentration']


    #drop first indexing row and Time(HH....)
    ctd_data = ctd_data.drop(ctd_data.columns[[0, 1]], 1)

    # Create empty lists
    channel = []
    index = []
    unit = []
    input_format = []
    output_format = []
    na_value = []
    #drop columns in the drop list
    ctd_data = ctd_data[ctd_data.columns[~ctd_data.columns.isin(drop_list)]]
    #ctd_data.reset_index(drop=True, inplace=True)

    col_list = ctd_data.columns.tolist()
    print(col_list)
    n_col_list = len(col_list)

    column_names = dict.fromkeys(ctd_data.columns, '') #set empty column names
    ctd_data = ctd_data.rename(columns=column_names)   #remove column names

    #append header information into the empty lists
    channel_list = ['Y', 'Y', 'N', 'Y', 'Y', 'Y', 'Y', 'Y', 'N', 'Y', 'Y', 'Y']
    index_list = ['Conductivity', 'Temperature', 'Pressure_Air', 'Fluorescence', 'Oxygen:Dissolved:Saturation', 'Pressure', 'Depth', 'Salinity:CTD', 'Cast_direction', 'Event_number', 'Date', 'TIME:UTC']
    unit_list = ['mS/cm', 'deg C(ITS90)', 'decibar', 'mg/m^3', '%', 'decibar', 'meters', 'PSS-78', 'n/a', 'n/a', 'n/a', 'n/a']
    input_format_list = ['R4', 'R4', 'R4', 'R4', 'R4', 'R4', 'R4', 'R4', ' ', 'I4', 'D:dd/mm/YYYY', 'T:HH:MM:SS']
    output_format_list = ['R4:F11.4', 'R4:F9.4', 'R4:F7.1', 'R4:F8.3', 'R4:F11.4', 'R4:F7.1', 'R4:F7.1', 'R4:F9.4', ' ', 'I:I4', 'D:YYYY/mm/dd', 'T:HH:MM:SS']
    na_value_list = ['-99', '-99', '-99', '-99', '-99', '-99', '-99', '-99', '','-99', '', '']

    for col in col_list:
        if col == 'conductivity(mS/cm)' or col == 'Conductivity ' or col == 'Conductivity':
            channel.append(channel_list[0]), index.append(index_list[0]), unit.append(unit_list[0]), input_format.append(input_format_list[0]), output_format.append(output_format_list[0]), na_value.append(na_value_list[0])
        elif col == 'temperature(Â°C)' or col == 'Temperature ' or col == 'Temperature':
            channel.append(channel_list[1]), index.append(index_list[1]), unit.append(unit_list[1]), input_format.append(input_format_list[1]), output_format.append(output_format_list[1]), na_value.append(na_value_list[1])
        elif col == 'pressure(dbar)' or col == 'Pressure ' or col == 'Pressure':
            channel.append(channel_list[2]), index.append(index_list[2]), unit.append(unit_list[2]), input_format.append(input_format_list[2]), output_format.append(output_format_list[2]), na_value.append(na_value_list[2])
        elif col == 'chlorophyll(Âµg/l)' or col == 'Chlorophyll a ' or col == 'Chlorophyll a':
            channel.append(channel_list[3]), index.append(index_list[3]), unit.append(unit_list[3]), input_format.append(input_format_list[3]), output_format.append(output_format_list[3]), na_value.append(na_value_list[3])
        elif col == 'oxygensaturation(%)' or col == 'Dissolved OÃ¢Â‚Â‚ saturation ' or col == 'Dissolved Oâ\x82\x82 saturation ' or col == 'Dissolved OÃ¢Â‚Â‚ saturation' or col == 'Dissolved Oâ\x82\x82 saturation':
            channel.append(channel_list[4]), index.append(index_list[4]), unit.append(unit_list[4]), input_format.append(input_format_list[4]), output_format.append(output_format_list[4]), na_value.append(na_value_list[4])
        elif col == 'seapressure(dbar)'  or col == 'Sea pressure' or col == 'Sea pressure ':
            channel.append(channel_list[5]), index.append(index_list[5]), unit.append(unit_list[5]), input_format.append(input_format_list[5]), output_format.append(output_format_list[5]), na_value.append(na_value_list[5])
        elif col == 'depth(m)' or col == 'Depth ' or col == 'Depth':
            channel.append(channel_list[6]), index.append(index_list[6]), unit.append(unit_list[6]), input_format.append(input_format_list[6]), output_format.append(output_format_list[6]), na_value.append(na_value_list[6])
        elif col == 'salinity(PSU)' or col == 'Salinity ' or col == 'Salinity':
            channel.append(channel_list[7]), index.append(index_list[7]), unit.append(unit_list[7]), input_format.append(input_format_list[7]), output_format.append(output_format_list[7]), na_value.append(na_value_list[7])
        elif col == 'Cast_direction':
            channel.append(channel_list[8]), index.append(index_list[8]), unit.append(unit_list[8]), input_format.append(input_format_list[8]), output_format.append(output_format_list[8]), na_value.append(na_value_list[8])
        elif col == 'Event':
            channel.append(channel_list[9]), index.append(index_list[9]), unit.append(unit_list[9]), input_format.append(input_format_list[9]), output_format.append(output_format_list[9]), na_value.append(na_value_list[9])
        elif col == 'Date':
            channel.append(channel_list[10]), index.append(index_list[10]), unit.append(unit_list[10]), input_format.append(input_format_list[10]), output_format.append(output_format_list[10]), na_value.append(na_value_list[10])
        elif col == 'TIME:UTC':
            channel.append(channel_list[11]), index.append(index_list[11]), unit.append(unit_list[11]), input_format.append(input_format_list[11]), output_format.append(output_format_list[11]), na_value.append(na_value_list[11])


    header = pd.DataFrame([channel, index, unit, input_format, output_format, na_value])
    #print(header[4])
    #print(ctd_data)
    column_names_header = dict.fromkeys(header.columns, '') #set empty column names
    header = header.rename(columns=column_names_header)

    ctd_data_header = header.append(ctd_data)
    ctd_data_header.to_csv(dest_dir + output_name, index = False, header = False)

#ADD_6LINEHEADER_2(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2021-038/', year = '2021', cruise_number = '038')



#define function to plot and check the location of CTD

#left_lon, right_lon, bot_lat, top_lat = [-131, -132, 53.5, 53.9]

#run function to plot the locations of RBR casts
#Plot_Track_Location(dest_dir='C:/Projects/RBR_CTD_Pycharm/2020-085/', year = '2020', cruise_number= '085', left_lon=-131.5, right_lon=-131.8, bot_lat=52.5, top_lat=52.7)


def Plot_Track_Location(dest_dir, year, cruise_number, left_lon, right_lon, bot_lat, top_lat):
    """
     Read in a csv file and output a map
     Inputs:
         - folder, year, cruise: csv file containing raw RBR data
     Outputs:
         - A map showing sampling locations
     """
    input_name = str(year) + '-' + str(cruise_number) + '_header-merge.csv'
    input_filename = dest_dir + input_name
    header = pd.read_csv(input_filename, header=0)
    header['lat_degree'] = header['LOC:LATITUDE'].str[:2].astype(int)
    header['lat_min'] = header['LOC:LATITUDE'].str[3:10].astype(float)
    header['lat'] = header['lat_degree'] + header['lat_min']/60
    header['lon_degree'] = header['LOC:LONGITUDE'].str[:3].astype(int)
    header['lon_min'] = header['LOC:LONGITUDE'].str[4:12].astype(float)
    header['lon'] = 0 - (header['lon_degree'] + header['lon_min']/60)
    event = header['LOC:STATION'].astype(str)

    lon = header['lon'].tolist()
    lat = header['lat'].tolist()
    event = event.tolist()

    m = Basemap(llcrnrlon=left_lon, llcrnrlat=bot_lat,
                urcrnrlon=right_lon, urcrnrlat=top_lat,
                projection='lcc',
                resolution='h', lat_0=0.5 * (bot_lat + top_lat),
                lon_0=0.5 * (left_lon + right_lon))  # lat_0=53.4, lon_0=-129.0)

    x, y = m(lon, lat)

    fig = plt.figure(num=None, figsize=(8, 6), dpi=100)
    m.drawcoastlines(linewidth=0.2)
    m.drawmapboundary(fill_color='white')
    #m.fillcontinents(color='0.8')
    m.drawrivers()



    m.scatter(x, y, marker='D', color='m', s=5)
    #m.plot(x, y, marker='D', color='m', markersize=4)
 #   for event, xpt, ypt in zip(event, x, y):
 #       plt.text(xpt, ypt, event)

    parallels = np.arange(bot_lat, top_lat, 0.5)  # parallels = np.arange(48., 54, 0.2), parallels = np.linspace(bot_lat, top_lat, 10)
    m.drawparallels(parallels, labels=[True, False, True, False])  # draw parallel lat lines
    meridians = np.arange(left_lon, right_lon, 0.5)
    m.drawmeridians(meridians, labels=[False, False, False, True])
    plt.title(year + '-' + cruise_number)
    plt.savefig(dest_dir + 'Fig_1.png')
    plt.show()

#---------------------  create a second map w/ Cartopy # just to double check - had an issue with Basemap once only

    map = plt.axes(projection=ccrs.PlateCarree())
    map.set_extent([left_lon, right_lon, bot_lat, top_lat]) # try left_lon, right_lon, bot_lat, top_lat
    x, y = (lon, lat)
    map.coastlines()
    gl = map.gridlines(crs=ccrs.PlateCarree(), linewidth = 0.5, color = 'black', alpha =0.5, linestyle='--', draw_labels =True)
    gl.top_labels = False
    gl.left_labels = True
    gl.bottom_labels = True
    gl.right_labels = False
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()

    gl.xlabel_style = {'color': 'black', 'weight': 'bold', 'size': 6}
    gl.ylabel_style = {'color': 'black', 'weight': 'bold', 'size': 6}


    cax = plt.scatter(x, y, transform=ccrs.PlateCarree(), marker='.', color='red', s=25)
    plt.title(year + '-' + cruise_number)
    plt.savefig(dest_dir + 'Figure_1.png')
    plt.show()

#Plot_Track_Location(dest_dir='C:/Projects/RBR_CTD_Pycharm/2021-038/', year = '2021', cruise_number= '038', left_lon=-115, right_lon=-126, bot_lat=67, top_lat=73)

#------------------------------------------------- Step 5.  create variable dictionaries  ------------------------------------------------------
#input: .csv file

#output variables: cast cast_d, cast_u



#---------------------------   step 6. Plot and Check and correct for zero-order hold   --------------------------------------------------

#PLOT_PRESSURE_DIFF(dest_dir = '/home/guanl/Desktop/Projects/RBR/Processing/2020-085/CTD data/', year= '2020', cruise_number= '085')

def PLOT_PRESSURE_DIFF(dest_dir, year, cruise_number, input_ext):
    """
     Read in a csv file and output a plot to check zero-order holds
     Inputs:
         - folder, year, cruise: csv file containing raw RBR data
     Outputs:
         - a plot showing the time derivative of raw pressure
     """

    input_name = str(year) + "-" + str(cruise_number) + input_ext # will depend on the need for zero order holds
    input_filename = dest_dir + input_name
    ctd_data = pd.read_csv(input_filename, header=None, low_memory=False)
    ctd_data = ctd_data.rename(columns=ctd_data.iloc[1]) #assign the second row as column names
    ctd_data = ctd_data.rename(columns={'Oxygen:Dissolved:Saturation': 'Oxygen', 'Salinity:CTD': 'Salinity', 'TIME:UTC': 'TIME'})
    ctd = ctd_data.iloc[6:]
    ctd.index = np.arange(0, len(ctd))
    #ctd = ctd[1000:4000] # to limit the number of records plotted -

    pressure = ctd['Pressure'].apply(pd.to_numeric)
    pressure_lag = pressure[1:]
    pressure_lag.index = np.arange(0, len(pressure_lag))
    pressure_diff = pressure_lag - pressure


    fig = plt.figure(num=None, figsize=(14, 6), dpi=100)
    plt.plot(pressure_diff, color='blue', linewidth= 0.5, label='Pressure_diff')
    plt.ylabel('Pressure (decibar)')
    plt.xlabel('Scans')
    plt.grid()
    plt.legend()
    plt.title(year + '-' + cruise_number + ' ' + input_ext)
    plt.savefig(dest_dir + 'zero_order_holds_' + input_ext + '.png')
    plt.show()
    #check the plot then save it as Fig_2 or Fig_3.

#PLOT_PRESSURE_DIFF(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2021-038/', year= '2021', cruise_number= '038', input_ext = '_CTD_DATA-6linehdr.csv')


# output variables: cast cast_d, cast_u

def CREATE_CAST_VARIABLES(year, cruise_number, dest_dir, input_ext):
    """
     Read in a csv file and output data dictionaries to hold profile data
     Inputs:
         - folder, year, cruise: csv file containing raw RBR data
     Outputs:
         - three dictionaries containing casts, downcasts and upcasts
     """
    input_name = str(year) + "-" + str(cruise_number) + input_ext
    input_filename = dest_dir + input_name
    ctd_data = pd.read_csv(input_filename, header=None, low_memory=False)  # read data without header
    ctd_data = ctd_data.rename(columns=ctd_data.iloc[1])  # assign the second row as column names
    ctd_data = ctd_data.rename(
        columns={'Oxygen:Dissolved:Saturation': 'Oxygen', 'Salinity:CTD': 'Salinity', 'TIME:UTC': 'TIME'})
    ctd = ctd_data.iloc[6:]
    # drop NaNs from Zero order holds correction (not including O or F in case they aren't present - but will capture)
    if input_ext == '_CTD_DATA-6linehdr.csv':
        ctd = ctd
    elif input_ext == '_CTD_DATA-6linehdr_corr_hold.csv':
        ctd = ctd = ctd.dropna(0,
                               subset=['Conductivity', 'Temperature', 'Pressure_Air', 'Pressure', 'Depth', 'Salinity'],
                               how='all')## I don't think this does anything - NaNs now dropped in Correct_Hold stage
    ctd = ctd.copy()
    cols = ctd.columns[
           0:-4]
    cols_names = ctd.columns.tolist()
    ctd[cols] = ctd[cols].apply(pd.to_numeric, errors='coerce', axis=1)
    ctd['Cast_direction'] = ctd['Cast_direction'].str.strip()


    n = ctd['Event_number'].nunique()

    var_holder = {}
    for i in range(1, n + 1, 1):
        var_holder['cast' + str(i)] = ctd.loc[(ctd['Event_number'] == str(i))]
    # var_holder['Processing_history'] = ""

    var_holder_d = {}
    for i in range(1, n + 1, 1):
        var_holder_d['cast' + str(i)] = ctd.loc[(ctd['Event_number'] == str(i)) & (ctd['Cast_direction'] == 'd')]
    # var_holder_d['Processing_history'] = ""

    var_holder_u = {}
    for i in range(1, n + 1, 1):
        var_holder_u['cast' + str(i)] = ctd.loc[(ctd['Event_number'] == str(i)) & (ctd['Cast_direction'] == 'u')]
    # var_holder_u['Processing_history'] = ""

    return var_holder, var_holder_d, var_holder_u


#cast, cast_d, cast_u = CREATE_CAST_VARIABLES(year = '2021' , cruise_number = '045', dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2021-045/IOS_SHELL_STEPS', input_ext = '_CTD_DATA-6linehdr.csv')

# --------------------------------------   plot data from all profiles  -----------------------------------------------------
# plot salinity of all cast together

def First_Plots(year, cruise_number, dest_dir, input_ext):


    """ Plot pre-processing and after Zero-order Holds if needed"""

    cast, cast_d, cast_u = CREATE_CAST_VARIABLES(year, cruise_number,
                                                 dest_dir, input_ext)

    number_of_colors = len(cast)
    color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                    for i in range(number_of_colors)]

    #get variables
    vars = list(dict.fromkeys(cast['cast1']))


    fig, ax = plt.subplots()
    for i in range (0, len(cast), 1):
        ax.plot(cast_d['cast' + str(i+1)].Salinity, cast_d['cast' + str(i+1)].Pressure, color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_u['cast' + str(i+1)].Salinity, cast_u['cast' + str(i+1)].Pressure, '--', color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_d['cast1'].Salinity, cast_d['cast1'].Pressure, color='blue', label='cast1')
        #ax.plot(cast_u['cast1'].Salinity, cast_u['cast1'].Pressure, '--', color='blue', label='cast1')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('Pre-Processing', fontsize=5)
    plt.savefig(dest_dir + 'Pre_Processing_S.png')
    ax.legend()



            #plot temperature of all cast together
    fig, ax = plt.subplots()
    for i in range (0, len(cast), 1):
        ax.plot(cast_d['cast' + str(i+1)].Temperature, cast_d['cast' + str(i+1)].Pressure, color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_u['cast' + str(i+1)].Temperature, cast_u['cast' + str(i+1)].Pressure, '--', color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_d['cast1'].Temperature, cast_d['cast1'].Pressure, color='blue', label='cast1')
    #ax.plot(cast_u['cast1'].Temperature, cast_u['cast1'].Pressure, '--', color='blue', label='cast1')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Temperature(C)')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('Pre-Processing', fontsize=5)
    plt.savefig(dest_dir + 'Pre_Processing_T.png')
    ax.legend()


        #plot Conductivity of all cast together
    fig, ax = plt.subplots()
    for i in range (0, len(cast), 1):
        ax.plot(cast_d['cast' + str(i+1)].Conductivity, cast_d['cast' + str(i+1)].Pressure, color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_u['cast' + str(i+1)].COnductivity, cast_u['cast' + str(i+1)].Pressure, '--', color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_d['cast1'].Conductivity, cast_d['cast1'].Pressure, color='blue', label='cast1')
        #ax.plot(cast_d['cast2'].Conductivity[1:], cast_d['cast2'].Pressure[1:], color='red', label='cast2')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Conductivity (S/cm)')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('Pre-Processing', fontsize=5)
    plt.savefig(dest_dir + 'Pre_Processing_C.png')
    ax.legend()

    for var in vars:
        if var == 'Oxygen':
            #plot Oxygen of all cast together
            fig, ax = plt.subplots()
            for i in range (0, len(cast), 1):
                ax.plot(cast_d['cast' + str(i+1)].Oxygen, cast_d['cast' + str(i+1)].Pressure, color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_u['cast' + str(i+1)].Oxygen, cast_u['cast' + str(i+1)].Pressure, '--', color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_d['cast1'].Oxygen, cast_d['cast1'].Pressure, color='blue', label='cast1')
        #ax.plot(cast_u['cast1'].Oxygen, cast_u['cast1'].Pressure, '--', color='blue', label='cast1')
            ax.invert_yaxis()
            ax.xaxis.set_label_position('top')
            ax.xaxis.set_ticks_position('top')
            ax.set_xlabel('Oxygen Saturation (%)')   # Check unit here
            ax.set_ylabel('Pressure (decibar)')
            ax.set_title('Pre-Processing', fontsize=5)
            plt.savefig(dest_dir + 'Pre_Processing_O.png')
            ax.legend()

    #plot Fluorescence of all cast together
    for var in vars:
        if var == 'Fluorescence':
            fig, ax = plt.subplots()
            for i in range (0, len(cast), 1):
                ax.plot(cast_d['cast' + str(i+1)].Fluorescence, cast_d['cast' + str(i+1)].Pressure, color=color[i], label= 'cast' + str(i+1))
            #ax.plot(cast_u['cast' + str(i+1)].Fluorescence, cast_u['cast' + str(i+1)].Pressure, '--', color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_d['cast1'].Fluorescence, cast_d['cast1'].Pressure, color='blue', label='cast1')
        #ax.plot(cast_u['cast1'].Fluorescence, cast_u['cast1'].Pressure, '--', color='blue', label='cast1')
            ax.invert_yaxis()
            ax.xaxis.set_label_position('top')
            ax.xaxis.set_ticks_position('top')
            ax.set_xlabel('Fluorescence (ug/L)')   # Check unit here
            ax.set_ylabel('Pressure (decibar)')
            ax.set_title('Pre-Processing', fontsize=5)
            plt.savefig(dest_dir + 'Pre_Processing_F.png')
            ax.legend()


        # TS Plot
    fig, ax = plt.subplots()
    for i in range (0, len(cast), 1):
        ax.plot(cast_d['cast' + str(i+1)].Salinity, cast_d['cast' + str(i+1)].Temperature, color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_d['cast11'].Salinity, cast_d['cast11'].Temperature, color='blue')
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Temperature (C)')
    ax.set_title('Pre-Processing T-S Plot')
    plt.savefig(dest_dir + 'Pre_Processing_T-S.png')
    ax.legend()


        #--------------------------------------  Plot profiles by group ------------------------------------------------------------------------
        # T, C, O, S, F for each profile in one plot
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, sharey=True)
        #Temperature
    ax1.plot(cast_d['cast1'].Temperature, cast_d['cast1'].Pressure, color='red', label='cast_down')
    ax1.plot(cast_u['cast1'].Temperature, cast_u['cast1'].Pressure, '--', color='red', label='cast_up')
    ax1.set_ylabel('Pressure(decibar)', fontsize=8)
    ax1.set_ylim(ax1.get_ylim()[::-1])
    ax1.set_xlabel('Temperature(C)', fontsize=8)
    ax1.xaxis.set_label_position('top')
    ax1.xaxis.set_ticks_position('top')
    ax.set_title('Pre-Processing', fontsize=5)
    ax1.legend()



        #Conductivity
    ax2.plot(cast_d['cast1'].Conductivity, cast_d['cast1'].Pressure, color='yellow', label='cast_down')
    ax2.plot(cast_u['cast1'].Conductivity, cast_u['cast1'].Pressure, '--', color='yellow', label='cast1_up')
    ax2.set_ylabel('Pressure(decibar)', fontsize=8)
    ax2.set_ylim(ax1.get_ylim()[::-1])
    ax2.set_xlabel('Conductivity (S/cm)', fontsize=8)
    ax2.xaxis.set_label_position('top')
    ax2.xaxis.set_ticks_position('top')
    ax.set_title('Pre-Processing', fontsize=5)
    ax2.legend()


        #Oxygen
    for var in vars:
        if var == 'Oxygen':
            ax3.plot(cast_d['cast1'].Oxygen, cast_d['cast1'].Pressure, color='black', label='cast_down')
            ax3.plot(cast_u['cast1'].Oxygen, cast_u['cast1'].Pressure, '--', color='black', label='cast_up')
            ax3.set_ylabel('Pressure(decibar)', fontsize=8)
            ax3.set_ylim(ax1.get_ylim()[::-1])
            ax3.set_xlabel('Oxygen Saturation (%)', fontsize=8)
            ax3.xaxis.set_label_position('top')
            ax3.xaxis.set_ticks_position('top')
            ax.set_title('Pre-Processing', fontsize=5)
            ax3.legend()

        #Salinity
    ax4.plot(cast_d['cast1'].Salinity, cast_d['cast1'].Pressure, color='blue', label='cast_down')
    ax4.plot(cast_u['cast1'].Salinity, cast_u['cast1'].Pressure, '--', color='blue', label='cast_up')
    ax4.set_ylabel('Pressure(decibar)', fontsize=8)
    ax4.set_ylim(ax1.get_ylim()[::-1])
    ax4.set_xlabel('Salinity', fontsize=8)
    ax4.xaxis.set_label_position('top')
    ax4.xaxis.set_ticks_position('top')
    ax.set_title('Pre-Processing', fontsize=5)
    ax4.legend()

        #Fluorescence
    for var in vars:
        if var == 'Fluorescence':
            ax5.plot(cast_d['cast1'].Fluorescence, cast_d['cast1'].Pressure, color='green', label='cast_down')
            ax5.plot(cast_u['cast1'].Fluorescence, cast_u['cast1'].Pressure, '--', color='green', label='cast1_up')
            ax5.set_ylabel('Pressure(decibar)', fontsize=8)
            ax5.set_ylim(ax1.get_ylim()[::-1])
            ax5.set_xlabel('Fluoresence(ug/L)', fontsize=8)
            ax5.xaxis.set_label_position('top')
            ax5.xaxis.set_ticks_position('top')
            ax.set_title('Pre-Processing', fontsize=5)
            ax5.legend()


        #--------------------------------   Plot by Index and by Profile--------------------------------------------------------------------------------------------
        # separate plot for T, C, O, S, F of each profile
    #Temperature
    fig, ax = plt.subplots()
    ax.plot(cast_d['cast1'].Temperature, cast_d['cast1'].Pressure, color='red', label='cast_down')
    ax.plot(cast_u['cast1'].Temperature, cast_u['cast1'].Pressure, '--', color='red', label='cast_up')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Temperature(C)')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('Pre-Processing', fontsize=5)
    plt.savefig(dest_dir + 'Pre_Cast1_T')
    ax.legend()

        #Salinity
    fig, ax = plt.subplots()
    ax.plot(cast_d['cast1'].Salinity, cast_d['cast1'].Pressure, color='blue', label='cast1_d')
    ax.plot(cast_u['cast1'].Salinity, cast_u['cast1'].Pressure, '--', color='blue', label='cast1_u')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('Pre-Processing', fontsize=5)
    plt.savefig(dest_dir + 'Pre_Cast1_S')
    ax.legend()

        #Conductivity
    fig, ax = plt.subplots()
    ax.plot(cast_d['cast1'].Conductivity, cast_d['cast1'].Pressure, color='yellow', label='cast1_d')
    ax.plot(cast_u['cast1'].Conductivity, cast_u['cast1'].Pressure, '--', color='yellow', label='cast1_u')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Conductivity (S/cm)')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('Pre-Processing', fontsize=5)
    plt.savefig(dest_dir + 'Pre_Cast1_C')
    ax.legend()


        #Oxygen
    for var in vars:
        if var == 'Oxygen':
            fig, ax = plt.subplots()
            ax.plot(cast_d['cast1'].Oxygen, cast_d['cast1'].Pressure, color='black', label='cast1_d')
            ax.plot(cast_u['cast1'].Oxygen, cast_u['cast1'].Pressure, '--', color='black', label='cast1_u')
            ax.invert_yaxis()
            ax.xaxis.set_label_position('top')
            ax.xaxis.set_ticks_position('top')
            ax.set_xlabel('Oxygen Saturation (%)')   # Check unit here
            ax.set_ylabel('Pressure (decibar)')
            ax.set_title('Pre-Processing', fontsize=5)
            plt.savefig(dest_dir + 'Pre_Cast1_O')
            ax.legend()

        #Fluoresence
    for var in vars:
        if var == 'Fluorescence':
            fig, ax = plt.subplots()
            ax.plot(cast_d['cast1'].Fluorescence, cast_d['cast1'].Pressure, color='green', label='cast1_d')
            ax.plot(cast_u['cast1'].Fluorescence, cast_u['cast1'].Pressure, '--', color='green', label='cast1_u')
            ax.invert_yaxis()
            ax.xaxis.set_label_position('top')
            ax.xaxis.set_ticks_position('top')
            ax.set_xlabel('Fluorescence (ug/L)')   # Check unit here
            ax.set_ylabel('Pressure (decibar)')
            ax.set_title('Pre-Processing', fontsize=5)
            plt.savefig(dest_dir + 'Pre_Cast1_F')
            ax.legend()

    fig, ax = plt.subplots()
    ax.plot(cast_d['cast1'].Salinity, cast_d['cast1'].Temperature, color='red',
                label='cast1_d')
    ax.plot(cast_u['cast1'].Salinity, cast_u['cast1'].Temperature, '--', color='blue', label='cast1_u')
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Temperature (C)')
    ax.set_title('Pre-Processing T-S Plot')
    plt.savefig(dest_dir + 'Pre_Cast1_T-S.png')
    ax.legend()

    number_of_colors = len(cast)
    color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

# pressure check

    fig, ax = plt.subplots()
    for i in range(0, len(cast), 1):
        ax.plot(cast_d['cast' + str(i + 1)].Conductivity[0:20], cast_d['cast' + str(i + 1)].Pressure[0:20],
                color=color[i], label='cast' + str(i + 1))
        ax.plot(cast_u['cast' + str(i + 1)].Conductivity[-20:-1], cast_u['cast' + str(i + 1)].Pressure[-20:-1], '--',
                color=color[i], label='cast' + str(i + 1))
    # ax.plot(cast_d['cast1'].Conductivity[0:10], cast_d['cast1'].Pressure[0:10], color='blue', label='cast1')
    # ax.plot(cast_d['cast2'].Conductivity[0:10], cast_d['cast2'].Pressure[0:10], color='red', label='cast2')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Conductivity (S/cm)')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('Checking need for Pressure correction')
    plt.savefig(dest_dir + 'PC_need_CvP.png')
    ax.legend()

    fig, ax = plt.subplots()
    for i in range(0, len(cast), 1):
        ax.plot(cast_d['cast' + str(i + 1)].Conductivity[0:20], cast_d['cast' + str(i + 1)].Depth[0:20], color=color[i],
                label='cast' + str(i + 1))
        ax.plot(cast_u['cast' + str(i + 1)].Conductivity[-20:-1], cast_u['cast' + str(i + 1)].Depth[-20:-1], '--', color=color[i],
                label='cast' + str(i + 1))
    # ax.plot(cast_d['cast1'].Conductivity[0:10], cast_d['cast1'].Depth[0:10], color='blue', label='cast1')
    # ax.plot(cast_d['cast2'].Conductivity[0:10], cast_d['cast2'].Depth[0:10], color='red', label='cast2')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Conductivity (S/cm)')
    ax.set_ylabel('Depth (m)')
    ax.set_title('Checking need for Pressure correction')
    plt.savefig(dest_dir + 'PC_need_CvD.png')
    ax.legend()
    plt.show()

##first_plots()

def CORRECT_HOLD(dest_dir, year, cruise_number, metadata_dict):
    """
    Read 6linehdr.csv and correct for zero order holds.  Look for repeat values in Pressure and replace with NaN, then
    look for repeats in the other sensors at the same place and replace
    with NaN..

    Adapted from RSKtools function RSKcorrecthold: 'This function identifies zero-hold points by looking
    for where consecutive differences for each channel are equal to zero, and replaces them with Nan or an
    interpolated value."  This function uses Nan. SH

    Output a new csv with the corrected values.
    """
    input_name = str(year) + "-" + str(cruise_number) + '_CTD_DATA-6linehdr.csv'
    output_name = str(year) + "-" + str(cruise_number) + '_CTD_DATA-6linehdr_corr_hold.csv'
    input_filename = dest_dir + input_name
    ctd_data = pd.read_csv(input_filename, header=None, low_memory=False)
    ctd_data = ctd_data.rename(columns=ctd_data.iloc[1]) #assign the second row as column names
    ctd_data = ctd_data.rename(columns={'Oxygen:Dissolved:Saturation': 'Oxygen', 'Salinity:CTD': 'Salinity', 'TIME:UTC': 'TIME'})
    header = ctd_data.iloc[0:6] # keep the header to use at the end
    ctd = ctd_data.iloc[6:]
    ctd = ctd.copy()
    vars = list(dict.fromkeys(ctd))
    #print(vars)
    #cols = ctd.columns[0:-4]
    #ctd[cols] = ctd[cols].apply(pd.to_numeric, errors='coerce', axis=1)
    ctd.index = np.arange(0, len(ctd))
    pressure = ctd['Pressure'].apply(pd.to_numeric)
    pressure_lag = pressure[1:]
    pressure_lag.index = np.arange(0, len(pressure_lag))
    air = ctd['Pressure_Air'].apply(pd.to_numeric)
    air_lag = air[1:]
    air_lag.index = np.arange(0, len(air_lag))
    conductivity = ctd['Conductivity'].apply(pd.to_numeric)
    conductivity_lag = conductivity[1:]
    conductivity_lag.index = np.arange(0, len(conductivity_lag))
    temperature = ctd['Temperature'].apply(pd.to_numeric)


    temperature_lag = temperature[1:]
    temperature_lag.index = np.arange(0, len(temperature_lag))
    for var in vars:
        if var == 'Fluorescence':
            fluorescence = ctd['Fluorescence'].apply(pd.to_numeric)
            fluorescence_lag = fluorescence[1:]
            fluorescence_lag.index = np.arange(0, len(fluorescence_lag))
    for var in vars:
        if var == 'Oxygen':
            oxygen = ctd['Oxygen'].apply(pd.to_numeric)
            oxygen_lag = oxygen[1:]
            oxygen_lag.index = np.arange(0, len(oxygen_lag))
    depth = ctd['Depth'].apply(pd.to_numeric)
    depth_lag = depth[1:]
    depth_lag.index = np.arange(0, len(depth_lag))
    salinity = ctd['Salinity'].apply(pd.to_numeric)
    salinity_lag = salinity[1:]
    salinity_lag.index = np.arange(0, len(salinity_lag))

    for i in range(0, len(ctd) - 1, 1):
        if pressure[i] == pressure_lag[i]:


            #pressure.iloc[i + 1] = np.nan
            if conductivity[i] == conductivity_lag[i]:
                conductivity.iloc[i + 1] = np.nan
            if air[i] == air_lag[i]:
                air.iloc[i + 1] = np.nan
            if temperature[i] == temperature_lag[i]:
                temperature.iloc[i + 1] = np.nan
            for var in vars:
                if var == 'Fluorescence':
                    if fluorescence[i] == fluorescence_lag[i]:
                        fluorescence.iloc[i + 1] = np.nan
            for var in vars:
                if var == 'Oxygen':
                    if oxygen[i] == oxygen_lag[i]:
                        oxygen.iloc[i + 1] = np.nan
            #if depth[i] == depth_lag[i]:
                #depth.iloc[i + 1] = np.nan
            if salinity[i] == salinity_lag[i]:
                salinity.iloc[i + 1] = np.nan


    #ctd['Pressure'] = pressure  # this worked when pressure was set to NaN
    ctd['Conductivity'] = conductivity
    ctd['Temperature'] = temperature
    ctd['Pressure_Air'] = air
    for var in vars:
        if var == 'Fluorescence':
            ctd['Fluorescence'] = fluorescence
    for var in vars:
        if var == 'Oxygen':
            ctd['Oxygen'] = oxygen
    #ctd['Depth'] = depth
    ctd['Salinity'] = salinity

    # drop the NaNs before they get into the CSV for IOS Shell Processing
    ctd = ctd.dropna(0, subset=['Conductivity', 'Temperature', 'Pressure_Air', 'Salinity'],
                     how='all') #'Pressure', 'Depth',  used to be in this list  #sometimes 'any' is required
    #ctd = ctd.reset_index(drop=True)

    metadata_dict['Processing_history'] = '-Zero-Order Holds Correction:|' \
                                          ' Correction type = Substitute with Nan' \
                                          ' Corrections applied:|' \
                                          ' All channels corrected where zero-order holds concur with Pressure Holds:|'
    metadata_dict['ZEROORDER_Time'] = datetime.now()

    #TODO: Improve metadata_dict entry here? Ask Lu

    columns = ctd.columns.tolist()

    #ctd = ctd.rename(columns={'Oxygen': 'Oxygen:Dissolved:Saturation', 'Salinity': 'Salinity:CTD', 'TIME': 'TIME:UTC'})
    column_names = dict.fromkeys(ctd.columns, '')  # set empty column names
    ctd = ctd.rename(columns=column_names)  # remove column names


    column_names_header = dict.fromkeys(header.columns, '')  # set empty column names
    header = header.rename(columns=column_names_header)
    #print(header)
    #print(ctd)

    ctd_header = header.append(ctd)
    ctd_header.to_csv(dest_dir + output_name, index=False, header=False)






def CALIB(var, var_downcast, var_upcast, correction_value, metadata_dict, zoh):
    """
     Correct pressure and depth data
     Inputs:
         - cast, downcast, upcast and metadata dictionaries
     Outputs:
         - cast, downcast, upcast and metadata dictionaries after pressure correction
     """
    n = len(var.keys())
    var1 = deepcopy(var)
    var2 = deepcopy(var_downcast)
    var3 = deepcopy(var_upcast)
    for i in range(1, n + 1, 1):
        var1['cast' + str(i)].Pressure = var1['cast' + str(i)].Pressure + correction_value
        var1['cast' + str(i)].Depth = var1['cast' + str(i)].Depth + correction_value
        var2['cast' + str(i)].Pressure = var2['cast' + str(i)].Pressure + correction_value
        var2['cast' + str(i)].Depth = var2['cast' + str(i)].Depth + correction_value
        var3['cast' + str(i)].Pressure = var3['cast' + str(i)].Pressure + correction_value
        var3['cast' + str(i)].Depth = var3['cast' + str(i)].Depth + correction_value

    # check if a correction was done - need to see if this is the first addition of "processing history'
    if zoh == 'no':
        metadata_dict['Processing_history'] = '-CALIB parameters:|' \
                                              ' Calibration type = Correct|' \
                                              ' Calibrations applied:|' \
                                              ' Pressure (decibar) = {}'.format(str(correction_value)) + '|' \
                                                                                                         ' Depth (meters) = {}'.format(
            str(correction_value)) + '|'

    elif zoh == 'yes':
        metadata_dict['Processing_history'] += '-CALIB parameters:|' \
                                               ' Calibration type = Correct|' \
                                               ' Calibrations applied:|' \
                                               ' Pressure (decibar) = {}'.format(str(correction_value)) + '|' \
                                                                                                          ' Depth (meters) = {}'.format(
            str(correction_value)) + '|'

    metadata_dict['CALIB_Time'] = datetime.now()


    return var1, var2, var3

                                                                # 0 if no neg pressures
#cast_pc, cast_d_pc, cast_u_pc = CALIB(cast, cast_d, cast_u, correction_value = 0, metadata_dict = metadata, zoh='yes')
# ------------------------------  Step 8: Data Despiking  --------------------------------------------------------------------------
# plot profiles to look for spikes






# ------------------------------  Step 9: Clip   -----------------------------------------------------------------------------------

def CLIP_DOWNCAST(var, metadata_dict, limit_drop):
    """
     CLIP the unstable measurement from sea surface and bottom
     Inputs:
         - Downcast, metadata dictionary, limit_drop,
     Outputs:
         - Downcast after removing records near surface and bottom
     """
    var_clip = deepcopy(var)
    for i in range(1, len(var_clip) + 1, 1):
        pressure = var_clip['cast' + str(i)].Pressure
        print(pressure.shape)
        diff = var_clip['cast' + str(i)].Pressure.diff()
        index_start = pressure.index[0]
        diff_drop = diff.loc[(diff > limit_drop)]
        for j in range(0, len(diff.loc[(diff > limit_drop)]), 1):
            index_1 = diff_drop.index[j]
            if (diff_drop.index[j + 1] == index_1 + 1) and (diff_drop.index[j + 2] == index_1 + 2) and (
                    diff_drop.index[j + 3] == index_1 + 3) and (diff_drop.index[j + 4] == index_1 + 4) and (
                    diff_drop.index[j + 5] == index_1 + 5) and (diff_drop.index[j + 6] == index_1 + 6) and (
                    diff_drop.index[j + 7] == index_1 + 7) and (diff_drop.index[j + 8] == index_1 + 8):
                index_end_1 = index_1 - 1
                break
        cut_start = index_end_1 - index_start
        for j in range(-1, -len(diff.loc[(diff > limit_drop)]), -1):
            index_2 = diff_drop.index[j]
            if (diff_drop.index[j - 1] == index_2 - 1) and (diff_drop.index[j - 2] == index_2 - 2) and (
                    diff_drop.index[j - 3] == index_2 - 3) and (diff_drop.index[j - 4] == index_2 - 4) and (
                    diff_drop.index[j - 5] == index_2 - 5):
                index_end_2 = index_2 + 1
                break
        cut_end = index_end_2 - index_start
        # list_start.append(cut_start)
        # list_end.append(cut_end)
        var_clip['cast' + str(i)] = var_clip['cast' + str(i)][cut_start:cut_end]
        metadata_dict['Processing_history'] += '-CLIP_downcast{}'.format(
            str(i)) + ': First Record = {}'.format(str(cut_start)) + ', Last Record = {}'.format(
            str(cut_end)) + '|'
        metadata_dict['CLIP_D_Time' + str(i)] = datetime.now()

    return var_clip


def CLIP_UPCAST(var, metadata_dict, limit_rise):
    """
     CLIP the unstable measurement from sea surface and bottom
     Inputs:
         - Upcast, metadata dictionary, limit_rise,
     Outputs:
         - Upcast after removing records near surface and bottom
     """
    var_clip = deepcopy(var)
    for i in range(1, len(var_clip) + 1, 1):
        pressure = var_clip['cast' + str(i)].Pressure
        diff = var_clip['cast' + str(i)].Pressure.diff()
        index_start = pressure.index[0]
        # index_end = pressure.index[-1]
        diff_rise = diff.loc[(diff < limit_rise)]
        for j in range(0, len(diff.loc[(diff < limit_rise)]), 1):
            index_1 = diff_rise.index[j]
            if (diff_rise.index[j + 1] == index_1 + 1) and (diff_rise.index[j + 2] == index_1 + 2) and (
                    diff_rise.index[j + 3] == index_1 + 3) and (diff_rise.index[j + 4] == index_1 + 4) and (
                    diff_rise.index[j + 5] == index_1 + 5) and (diff_rise.index[j + 6] == index_1 + 6) and (
                    diff_rise.index[j + 7] == index_1 + 7) and (diff_rise.index[j + 8] == index_1 + 8):
                index_end_1 = index_1 - 1
                break
        cut_start = index_end_1 - index_start

        for j in range(-1, -len(diff.loc[(diff < limit_rise)]), -1):
            index_2 = diff_rise.index[j]
            if (diff_rise.index[j - 1] == index_2 - 1) and (diff_rise.index[j - 2] == index_2 - 2) and (
                    diff_rise.index[j - 3] == index_2 - 3) and (diff_rise.index[j - 4] == index_2 - 4) and (
                    diff_rise.index[j - 5] == index_2 - 5):
                index_end_2 = index_2 + 1
                break
        cut_end = index_end_2 - index_start
        var_clip['cast' + str(i)] = var_clip['cast' + str(i)][cut_start:cut_end]
        metadata_dict['Processing_history'] += '-CLIP_upcast{}'.format(str(i)) + ': First Record = {}'.format(
            str(cut_start)) + ', Last Record = {}'.format(str(cut_end)) + '|'
        metadata_dict['CLIP_U_Time' + str(i)] = datetime.now()
    return var_clip


# run both functions
# cast_d_clip = CLIP_DOWNCAST(cast_d_pc, metadata, limit_drop = 0.02)
# cast_u_clip = CLIP_UPCAST(cast_u_pc, metadata, limit_rise = -0.02)


# Plot to check the profiles after clip by cast
def plot_clip(cast, cast_d_clip, cast_d_pc):
    """ plot the clipped casts to make sure they are OK"""

    fig, ax = plt.subplots()

    ax.plot(cast_d_pc['cast1'].TIME, cast_d_pc['cast1'].Pressure, color='blue', label='cast1 Pre-clip')

    ax.plot(cast_d_clip['cast1'].TIME, cast_d_clip['cast1'].Pressure, '--', color='red', label='cast1 Post_clip')
    # ax.plot(cast_u_clip['cast5'].TIME, cast_u_clip['cast5'].Pressure, color='blue', label='cast1')
    # ax.plot(cast_u_clip['cast1'].TIME, cast_u_clip['cast1'].Pressure, color='blue', label='cast1')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Time')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('After Clip')
    ax.legend()

    # plot all cast together

    number_of_colors = len(cast)
    color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

    fig, ax = plt.subplots()
    for i in range(0, len(cast), 1):
        ax.plot(cast_d_clip['cast' + str(i + 1)].TIME, cast_d_clip['cast' + str(i + 1)].Pressure, color=color[i],
                label='cast' + str(i + 1))
        # ax.plot(cast_u['cast' + str(i+1)].Salinity, cast_u['cast' + str(i+1)].Pressure, '--', color=color[i], label= 'cast' + str(i+1))
    # ax.plot(cast_d['cast1'].Salinity, cast_d['cast1'].Pressure, color='blue', label='cast1')
    # ax.plot(cast_u['cast1'].Salinity, cast_u['cast1'].Pressure, '--', color='blue', label='cast1')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Time')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('After Clip')
    ax.legend()
    plt.show()


# plot_clip(cast, cast_d_clip, cast_d_pc)
# ------------------------------  Step 10: Filter  -----------------------------------------------------------------------------------
# apply a moving average FIR filter (a simple low pass )

# def filter(x, n):# n -  filter size, 9 suggest by RBR manual, choose the smallest one which can do the job
#    b = (np.ones(n))/n #numerator co-effs of filter transfer function
#    #b = repeat(1.0/n, n)
#    a = np.ones(1)  #denominator co-effs of filter transfer function
#    #y = signal.convolve(x,b) #filter output using convolution
#    #y = signal.lfilter(b, a, x) #filter output using lfilter function
#    y = signal.filtfilt(b, a, x)  # Apply a digital filter forward and backward to a signal.
#    return y

def FILTER(var_downcast, var_upcast, window_width, sample_rate, time_constant, filter_type, metadata_dict):
    """
     Filter the profile data using a low pass filter: moving average
     Inputs:
         - downcast and upcast data dictionaries
     Outputs:
         - two dictionaries containing downcast and upcast profiles after applying filter
     """
    # filter type: 0 - FIR, 1 - moving average

    cast_number = len(var_downcast.keys())
    if filter_type == 0:
        Wn = (1.0 / time_constant) / (sample_rate * 2)
        b, a = signal.butter(2, Wn, "low")
        filter_name = "FIR"
    elif filter_type == 1:
        b = (np.ones(window_width)) / window_width  # numerator co-effs of filter transfer function
        a = np.ones(1)  # denominator co-effs of filter transfer function
        filter_name = "Moving average filter"

    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)
    for i in range(1, cast_number + 1, 1):
        var1['cast' + str(i)].Temperature = signal.filtfilt(b, a, var1['cast' + str(i)].Temperature)
        var1['cast' + str(i)].Conductivity = signal.filtfilt(b, a, var1['cast' + str(i)].Conductivity)
        var1['cast' + str(i)].Pressure = signal.filtfilt(b, a, var1['cast' + str(i)].Pressure)  # !!!!!  might need!
        # var1['cast' + str(i)].Fluorescence = signal.filtfilt(b, a, var1['cast' + str(i)].Fluorescence)
        var2['cast' + str(i)].Temperature = signal.filtfilt(b, a, var2['cast' + str(i)].Temperature)
        var2['cast' + str(i)].Conductivity = signal.filtfilt(b, a, var2['cast' + str(i)].Conductivity)
        var2['cast' + str(i)].Pressure = signal.filtfilt(b, a, var2['cast' + str(i)].Pressure)
        # var2['cast' + str(i)].Fluorescence = signal.filtfilt(b, a, var2['cast' + str(i)].Fluorescence)      #!!!!!  might need!
    metadata_dict['Processing_history'] += '-FILTER parameters:|' \
                                           ' ' + filter_name + ' was used.|' \
                                                               ' Filter width = {}'.format(str(window_width)) + '.|' \
                                                                                                                ' The following channel(s) were filtered.|' \
                                                                                                                ' Pressure|' \
                                                                                                                ' Temperature|' \
                                                                                                                ' Conductivity|'
    metadata_dict['FILTER_Time'] = datetime.now()

    return var1, var2


# cast_d_filtered, cast_u_filtered = FILTER(cast_d_clip, cast_u_clip, window_width = 6, sample_rate = 8, time_constant = 1/8, filter_type = 1, metadata_dict = metadata) # n = 5 should be good.


# plot to check values before and after filtering
def plot_filter(cast_d_filtered, cast_d_clip):
    """ check the filter plots"""

    vars = list(dict.fromkeys(cast_d_filtered['cast1']))

    fig, ax = plt.subplots()
    ax.plot(cast_d_clip['cast1'].Conductivity, cast_d_clip['cast1'].Pressure, color='blue', label='Pre-filtering')
    ax.plot(cast_d_filtered['cast1'].Conductivity, cast_d_filtered['cast1'].Pressure, '--', color='red',
            label='Post-filtering')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Conductivity')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('After Filter')
    ax.legend()
    plt.savefig(dest_dir + 'After_Filter_C')

    fig, ax = plt.subplots()
    ax.plot(cast_d_clip['cast1'].Salinity, cast_d_clip['cast1'].Pressure, color='blue', label='Pre-filtering')
    ax.plot(cast_d_filtered['cast1'].Salinity, cast_d_filtered['cast1'].Pressure, '--', color='red',
            label='Post-filtering')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('After Filter')
    ax.legend()
    plt.savefig(dest_dir + 'After_Filter_S')

    for var in vars:
        if var == 'Oxygen':
            fig, ax = plt.subplots()
            ax.plot(cast_d_clip['cast1'].Oxygen, cast_d_clip['cast1'].Pressure, color='blue', label='Pre-filtering')
            ax.plot(cast_d_filtered['cast1'].Oxygen, cast_d_filtered['cast1'].Pressure, '--', color='red',
                    label='Post-filtering')
            ax.invert_yaxis()
            ax.xaxis.set_label_position('top')
            ax.xaxis.set_ticks_position('top')
            ax.set_xlabel('Oxygen')
            ax.set_ylabel('Pressure (decibar)')
            ax.set_title('After Filter')
            ax.legend()
            plt.savefig(dest_dir + 'After_Filter_O')

    for var in vars:
        if var == 'Fluorescence':
            fig, ax = plt.subplots()
            ax.plot(cast_d_clip['cast1'].Fluorescence, cast_d_clip['cast1'].Pressure, color='blue', label='Pre-filtering')
            ax.plot(cast_d_filtered['cast1'].Fluorescence, cast_d_filtered['cast1'].Pressure, '--', color='red',
                    label='Post-filtering')
            ax.invert_yaxis()
            ax.xaxis.set_label_position('top')
            ax.xaxis.set_ticks_position('top')
            ax.set_xlabel('Fluorescence')
            ax.set_ylabel('Pressure (decibar)')
            ax.set_title('After Filter')
            ax.legend()
            plt.savefig(dest_dir + 'After_Filter_F')


    fig, ax = plt.subplots()
    ax.plot(cast_d_clip['cast1'].Temperature, cast_d_clip['cast1'].Pressure, color='blue', label='Pre-filtering')
    ax.plot(cast_d_filtered['cast1'].Temperature, cast_d_filtered['cast1'].Pressure, '--', color='red',
            label='Post-filtering')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Temperature')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('After Filter')
    ax.legend()
    plt.savefig(dest_dir + 'After_Filter_T')
    plt.show()



# plot_filter(cast_d_filtered)


# ------------------------------  Step 11: Shift conductivity and recalculate salinity ----------------------------------------------------
# input variable: cast_d_filtered, cast_u_filtered, try 2-3s (12-18 scans)

def SHIFT_CONDUCTIVITY(var_downcast, var_upcast, shifted_scan_number,
                       metadata_dict):  # n: number of scans shifted. +: delay; -: advance
    """
     Delay the conductivity signal, and recalculate salinity
     Inputs:
         - downcast and upcast data dictionaries, metadata dictionary
     Outputs:
         - two dictionaries containing downcast and upcast profiles
     """
    cast_number = len(var_downcast.keys())
    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)
    for i in range(1, cast_number + 1, 1):
        index_1 = var1['cast' + str(i)].Conductivity.index[0]
        v1 = var1['cast' + str(i)].Conductivity[index_1]
        index_2 = var2['cast' + str(i)].Conductivity.index[0]
        v2 = var2['cast' + str(i)].Conductivity[index_2]
        # shift C for n scans
        var1['cast' + str(i)].Conductivity = var1['cast' + str(i)].Conductivity.shift(periods=shifted_scan_number,
                                                                                      fill_value=v1)
        # calculates SP from C using the PSS-78 algorithm (2 < SP < 42)
        var1['cast' + str(i)].Salinity = gsw.SP_from_C(var1['cast' + str(i)].Conductivity,
                                                       var1['cast' + str(i)].Temperature,
                                                       var1['cast' + str(i)].Pressure)
        var2['cast' + str(i)].Conductivity = var2['cast' + str(i)].Conductivity.shift(periods=shifted_scan_number,
                                                                                      fill_value=v2)
        var2['cast' + str(i)].Salinity = gsw.SP_from_C(var2['cast' + str(i)].Conductivity,
                                                       var2['cast' + str(i)].Temperature,
                                                       var2['cast' + str(i)].Pressure)
    metadata_dict['Processing_history'] += '-SHIFT parameters:|' \
                                           ' Shift Channel: Conductivity|' \
                                           ' # of Records to Delay (-ve for Advance):|' \
                                           ' Shift = {}'.format(str(shifted_scan_number)) + '|' \
                                                                                            ' Salinity was recalculated after shift|'
    metadata_dict['SHIFT_Conductivity_Time'] = datetime.now()

    return var1, var2


# cast_d_shift_c, cast_u_shift_c = SHIFT_CONDUCTIVITY(cast_d_filtered, cast_u_filtered, shifted_scan_number = 2, metadata_dict = metadata)  # delay conductivity data by 2 scans

# Plot Salinity and T-S to check the index after shift
# Salinity

def plot_shift_c(cast_d_shift_c, cast_d_filtered):
    fig, ax = plt.subplots()
    ax.plot(cast_d_filtered['cast1'].Salinity, cast_d_filtered['cast1'].Pressure, color='blue', label='Pre-shift')
    #ax.plot(cast_u_filtered['cast1'].Salinity, cast_u_filtered['cast1'].Pressure, '--', color='blue', label='Pre-shift')
    ax.plot(cast_d_shift_c['cast1'].Salinity, cast_d_shift_c['cast1'].Pressure, color='red', label='Post-shift')
    #ax.plot(cast_u_shift_c['cast1'].Salinity, cast_u_shift_c['cast1'].Pressure, '--', color='red', label='Post-shift')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('After Shift Conductivity')
    plt.savefig(dest_dir + 'After_Shift_Conductivity_S.png')
    ax.legend()

    # TS Plot
    fig, ax = plt.subplots()
    ax.plot(cast_d_filtered['cast1'].Salinity, cast_d_filtered['cast1'].Temperature, color='blue', label='Pre-shift')
    #ax.plot(cast_u_filtered['cast1'].Salinity, cast_u_filtered['cast1'].Temperature, '--', color='blue',
            #label='Pre-shift')
    ax.plot(cast_d_shift_c['cast1'].Salinity, cast_d_shift_c['cast1'].Temperature, color='red', label='Post-shift')
    #ax.plot(cast_u_shift_c['cast1'].Salinity, cast_u_shift_c['cast1'].Temperature, '--', color='red',
            #label='Post-shift')
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Temperature (C)')
    ax.set_title('After Shift Conductivity T-S Plot')
    ax.legend()
    plt.savefig(dest_dir + 'After_Shift_Conductivity_T-S.png')
    plt.show()

#plot before shift Oxygen T-O

# plot_shift_c(cast_u_shift_c)

# ------------------------------    Step 12: Shift Oxygen   ----------------------------------------------------
# input variable: cast_d_filtered, cast_u_filtered, try 2-3s (12-18 scans)

def SHIFT_OXYGEN(var_downcast, var_upcast, shifted_scan_number,
                 metadata_dict):  # n: number of scans shifted. +: delay; -: advance
    """
     Advance oxygen data by 2-3s
     Inputs:
         - downcast and upcast data dictionaries, metadata dictionary
     Outputs:
         - two dictionaries containing downcast and upcast profiles
     """
    cast_number = len(var_downcast.keys())
    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)
    for i in range(1, cast_number + 1, 1):
        index_1 = var1['cast' + str(i)].Oxygen.index[-1]
        v1 = var1['cast' + str(i)].Oxygen[index_1]
        index_2 = var2['cast' + str(i)].Oxygen.index[-1]
        v2 = var2['cast' + str(i)].Oxygen[index_2]
        # shift C for n scans
        var1['cast' + str(i)].Oxygen = var1['cast' + str(i)].Oxygen.shift(periods=shifted_scan_number, fill_value=v1)
        var2['cast' + str(i)].Oxygen = var2['cast' + str(i)].Oxygen.shift(periods=shifted_scan_number, fill_value=v2)
    metadata_dict['Processing_history'] += '-SHIFT parameters:|' \
                                           ' Shift Channel: Oxygen:Dissolved:Saturation|' \
                                           ' # of Records to Delay (-ve for Advance):|' \
                                           ' Shift = {}'.format(str(shifted_scan_number)) + '|'
    metadata_dict['SHIFT_Oxygen_Time'] = datetime.now()

    return var1, var2


# try:
# oxy = cast_d['Oxygen'].values
# cast_d_shift_o, cast_u_shift_o = SHIFT_OXYGEN(cast_d_shift_c, cast_u_shift_c, shifted_scan_number = -11, metadata_dict = metadata)  # advance oxygen data by 11 scans
# except KeyError:
# cast_d_shift_o, cast_u_shift_o = cast_d_shift_c, cast_u_shift_c

# plot T-O2 before and after shift to check the results
# T-O2 Plot to check whether the shift bring the downcas and upcast together

# for var in vars:
# if var == 'Oxygen':
def plot_shift_o(cast_d_shift_o, cast_d_shift_c):
    """Check Oxy plots after shift """

    fig, ax = plt.subplots()
    ax.plot(cast_d_shift_c['cast1'].Temperature, cast_d_shift_c['cast1'].Oxygen, color='blue', label='Pre-shift')
    #ax.plot(cast_u_shift_c['cast2'].Temperature, cast_u_shift_c['cast2'].Oxygen, '--', color='blue', label='Pre-shift')
    ax.plot(cast_d_shift_o['cast1'].Temperature, cast_d_shift_o['cast1'].Oxygen, color='red', label='Post-shift')
    #ax.plot(cast_u_shift_o['cast2'].Temperature, cast_u_shift_o['cast2'].Oxygen, '--', color='red', label='Post-shift')
    ax.set_ylabel('Oxygen Saturation (%)')
    ax.set_xlabel('Temperature (C)')
    ax.set_title('After Shift Oxygen T-O Plot')
    ax.legend()
    plt.savefig(dest_dir + 'After_Shift_Oxygen_T-O.png')
    plt.show()


# plot_shift_o(cast_d_shift_o)

# ------------------------------    Step 13: Delete (swells/slow drop)  ----------------------------------------------------
# correct for the wake effect, remove the pressure reversal

def DELETE_PRESSURE_REVERSAL(var_downcast, var_upcast, metadata_dict):
    """
     Detect and delete pressure reversal
     Inputs:
         - downcast and upcast data dictionaries, metadata dictionary
     Outputs:
         - two dictionaries containing downcast and upcast profiles
     """

    cast_number = len(var_downcast.keys())
    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)
    for i in range(1, cast_number + 1, 1):
        press = var1['cast' + str(i)].Pressure.values
        ref = press[0]
        inversions = np.diff(np.r_[press, press[-1]]) < 0
        mask = np.zeros_like(inversions)
        for k, p in enumerate(inversions):
            if p:
                ref = press[k]
                cut = press[k + 1:] < ref
                mask[k + 1:][cut] = True
        var1['cast' + str(i)][mask] = np.NaN

    for i in range(1, cast_number + 1, 1):
        press = var2['cast' + str(i)].Pressure.values
        ref = press[0]
        inversions = np.diff(np.r_[press, press[-1]]) > 0
        mask = np.zeros_like(inversions)
        for k, p in enumerate(inversions):
            if p:
                ref = press[k]
                cut = press[k + 1:] > ref
                mask[k + 1:][cut] = True
        var2['cast' + str(i)][mask] = np.NaN
    metadata_dict['Processing_history'] += '-DELETE_PRESSURE_REVERSAL parameters:|' \
                                           ' Remove pressure reversals|'
    metadata_dict['DELETE_PRESSURE_REVERSAL_Time'] = datetime.now()

    return var1, var2


# cast_d_wakeeffect, cast_u_wakeeffect = DELETE_PRESSURE_REVERSAL(cast_d_shift_o, cast_u_shift_o, metadata_dict = metadata)

# Plot Salinity and T-S to check the index after shift
# Salinity

def plot_processed(cast_d_wakeeffect, cast_d_shift_o, cast_d, cast):
    """Plot the processed casts"""
    vars = list(dict.fromkeys(cast['cast1']))
    fig, ax = plt.subplots()

    ax.plot(cast_d_shift_o['cast1'].Salinity, cast_d_shift_o['cast1'].Pressure, color='blue', label='Pre-Delete')
    # ax.plot(cast_u_filtered['cast1'].Salinity, cast_u_filtered['cast1'].Pressure, '--', color='blue', label='Pre-shift')
    ax.plot(cast_d_wakeeffect['cast1'].Salinity, cast_d_wakeeffect['cast1'].Pressure, color='red', label='Post-Delete')
    # ax.plot(cast_u_wakeeffect['cast1'].Salinity, cast_u_wakeeffect['cast1'].Pressure, '--', color='red', label='Post-shift')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel(' ')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('After Delete Pressure Reversal')
    ax.legend()
    plt.savefig(dest_dir + 'After_Delete_S.png')

    fig, ax = plt.subplots()
    ax.plot(cast_d_shift_o['cast1'].Conductivity, cast_d_shift_o['cast1'].Pressure, color='blue', label='Pre-Delete')
    # ax.plot(cast_u_filtered['cast1'].Salinity, cast_u_filtered['cast1'].Pressure, '--', color='blue', label='Pre-shift')
    ax.plot(cast_d_wakeeffect['cast1'].Conductivity, cast_d_wakeeffect['cast1'].Pressure, color='red', label='Post-Delete')
    # ax.plot(cast_u_wakeeffect['cast1'].Salinity, cast_u_wakeeffect['cast1'].Pressure, '--', color='red', label='Post-shift')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Conductivity')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('After Delete Pressure Reversal')
    ax.legend()
    plt.savefig(dest_dir + 'After_Delete_C.png')

    # TS Plot
    fig, ax = plt.subplots()
    ax.plot(cast_d_shift_o['cast1'].Salinity, cast_d_shift_o['cast1'].Temperature, color='blue', label='Pre-Delete')
    # ax.plot(cast_u_filtered['cast1'].Salinity, cast_u_filtered['cast1'].Temperature, '--', color='blue', label='Pre-shift')
    ax.plot(cast_d_wakeeffect['cast1'].Salinity, cast_d_wakeeffect['cast1'].Temperature, color='red',
            label='Post-Delete')
    # ax.plot(cast_u_wakeeffect['cast1'].Salinity, cast_u_wakeeffect['cast1'].Temperature, '--', color='red', label='Post-shift')
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Temperature (C)')
    ax.set_title('T-S Plot (after delete pressure reversal)')
    ax.legend()
    plt.savefig(dest_dir + 'After_Delete_T-S.png')

    # -------------------------  Plot processed profiles    ------------------------------------------------------------
    number_of_colors = len(cast)
    color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

    # plot salinity of all cast together
    fig, ax = plt.subplots()
    # for i in range (0, len(cast_d), 1):
    # ax.plot(cast_d['cast' + str(i+1)].Salinity, cast_d['cast' + str(i+1)].Pressure, color=color[i], label= 'cast' + str(i+1))
    # ax.plot(cast_d_wakeeffect['cast' + str(i+1)].Salinity, cast_d_wakeeffect['cast' + str(i+1)].Pressure, '--', color=color[i], label='cast' + str(i+1))
    ax.plot(cast_d['cast1'].Salinity, cast_d['cast1'].Pressure, color='blue', label='Pre_Processing')
    # ax.plot(cast_u['cast1'].Salinity, cast_u['cast1'].Pressure, '--', color='blue', label='Before_processing')
    ax.plot(cast_d_wakeeffect['cast1'].Salinity, cast_d_wakeeffect['cast1'].Pressure, color='red',
            label='After_processing')
    # ax.plot(cast_u_wakeeffect['cast1'].Salinity, cast_u_wakeeffect['cast1'].Pressure, '--', color='red', label='cast1')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('Pre and Post Processing')
    ax.legend()
    plt.savefig(dest_dir + 'pre_post_S.png')

    # plot temperature of all cast together
    fig, ax = plt.subplots()
    # for i in range (0, len(cast_d), 1):
    #    #ax.plot(cast_d['cast' + str(i+1)].Temperature, cast_d['cast' + str(i+1)].Pressure, color=color[i], label= 'cast' + str(i+1))
    #    ax.plot(cast_d_wakeeffect['cast' + str(i+1)].Temperature, cast_d_wakeeffect['cast' + str(i+1)].Pressure, '--', color=color[i], label='cast' + str(i+1))
    ax.plot(cast_d['cast1'].Temperature, cast_d['cast1'].Pressure, color='blue', label='Pre_Processing')
    # ax.plot(cast_u['cast1'].Temperature, cast_u['cast1'].Pressure, '--', color='blue', label='cast1')
    ax.plot(cast_d_wakeeffect['cast1'].Temperature, cast_d_wakeeffect['cast1'].Pressure, color='red', label='Post_Processing')
    # ax.plot(cast_u_wakeeffect['cast1'].Temperature, cast_u_wakeeffect['cast1'].Pressure, '--', color='red', label='cast1')
    ax.invert_yaxis()
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Temperature(C)')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('Pre and Post Processing')
    ax.legend()
    plt.savefig(dest_dir + 'pre_post_T.png')

    # plot Conductivity of all cast together
    fig, ax = plt.subplots()
    #for i in range(0, len(cast_d), 1):
        # ax.plot(cast_d['cast' + str(i+1)].Conductivity, cast_d['cast' + str(i+1)].Pressure, color=color[i], label= 'cast' + str(i+1))
        #ax.plot(cast_d_wakeeffect['cast' + str(i + 1)].Conductivity, cast_d_wakeeffect['cast' + str(i + 1)].Pressure,
                #color=color[i], label='cast' + str(i + 1))
    ax.plot(cast_d['cast1'].Conductivity, cast_d['cast1'].Pressure, color='blue', label='Pre_Processing')
    # ax.plot(cast_u['cast2'].Conductivity[1:], cast_u['cast2'].Pressure[1:], '--', color='red', label='cast2')
    ax.plot(cast_d_wakeeffect['cast1'].Conductivity, cast_d_wakeeffect['cast1'].Pressure, color='red', label='Post_Processing')
    # ax.plot(cast_u_wakeeffect['cast1'].Conductivity, cast_u_wakeeffect['cast1'].Pressure, color='red', label='cast1')
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_xlabel('Conductivity (S/cm)')
    ax.set_ylabel('Pressure (decibar)')
    ax.set_title('Pre and Post Processing')
    ax.legend()
    plt.savefig(dest_dir + 'pre_post_C.png')

    # plot Oxygen of all cast together
    fig, ax = plt.subplots()
    for var in vars:
        if var == 'Oxygen':
            #for i in range(0, len(cast_d), 1):
                # ax.plot(cast_d['cast' + str(i+1)].Oxygen, cast_d['cast' + str(i+1)].Pressure, color=color[i], label= 'cast' + str(i+1))
                #ax.plot(cast_d_wakeeffect['cast' + str(i + 1)].Oxygen, cast_d_wakeeffect['cast' + str(i + 1)].Pressure,
                        #color=color[i], label='cast' + str(i + 1))
            ax.plot(cast_d['cast1'].Oxygen, cast_d['cast1'].Pressure, color='blue', label='Pre_Processing')
                # ax.plot(cast_u['cast1'].Oxygen, cast_u['cast1'].Pressure, '--', color='blue', label='cast1')
            ax.plot(cast_d_wakeeffect['cast1'].Oxygen, cast_d_wakeeffect['cast1'].Pressure, color='red', label='Post_Processing')
                # ax.plot(cast_d_wakeeffect['cast1'].Oxygen, cast_d_wakeeffect['cast1'].Pressure, color='red', label='cast1')
            ax.invert_yaxis()
            ax.xaxis.set_label_position('top')
            ax.xaxis.set_ticks_position('top')
            ax.set_xlabel('Oxygen Saturation (%)')  # Check unit here
            ax.set_ylabel('Pressure (decibar)')
            ax.set_title('Pre and Post Processing')
            ax.legend()
            plt.savefig(dest_dir + 'pre_post_O.png')

    # plot Fluorescence of all cast together
    fig, ax = plt.subplots()
    for var in vars:
        if var == 'Fluorescence':
            #for i in range(0, len(cast_d), 1):
                # ax.plot(cast_d['cast' + str(i+1)].Fluorescence, cast_d['cast' + str(i+1)].Pressure, color=color[i], label= 'cast' + str(i+1))
               # ax.plot(cast_d_wakeeffect['cast' + str(i + 1)].Fluorescence,
                        #cast_d_wakeeffect['cast' + str(i + 1)].Pressure, '--', color=color[i],
                        #label='cast' + str(i + 1))
            ax.plot(cast_d['cast1'].Fluorescence, cast_d['cast1'].Pressure, color='blue', label='Pre Processing')
                # ax.plot(cast_u['cast1'].Fluorescence, cast_u['cast1'].Pressure, '--', color='blue', label='cast1')
            ax.plot(cast_d_wakeeffect['cast1'].Fluorescence, cast_d_wakeeffect['cast1'].Pressure, color='red', label='Post Processing')
                # ax.plot(cast_u_wakeeffect['cast1'].Fluorescence, cast_u_wakeeffect['cast1'].Pressure, color='red', label='cast1')
            ax.invert_yaxis()
            ax.xaxis.set_label_position('top')
            ax.xaxis.set_ticks_position('top')
            ax.set_xlabel('Fluorescence (ug/L)')  # Check unit here
            ax.set_ylabel('Pressure (decibar)')
            ax.set_title('Pre and Post processing')
            ax.legend()
            plt.savefig(dest_dir + 'pre_post_F.png')

    fig, ax = plt.subplots()
    # for i in range (0, len(cast_d), 1):
    #    #ax.plot(cast_d['cast' + str(i+1)].Salinity, cast_d['cast' + str(i+1)].Temperature, color=color[i], label= 'cast' + str(i+1))
    #    ax.plot(cast_d_wakeeffect['cast' + str(i+1)].Salinity, cast_d_wakeeffect['cast' + str(i+1)].Temperature, color=color[i], label='cast' + str(i+1))
    ax.plot(cast_d['cast1'].Salinity, cast_d['cast1'].Temperature, color='blue', label='Pre-Processing')
    # ax.plot(cast_d_shift_o['cast1'].Salinity, cast_d_shift_o['cast1'].Temperature, color='blue', label='Pre-shift')
    # ax.plot(cast_u_filtered['cast1'].Salinity, cast_u_filtered['cast1'].Temperature, '--', color='blue', label='Pre-shift')
    ax.plot(cast_d_wakeeffect['cast1'].Salinity, cast_d_wakeeffect['cast1'].Temperature, color='red',
            label='Post-Processing')
    # ax.plot(cast_u_wakeeffect['cast1'].Salinity, cast_u_wakeeffect['cast1'].Temperature, '--', color='red', label='Post-shift')
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Temperature (C)')
    ax.set_title('Pre and Post Processing T-S Plot')
    ax.legend()
    plt.savefig(dest_dir + 'pre_post_T-S.png')
    plt.show()


# plot_processed(cast_d_wakeeffect)

# ------------------------------    Step 14: bin averages  ----------------------------------------------------
# input variables: cast_d_wakeeffect, cast_u_wakeeffect
def BINAVE(var_downcast, var_upcast, interval, metadata_dict):
    """
     Bin average the profiles
     Note: Bin width and spacing are both universally chosen to be 1m in coastal waters
     Inputs:
         - downcast and upcast data dictionaries, metadata dictionary
     Outputs:
         - two dictionaries containing downcast and upcast profiles
     """
    cast_number = len(var_downcast.keys())
    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)
    for i in range(1, cast_number + 1, 1):
        start_d = np.floor(np.nanmin(var1['cast' + str(i)].Pressure.values))
        #start_d = np.round(start_d)
        stop_d = np.ceil(np.nanmax(var1['cast' + str(i)].Pressure.values))
        #stop_d = np.round(stop_d)
        new_press_d = np.arange(start_d - 0.5, stop_d + 1.5, interval)
        binned_d = pd.cut(var1['cast' + str(i)].Pressure, bins=new_press_d)
        obs_count_d = var1['cast' + str(i)].groupby(binned_d).size()
        var1['cast' + str(i)] = var1['cast' + str(i)].groupby(binned_d).mean()
        var1['cast' + str(i)]['Observation_counts'] = obs_count_d
        # Potential for whole row Nan values at top and bottom of output files
        var1['cast' + str(i)] = var1['cast' + str(i)].dropna(axis=0, how='any')  #drop the nans - ask if this is OK?
        var1['cast' + str(i)].reset_index(drop=True, inplace=True)

        start_u = np.ceil(np.nanmax(var2['cast' + str(i)].Pressure.values))
        stop_u = np.floor(np.nanmin(var2['cast' + str(i)].Pressure.values))
        new_press_u = np.arange(start_u + 0.5, stop_u - 1.5, -interval)
        binned_u = pd.cut(var2['cast' + str(i)].Pressure, bins=new_press_u[::-1])
        obs_count_u = var2['cast' + str(i)].groupby(binned_u).size()
        var2['cast' + str(i)] = var2['cast' + str(i)].groupby(binned_u).mean()
        var2['cast' + str(i)] = var2['cast' + str(i)].sort_values('Depth', ascending=False)
        var2['cast' + str(i)]['Observation_counts'] = obs_count_u
        var2['cast' + str(i)] = var2['cast' + str(i)].dropna(axis=0, how='any')
        var2['cast' + str(i)].reset_index(drop=True, inplace=True)
    metadata_dict['Processing_history'] += '-BINAVE parameters:' \
                                           ' Bin channel = Pressure|' \
                                           ' Averaging interval = 1.00|' \
                                           ' Minimum bin value = 0.000|' \
                                           ' Average value were used|' \
                                           ' Interpolated values were NOT used for empty bins|' \
                                           ' Channel NUMBER_OF_BIN_RECORDS was added to file|'
    metadata_dict['BINAVE_Time'] = datetime.now()
    return var1, var2


# cast_d_binned, cast_u_binned = BINAVE(cast_d_wakeeffect, cast_u_wakeeffect, interval = 1, metadata_dict = metadata)


# ------------------------------    Step 15: Final edits  ----------------------------------------------------
# input variables: cast_d_binned, cast_u_binned

def FINAL_EDIT(var_cast, oxy, metadata_dict):
    """
     Final editing the profiles: edit header information, correct the unit of conductivity
     Inputs:
         - downcast and upcast data dictionaries, metadata dictionary
     Outputs:
         - two dictionaries containing downcast and upcast profiles
     """

    vars = list(dict.fromkeys(var_cast['cast1']))
    cast_number = len(var_cast.keys())
    var = deepcopy(var_cast)

    if oxy == 'yes':
        col_list = ['Pressure', 'Depth', 'Temperature', 'Salinity', 'Fluorescence', 'Oxygen', 'Conductivity',
                    'Observation_counts']
    elif oxy == 'no':
        col_list = ['Pressure', 'Depth', 'Temperature', 'Salinity', 'Conductivity', 'Observation_counts']
    for i in range(1, cast_number + 1, 1):
        var['cast' + str(i)] = var['cast' + str(i)].reset_index(drop=True)  # drop index column
        var['cast' + str(i)] = var['cast' + str(i)][col_list]  # select columns
        var['cast' + str(i)].Conductivity = var['cast' + str(i)].Conductivity * 0.1  # convert Conductivity to S/m

        var['cast' + str(i)].Pressure = var['cast' + str(i)].Pressure.apply('{:,.1f}'.format)
        var['cast' + str(i)].Depth = var['cast' + str(i)].Depth.apply('{:,.1f}'.format)
        var['cast' + str(i)].Temperature = var['cast' + str(i)].Temperature.apply('{:,.4f}'.format)
        var['cast' + str(i)].Salinity = var['cast' + str(i)].Salinity.apply('{:,.4f}'.format)
        for var_item in vars:
            if var_item == 'Fluorescence':
                var['cast' + str(i)].Fluorescence = var['cast' + str(i)].Fluorescence.apply('{:,.3f}'.format)
        for var_item in vars:
            if var_item == 'Oxygen':
                var['cast' + str(i)].Oxygen = var['cast' + str(i)].Oxygen.apply('{:,.2f}'.format)
        var['cast' + str(i)].Conductivity = var['cast' + str(i)].Conductivity.apply('{:,.6f}'.format)

        if oxy == 'yes':
            var['cast' + str(i)].columns = ['Pressure', 'Depth', 'Temperature', 'Salinity', 'Fluorescence:URU',
                                            'Oxygen:Dissolved:Saturation:RBR', 'Conductivity', 'Number_of_bin_records']
        elif oxy == 'no':
            var['cast' + str(i)].columns = ['Pressure', 'Depth', 'Temperature', 'Salinity', 'Conductivity',
                                            'Number_of_bin_records']
    metadata_dict['Processing_history'] += '-Remove Channels:|' \
                                           ' The following CHANNEL(S) were removed:|' \
                                           ' Date|' \
                                           ' TIME:UTC|' \
                                           '-CALIB parameters:|' \
                                           ' Calobration type = Correct|' \
                                           ' Calibration applied:|' \
                                           ' Conductivity (S/m) = 0.1* Conductivity (mS/cm)|'
    metadata_dict['FINALEDIT_Time'] = datetime.now()
    return var


# cast_d_final = FINAL_EDIT(cast_d_binned, metadata_dict=metadata)


# ----------------------------    IOS Header File   ---------------------------------------------------------------------

# define function to write file section
def write_file(cast_number, cast_original, cast_final, oxy, metadata_dict):
    """
     Bin average the profiles
     Inputs:
         - cast_number, cast_original = cast, cast_final = cast_d_final, metadata_dict = metadata
     Outputs:
         - two dictionaries containing downcast and upcast profiles
     """

    vars = list(dict.fromkeys(cast_original['cast1']))

    start_time = pd.to_datetime(cast_original['cast' + str(cast_number)].Date.values[0] + ' ' +
                                cast_original['cast' + str(cast_number)].TIME.values[0]).strftime(
        "%Y/%m/%d %H:%M:%S.%f")[0:-3]
    end_time = pd.to_datetime(cast_original['cast' + str(cast_number)].Date.values[-1] + ' ' +
                              cast_original['cast' + str(cast_number)].TIME.values[-1]).strftime(
        "%Y/%m/%d %H:%M:%S.%f")[0:-3]


    sample_interval = metadata_dict['Sampling_Interval']
    time_increment = '0 0 0 ' + sample_interval + ' 0  ! (day hr min sec ms)'

    number_of_records = str(cast_final['cast' + str(cast_number)].shape[0])  # number of ensumbles
    data_description = metadata_dict['Data_description']
    number_of_channels = str(cast_final['cast' + str(cast_number)].shape[1])
    nan = -99
    file_type = "ASCII"

    print("*FILE")
    print("    " + '{:20}'.format('START TIME') + ": UTC " + start_time)
    print("    " + '{:20}'.format('END TIME') + ": UTC " + end_time)
    print("    " + '{:20}'.format('TIME INCREMENT') + ": " + time_increment)
    print("    " + '{:20}'.format('NUMBER OF RECORDS') + ": " + number_of_records)
    print("    " + '{:20}'.format('DATA DESCRIPTION') + ": " + data_description)
    print("    " + '{:20}'.format('FILE TYPE') + ": " + file_type)
    print("    " + '{:20}'.format('NUMBER OF CHANNELS') + ": " + number_of_channels)
    print()
    print('{:>20}'.format('$TABLE: CHANNELS'))
    print('    ' + '! No Name                             Units          Minimum        Maximum')
    print('    ' + '!--- -------------------------------- -------------- -------------- --------------')

    print('{:>8}'.format('1') + " " + '{:33}'.format(
        list(cast_final['cast' + str(cast_number)].columns)[0]) + '{:15}'.format(
        "decibar") + '{:15}'.format(
        str(np.nanmin(cast_final['cast' + str(cast_number)].Pressure.astype(np.float)))) + '{:14}'.format(
        str(np.nanmax(cast_final['cast' + str(cast_number)].Pressure.astype(np.float)))))

    print('{:>8}'.format('2') + " " + '{:33}'.format(
        list(cast_final['cast' + str(cast_number)].columns)[1]) + '{:15}'.format(
        "meters") + '{:15}'.format(
        str(np.nanmin(cast_final['cast' + str(cast_number)].Depth.astype(np.float)))) + '{:14}'.format(
        str(np.nanmax(cast_final['cast' + str(cast_number)].Depth.astype(np.float)))))

    print('{:>8}'.format('3') + " " + '{:33}'.format(
        list(cast_final['cast' + str(cast_number)].columns)[2]) + '{:15}'.format(
        "'deg C(ITS90)'") + '{:15}'.format(
        str(np.nanmin(cast_final['cast' + str(cast_number)].Temperature.astype(np.float)))) + '{:14}'.format(
        str(np.nanmax(cast_final['cast' + str(cast_number)].Temperature.astype(np.float)))))

    print('{:>8}'.format('4') + " " + '{:33}'.format(
        list(cast_final['cast' + str(cast_number)].columns)[3]) + '{:15}'.format(
        "PSS-78") + '{:15}'.format(
        str(np.nanmin(cast_final['cast' + str(cast_number)].Salinity.astype(np.float)))) + '{:14}'.format(
        str(float('%.04f' % np.nanmax(cast_final['cast' + str(cast_number)].Salinity.astype(np.float))))))

    if oxy == 'yes':
        print('{:>8}'.format('5') + " " + '{:33}'.format(
            list(cast_final['cast' + str(cast_number)].columns)[4]) + '{:15}'.format(
            "mg/m^3") + '{:15}'.format(str(np.nanmin(
            cast_final['cast' + str(cast_number)]['Fluorescence:URU'].astype(np.float)))) + '{:14}'.format(str(float(
            '%.03f' % np.nanmax(cast_final['cast' + str(cast_number)]['Fluorescence:URU'].astype(np.float))))))

        print('{:>8}'.format('6') + " " + '{:33}'.format(
            list(cast_final['cast' + str(cast_number)].columns)[5]) + '{:15}'.format(
            "%") + '{:15}'.format(str(np.nanmin(
            cast_final['cast' + str(cast_number)]['Oxygen:Dissolved:Saturation:RBR'].astype(
                np.float)))) + '{:14}'.format(str(float('%.04f' % np.nanmax(
            cast_final['cast' + str(cast_number)]['Oxygen:Dissolved:Saturation:RBR'].astype(np.float))))))

        print('{:>8}'.format('7') + " " + '{:33}'.format(
            list(cast_final['cast' + str(cast_number)].columns)[6]) + '{:15}'.format(
            "S/m") + '{:15}'.format(
            str(np.nanmin(cast_final['cast' + str(cast_number)].Conductivity.astype(np.float)))) + '{:14}'.format(
            str(float('%.05f' % np.nanmax(cast_final['cast' + str(cast_number)].Conductivity.astype(np.float))))))

        print('{:>8}'.format('8') + " " + '{:33}'.format(
            list(cast_final['cast' + str(cast_number)].columns)[7]) + '{:15}'.format(
            "n/a") + '{:15}'.format(str(np.nanmin(
            cast_final['cast' + str(cast_number)]['Number_of_bin_records'].astype(np.float)))) + '{:14}'.format(
            str(np.nanmax(cast_final['cast' + str(cast_number)]['Number_of_bin_records'].astype(np.float)))))

    elif oxy == 'no':
        print('{:>8}'.format('5') + " " + '{:33}'.format(
            list(cast_final['cast' + str(cast_number)].columns)[4]) + '{:15}'.format(
            "S/m") + '{:15}'.format(
            str(np.nanmin(cast_final['cast' + str(cast_number)].Conductivity.astype(np.float)))) + '{:14}'.format(
            str(float('%.05f' % np.nanmax(cast_final['cast' + str(cast_number)].Conductivity.astype(np.float))))))

        print('{:>8}'.format('6') + " " + '{:33}'.format(
            list(cast_final['cast' + str(cast_number)].columns)[5]) + '{:15}'.format(
            "n/a") + '{:15}'.format(str(np.nanmin(
            cast_final['cast' + str(cast_number)]['Number_of_bin_records'].astype(np.float)))) + '{:14}'.format(
            str(np.nanmax(cast_final['cast' + str(cast_number)]['Number_of_bin_records'].astype(np.float)))))

    # Add in table of Channel summary
    print('{:>8}'.format('$END'))
    print()
    print('{:>26}'.format('$TABLE: CHANNEL DETAIL'))
    print('    ' + '! No  Pad   Start  Width  Format  Type  Decimal_Places')
    print('    ' + '!---  ----  -----  -----  ------  ----  --------------')
    # print('{:>8}'.format('1') + "  " + '{:15}'.format("' '") + '{:7}'.format(' ') + '{:7}'.format("' '") + '{:22}'.format('YYYY-MM-DDThh:mm:ssZ') + '{:6}'.format('D, T') + '{:14}'.format("' '"))
    print(
        '{:>8}'.format('1') + "  " + '{:6}'.format(str(nan)) + '{:7}'.format("' '") + '{:7}'.format(
            str(7)) + '{:8}'.format(
            'F') + '{:6}'.format('R4') + '{:3}'.format(1))
    print(
        '{:>8}'.format('2') + "  " + '{:6}'.format(str(nan)) + '{:7}'.format("' '") + '{:7}'.format(
            str(7)) + '{:8}'.format(
            'F') + '{:6}'.format('R4') + '{:3}'.format(1))
    print(
        '{:>8}'.format('3') + "  " + '{:6}'.format(str(nan)) + '{:7}'.format("' '") + '{:7}'.format(
            str(9)) + '{:8}'.format(
            'F') + '{:6}'.format('R4') + '{:3}'.format(4))
    print(
        '{:>8}'.format('4') + "  " + '{:6}'.format(str(nan)) + '{:7}'.format("' '") + '{:7}'.format(
            str(9)) + '{:8}'.format(
            'F') + '{:6}'.format('R4') + '{:3}'.format(4))
    for var_item in vars:
        if var_item == 'Fluorescence':
            print(
                '{:>8}'.format('5') + "  " + '{:6}'.format(str(nan)) + '{:7}'.format("' '") + '{:7}'.format(
                    str(8)) + '{:8}'.format(
                    'F') + '{:6}'.format('R4') + '{:3}'.format(3))
    for var_item in vars:
        if var_item == 'Oxygen':
            print(
                '{:>8}'.format('6') + "  " + '{:6}'.format(str(nan)) + '{:7}'.format("' '") + '{:7}'.format(
                    str(8)) + '{:8}'.format(
                    'F') + '{:6}'.format('R4') + '{:3}'.format(2))
    print(
        '{:>8}'.format('7') + "  " + '{:6}'.format(str(nan)) + '{:7}'.format("' '") + '{:7}'.format(
            str(10)) + '{:8}'.format(
            'F') + '{:6}'.format('R4') + '{:3}'.format(6))
    print(
        '{:>8}'.format('8') + "  " + '{:6}'.format("' '") + '{:7}'.format("' '") + '{:7}'.format(
            str(5)) + '{:8}'.format(
            'I') + '{:6}'.format('I') + '{:3}'.format(0))
    # Add in table of Channel detail summary
    print('{:>8}'.format('$END'))
    print()


# define function to write administation section
def write_admin(metadata_dict):
    mission = metadata_dict["Mission"]
    agency = metadata_dict["Agency"]
    country = metadata_dict["Country"]
    project = metadata_dict["Project"]
    scientist = metadata_dict["Scientist"]
    platform = metadata_dict["Platform"]
    print("*ADMINISTRATION")
    print("    " + '{:20}'.format('MISSION') + ": " + mission)
    print("    " + '{:20}'.format('AGENCY') + ": " + agency)
    print("    " + '{:20}'.format('COUNTRY') + ": " + country)
    print("    " + '{:20}'.format('PROJECT') + ": " + project)
    print("    " + '{:20}'.format('SCIENTIST') + ": " + scientist)
    print("    " + '{:20}'.format('PLATFORM ') + ": " + platform)
    print()


def write_location(cast_number, metadata_dict):
    """
     write location part in IOS header file
     Inputs:
         - cast_number, metadata_list
     Outputs:
         - part of txt file
     """
    station_number = metadata_dict['Location']['LOC:STATION'].tolist()
    event_number = metadata_dict['Location']['LOC:Event Number'].tolist()
    lon = metadata_dict['Location']['LOC:LONGITUDE'].tolist()
    lat = metadata_dict['Location']['LOC:LATITUDE'].tolist()
    print("*LOCATION")
    print("    " + '{:20}'.format('STATION') + ": " + str(station_number[cast_number - 1]))
    print("    " + '{:20}'.format('EVENT NUMBER') + ": " + str(event_number[cast_number - 1]))
    print("    " + '{:20}'.format('LATITUDE') + ":  " + lat[cast_number - 1][0:10] + "0 " + lat[cast_number - 1][
                                                                                            -14:-1] + ")")
    print("    " + '{:20}'.format('LONGITUDE') + ": " + lon[cast_number - 1])
    print()


# define function to write instrument info
def write_instrument(metadata_dict):
    model = metadata_dict['Instrument_Model']
    serial_number = f'{0:0}' + metadata_dict['Serial_number']
    data_description = metadata_dict['Data_description']
    instrument_type = metadata_dict['Instrument_type']
    print("*INSTRUMENT")
    print("    MODEL               : " + model)
    print("    SERIAL NUMBER       : " + serial_number)
    print("    INSTRUMENT TYPE     : " + instrument_type + "                           ! custom item")
    print("    DATA DESCRIPTION    : " + data_description + "                               ! custom item")
    print()


# define function to write raw info
# def write_history(cast_original, cast_clip, cast_filtered, cast_shift_c, cast_shift_o, cast_wakeeffect, cast_binned, cast_final, metadata_dict, cast_number):
def write_history(cast_original, cast_clip, cast_filtered, cast_shift_c, cast_shift_o, cast_wakeeffect, cast_binned,
                  cast_final, metadata_dict, cast_number):
    vars = list(dict.fromkeys(cast_original['cast1']))
    print("*HISTORY")
    print()
    print("    $TABLE: PROGRAMS")
    print("    !   Name     Vers   Date       Time     Recs In   Recs Out")
    print("    !   -------- ------ ---------- -------- --------- ---------")
    print("        Z ORDER  " + '{:7}'.format(str(1.0))
          + '{:11}'.format(metadata['ZEROORDER_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[0])
          + '{:9}'.format(metadata['ZEROORDER_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[1])
          + '{:>9}'.format(str(cast_original['cast' + str(cast_number)].shape[0]))
          + '{:>10}'.format(str(cast_original['cast' + str(cast_number)].shape[0])))
    print("        CALIB    " + '{:7}'.format(str(1.0))
          + '{:11}'.format(metadata['CALIB_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[0])
          + '{:9}'.format(metadata['CALIB_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[1])
          + '{:>9}'.format(str(cast_original['cast' + str(cast_number)].shape[0]))
          + '{:>10}'.format(str(cast_original['cast' + str(cast_number)].shape[0])))
    print("        CLIP     " + '{:7}'.format(str(1.0))
          + '{:11}'.format(
        metadata['CLIP_D_Time' + str(cast_number)].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[0])
          + '{:9}'.format(
        metadata['CLIP_D_Time' + str(cast_number)].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[1])
          + '{:>9}'.format(str(cast_original['cast' + str(cast_number)].shape[0]))
          + '{:>10}'.format(str(cast_clip['cast' + str(cast_number)].shape[0])))
    print("        FILTER   " + '{:7}'.format(str(1.0))
          + '{:11}'.format(metadata['FILTER_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[0])
          + '{:9}'.format(metadata['FILTER_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[1])
          + '{:>9}'.format(str(cast_clip['cast' + str(cast_number)].shape[0]))
          + '{:>10}'.format(str(cast_filtered['cast' + str(cast_number)].shape[0])))
    print("        SHIFT    " + '{:7}'.format(str(1.0))
          + '{:11}'.format(metadata['SHIFT_Conductivity_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[0])
          + '{:9}'.format(metadata['SHIFT_Conductivity_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[1])
          + '{:>9}'.format(str(cast_filtered['cast' + str(cast_number)].shape[0]))
          + '{:>10}'.format(str(cast_shift_c['cast' + str(cast_number)].shape[0])))
    for var_item in vars:
        if var_item == 'Oxygen':
            print("        SHIFT    " + '{:7}'.format(str(1.0))
                  + '{:11}'.format(metadata['SHIFT_Oxygen_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[0])
                  + '{:9}'.format(metadata['SHIFT_Oxygen_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[1])
                  + '{:>9}'.format(str(cast_shift_c['cast' + str(cast_number)].shape[0]))
                  + '{:>10}'.format(str(cast_shift_o['cast' + str(cast_number)].shape[0])))
    print("        DELETE   " + '{:7}'.format(str(1.0))
          + '{:11}'.format(
        metadata['DELETE_PRESSURE_REVERSAL_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[0])
          + '{:9}'.format(
        metadata['DELETE_PRESSURE_REVERSAL_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[1])
          + '{:>9}'.format(str(cast_shift_o['cast' + str(cast_number)].shape[0]))
          + '{:>10}'.format(str(cast_wakeeffect['cast' + str(cast_number)].shape[0] -
                                list(cast_wakeeffect['cast' + str(cast_number)].isna().sum())[0])))
    print("        BINAVE   " + '{:7}'.format(str(1.0))
          + '{:11}'.format(metadata['BINAVE_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[0])
          + '{:9}'.format(metadata['BINAVE_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[1])
          + '{:>9}'.format(str(cast_wakeeffect['cast' + str(cast_number)].shape[0]))
          + '{:>10}'.format(str(cast_binned['cast' + str(cast_number)].shape[0])))
    print("        EDIT     " + '{:7}'.format(str(1.0))
          + '{:11}'.format(metadata['FINALEDIT_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[0])
          + '{:9}'.format(metadata['FINALEDIT_Time'].strftime("%Y/%m/%d %H:%M:%S.%f")[0:-7].split(" ")[1])
          + '{:>9}'.format(str(cast_binned['cast' + str(cast_number)].shape[0]))
          + '{:>10}'.format(str(cast_final['cast' + str(cast_number)].shape[0])))

    print("    $END")
    print(" $REMARKS")

    list_number = len(metadata['Processing_history'].split("|"))
    for i in range(0, list_number, 1):
        print("     " + metadata['Processing_history'].split("|")[i])
    print("$END")
    print()


def write_comments(oxy, metadata_dict, cast_d):
    cruise_ID = metadata["Mission"]
    print("*COMMENTS")
    print("    " + "-" * 85)
    print()
    print("    Data Processing Notes:")
    print("    " + "-" * 22)
    print("       " + "No calibration sampling was available.")
    print()
    print("       " + "For details on the processing see document: " + cruise_ID + "_Processing_Report.doc.")

    if oxy == 'yes':
        print("!--1--- --2--- ---3---- ---4---- ---5--- ---6--- ----7---- -8--")
        print("!Pressu Depth  Temperat Salinity Fluores Oxygen: Conductiv Numb")
        print("!re            ure               cence:  Dissolv ity       er_o")
        print("!                                URU     ed:               ~bin")
        print("!                                        Saturati          _rec")
        print("!                                        on:RBR            ords")
        print("!------ ------ -------- -------- ------- ------- --------- ----")
        print("*END OF HEADER")

    elif oxy == 'no':
        print("!--1--- --2--- ---3---- ---4---- ----5---- -6--")
        print("!Pressu Depth  Temperat Salinity Conductiv Numb")
        print("!re            ure               ity       er_o")
        print("!                                          ~bin")
        print("!                                          _rec")
        print("!                                          ords")
        print("!------ ------ -------- -------- --------- ----")
        print("*END OF HEADER")


def write_data(oxy, cast_data, cast_number, cast_d):
    # try:
    # check_for_oxy = cast_d.Oxygen.values
    if oxy == 'yes':
        for i in range(len(cast_data['cast' + str(cast_number)])):
            # print(cast_data['cast' + str(cast_number)]['Pressure'][i] + cast_data['cast' + str(cast_number)]['Depth'][i] + "  ")
            print('{:>7}'.format(cast_data['cast' + str(cast_number)].Pressure[i]) + " "
                  + '{:>6}'.format(cast_data['cast' + str(cast_number)].Depth[i]) + " "
                  + '{:>8}'.format(cast_data['cast' + str(cast_number)].Temperature[i]) + " "
                  + '{:>8}'.format(cast_data['cast' + str(cast_number)].Salinity[i]) + " "
                  + '{:>7}'.format(cast_data['cast' + str(cast_number)]['Fluorescence:URU'][i]) + " "
                  + '{:>7}'.format(cast_data['cast' + str(cast_number)]['Oxygen:Dissolved:Saturation:RBR'][i]) + " "
                  + '{:>9}'.format(cast_data['cast' + str(cast_number)]['Conductivity'][i]) + " "
                  + '{:>4}'.format(cast_data['cast' + str(cast_number)]['Number_of_bin_records'][i]) + " ")
    elif oxy == 'no':
        for i in range(len(cast_data['cast' + str(cast_number)])):
            # print(cast_data['cast' + str(cast_number)]['Pressure'][i] + cast_data['cast' + str(cast_number)]['Depth'][i] + "  ")
            print('{:>7}'.format(cast_data['cast' + str(cast_number)].Pressure[i]) + " "
                  + '{:>6}'.format(cast_data['cast' + str(cast_number)].Depth[i]) + " "
                  + '{:>8}'.format(cast_data['cast' + str(cast_number)].Temperature[i]) + " "
                  + '{:>8}'.format(cast_data['cast' + str(cast_number)].Salinity[i]) + " "
                  + '{:>9}'.format(cast_data['cast' + str(cast_number)]['Conductivity'][i]) + " "
                  + '{:>4}'.format(cast_data['cast' + str(cast_number)]['Number_of_bin_records'][i]) + " ")


def main_header(dest_dir, n_cast, meta_data, cast, cast_d, cast_d_clip, cast_d_filtered, cast_d_shift_c, cast_d_shift_o,
                cast_d_wakeeffect, cast_d_binned, cast_d_final, oxy):
    f_name = dest_dir.split("/")[-2]
    f_output = f_name.split("_")[0] + '-' + f'{n_cast:04}' + ".CTD"
    output = dest_dir + "CTD/" + f_output
    newnc_dir = '{}CTD/'.format(dest_dir)
    if not os.path.exists(newnc_dir):
        os.makedirs(newnc_dir)
    # Start
    # datetime object containing current date and time
    now = datetime.now()

    # dd/mm/YY H:M:S
    dt_string = now.strftime("%Y/%m/%d %H:%M:%S.%f")[0:-4]


    IOS_string = '*IOS HEADER VERSION 2.0      2020/03/01 2020/04/15 PYTHON'

    orig_stdout = sys.stdout
    file_handle = open(output, 'wt')
    try:
        sys.stdout = file_handle
        print("*" + dt_string)
        print(IOS_string)
        print()  # print("\n") pring("\n" * 40)
        write_file(n_cast, cast, cast_d_final, oxy, metadata_dict=meta_data)

        # def write_file(cast_number, cast_original, cast_final, metadata_dict):
        write_admin(metadata_dict=meta_data)
        write_location(n_cast, metadata_dict=metadata)
        write_instrument(metadata_dict=meta_data)
        write_history(cast_d, cast_d_clip, cast_d_filtered,
                     cast_d_shift_c, cast_d_shift_o, cast_d_wakeeffect,
                     cast_d_binned, cast_d_final, metadata_dict=meta_data, cast_number=n_cast)
        write_comments(oxy, metadata_dict=meta_data, cast_d=cast_d)
        write_data(oxy, cast_d_final, cast_number=n_cast, cast_d=cast_d)
        sys.stdout.flush()  # Recommended by Tom
    finally:
        sys.stdout = orig_stdout
    return os.path.abspath(output)


# main_header(dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2021-039_Python_Run/', n_cast = 1,
# meta_data = metadata, cast_d = cast_d, cast_d_clip = cast_d_clip, cast_d_filtered = cast_d_filtered,
# cast_d_shift_c = cast_d_shift_c, cast_d_shift_o = cast_d_shift_o,
# cast_d_wakeeffect = cast_d_wakeeffect, cast_d_binned = cast_d_binned, cast_d_final = cast_d_final)


# for i in range(1, len(cast)+1, 1):
# main_header(dest_dir='C:/Projects/RBR_CTD_Pycharm/2021-039_Python_Run/', n_cast=i,
# meta_data=metadata, cast_d=cast_d, cast_d_clip=cast_d_clip, cast_d_filtered=cast_d_filtered,
# cast_d_shift_c=cast_d_shift_c, cast_d_shift_o=cast_d_shift_o,
# cast_d_wakeeffect=cast_d_wakeeffect, cast_d_binned=cast_d_binned, cast_d_final=cast_d_final)


dest_dir = 'C:/Projects/RBR_CTD_Pycharm/2021-080/'  ####  DON"T FORGET TO CHANGE THE SHEETNAME TO Profile_annotation or openpyxl won't read it.
year = '2021'
cruise_number = '080'


def get_started(dest_dir):
    """ Start by opening the RSK files, find out how many channels and profiles there are

        Compare this to the cast list given by chief scientist

        prep metadata .csv

        have the header-merge.csv ready

        """

    files = os.listdir(dest_dir)  # list all the files in dest_dir
    files = list(filter(lambda f: f.endswith('.rsk'), files))  # keep the rsk files only
    n_files = len(files)  # get the number of files
    print(n_files)

    for k in range(0, n_files, 1):
        filename = str(dest_dir) + str(files[k])  # full path and name of .rsk file
        print(filename)
        rsk = pyrsktools.open(filename)  # load up an RSK

        # check the number of profiles
        n_profiles = len(list(rsk.profiles()))  # get the number of profiles recorded

        # last_profile = n_profiles[:-1] # get the last profile only
        last_profile = list(rsk.profiles())[:-1]
        print(n_profiles)
        print(rsk.samples)
        print(list(rsk.channels.items()))



get_started(dest_dir='C:/Projects/RBR_CTD_Pycharm/2021-080/')

def first_step(dest_dir, year, cruise_number, input_ext, rsk_file, event_start,all_last, left_lon, right_lon, bot_lat, top_lat, file_option): # file_options: single, multi, rsk_excel

    """Choose how to export the csv files from the rsk files

     Plot cruise, plot pre-processing plots, determine need for zero-order holds correction

     if multi_file choose an rsk_file for metadata"""

    print('checking need for zero-order holds correction')
    if file_option == 'single':
        EXPORT_FILES(dest_dir, rsk_file, year, cruise_number, event_start)
    elif file_option == 'multi':
        EXPORT_MULTIFILES(dest_dir, year, cruise_number, all_last, event_start)
    elif file_option == 'from_excel':
        ### input file =  # not needed, filtered to keep only .xls
        READ_EXCELrsk(dest_dir, year, cruise_number, event_start, all_last)
    MERGE_FILES(dest_dir, year, cruise_number)
    print('files merged')
    ADD_6LINEHEADER_2(dest_dir, year , cruise_number)
    Plot_Track_Location(dest_dir, year, cruise_number, left_lon, right_lon, bot_lat, top_lat)
    First_Plots(year, cruise_number, dest_dir, input_ext)
    PLOT_PRESSURE_DIFF(dest_dir, year, cruise_number, input_ext)

# input_ext '_CTD_DATA-6linehdr.csv' for first steps, all_last 'ALL' or 'LAST', event_start comes from cast list or cruise log
# file_option 'single', 'multi', or 'from_excel' !!! if 'from_excel enter the file in the function

#first_step(dest_dir=dest_dir, year=year, cruise_number=cruise_number, input_ext='_CTD_DATA-6linehdr.csv', rsk_file='204694_20211216_1628.rsk', event_start=1,
          # all_last='ALL', left_lon=-129, right_lon=-127, bot_lat=49.5, top_lat=51.5, file_option='from_excel')


metadata = CREATE_META_DICT(dest_dir=dest_dir, file='204694_20211216_1628.rsk', year=year,
                            cruise_number=cruise_number)


def second_step(dest_dir, year, cruise_number, metadata_dict, zoh, correction_value, oxy,
                n_cast):  # zoh is zero-order holds correction 'yes', or 'no'

    """Make corrections for zero-order holds and Pressure if needed"""
    #initialize casts for if statement
    cast, cast_d, cast_u = 0, 0, 0
    cast_pc, cast_d_pc, cast_u_pc = 0, 0, 0

    if zoh == 'yes':
        CORRECT_HOLD(dest_dir, year, cruise_number, metadata_dict)
        # check the plot then save it as Fig_3
        check_plot = PLOT_PRESSURE_DIFF(dest_dir, year, cruise_number, input_ext='_CTD_DATA-6linehdr_corr_hold.csv')
        # re_create cast variables
        t_cast, t_cast_d, t_cast_u = CREATE_CAST_VARIABLES(year, cruise_number, dest_dir,
                                                           input_ext='_CTD_DATA-6linehdr_corr_hold.csv')
        t_cast_pc, t_cast_d_pc, t_cast_u_pc = CALIB(t_cast, t_cast_d, t_cast_u, correction_value, metadata_dict,
                                                    zoh='yes')  # 0 if no neg pressures
        cast, cast_d, cast_u = t_cast, t_cast_d, t_cast_u
        cast_pc, cast_d_pc, cast_u_pc = t_cast_pc, t_cast_d_pc, t_cast_u_pc  # probably don't need the t_....
        print('using zero-order holds corrected variables')
        print('The following correction value has been applied to Pressure:')
        print(correction_value)
    elif zoh == 'no':
        cast, cast_d, cast_u = CREATE_CAST_VARIABLES(year, cruise_number, dest_dir, input_ext='_CTD_DATA-6linehdr.csv')
        cast_pc, cast_d_pc, cast_u_pc = CALIB(cast, cast_d, cast_u, correction_value, metadata_dict,
                                              zoh='no')  # 0 if no neg pressures
        print('using original variables')
        print('The following correction value has been applied to Pressure:')
        print(correction_value)

    pre_plots = First_Plots(year, cruise_number, dest_dir, input_ext='_CTD_DATA-6linehdr_corr_hold.csv')
    print('finished plotting first plots')


    vars = list(dict.fromkeys(cast['cast1']))

    # clip the casts

    cast_d_clip = CLIP_DOWNCAST(cast_d_pc, metadata, limit_drop=0.02)
    cast_u_clip = CLIP_UPCAST(cast_u_pc, metadata, limit_rise=-0.02)

    check_clip = plot_clip(cast, cast_d_clip, cast_d_pc)

    cast_d_filtered, cast_u_filtered = FILTER(cast_d_clip, cast_u_clip, window_width=6, sample_rate=8,
                                              time_constant=1 / 8, filter_type=1,
                                              metadata_dict=metadata)  # n = 5 should be good.

    check_filter = plot_filter(cast_d_filtered, cast_d_clip)

    cast_d_shift_c, cast_u_shift_c = SHIFT_CONDUCTIVITY(cast_d_filtered, cast_u_filtered, shifted_scan_number=2,
                                                        metadata_dict=metadata)  # delay conductivity data by 2 scans

    check_shift_c = plot_shift_c(cast_d_shift_c, cast_d_filtered)

    if oxy == 'yes':
        cast_d_shift_o, cast_u_shift_o = SHIFT_OXYGEN(cast_d_shift_c, cast_u_shift_c, shifted_scan_number=-11,
                                                      metadata_dict=metadata)  # advance oxygen data by 11 scans
        check_shift_o = plot_shift_o(cast_d_shift_o, cast_d_shift_c)
    elif oxy == 'no':
        cast_d_shift_o, cast_u_shift_o = cast_d_shift_c, cast_u_shift_c

    cast_d_wakeeffect, cast_u_wakeeffect = DELETE_PRESSURE_REVERSAL(cast_d_shift_o, cast_u_shift_o,
                                                                    metadata_dict=metadata)

    check_process = plot_processed(cast_d_wakeeffect, cast_d_shift_o, cast_d, cast)

    cast_d_binned, cast_u_binned = BINAVE(cast_d_wakeeffect, cast_u_wakeeffect, interval=1, metadata_dict=metadata)

    cast_d_final = FINAL_EDIT(cast_d_binned, oxy, metadata_dict)

    main_header(dest_dir, n_cast,
                metadata, cast, cast_d, cast_d_clip, cast_d_filtered,
                cast_d_shift_c, cast_d_shift_o,
                cast_d_wakeeffect, cast_d_binned, cast_d_final, oxy)

    for i in range(1, len(cast), 1):
        main_header(dest_dir, i + 1,
                    metadata, cast, cast_d, cast_d_clip, cast_d_filtered,
                    cast_d_shift_c, cast_d_shift_o,
                    cast_d_wakeeffect, cast_d_binned, cast_d_final, oxy)


second_step(dest_dir, year, cruise_number, metadata_dict=metadata, zoh='yes', correction_value=0, oxy='yes',
            n_cast=1)  # oxy = 'no' or 'yes', zoh = look at plots from first step and determine 'yes' or 'no'.