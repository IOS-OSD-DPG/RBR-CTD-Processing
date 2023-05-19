"""
author: Lu Guan
date: Oct. 06, 2020
about: This script is for processing RBR CTD data and producing .ctd files in IOS Header format.

Modified July 2021 - September 2021 by Samantha Huntington

Modified Jan. 2023 - ___ by Hana Hourston @hhourston
"""
# globals().clear()
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter, LatitudeLocator
import sys
import os
import pyrsktools
import numpy as np
import pandas as pd
from mpl_toolkits.basemap import Basemap
from copy import deepcopy  # copy,
from scipy import signal
import gsw
from matplotlib import pyplot as plt
from datetime import datetime
import random
from ocean_data_parser.convert.oxygen import O2stoO2c
from seawater import eos80

# import warnings
# import itertools
# from datetime import datetime, timezone
# import pyproj
# import glob
# from datetime import timedelta
# from decimal import Decimal
# import openpyxl
# from openpyxl import load_workbook


# Global variables
VARIABLES_POSSIBLE = [
    "Salinity",
    "Temperature",
    "Conductivity",
    "Oxygen",
    "Fluorescence",
    "Oxygen_mL_L",
    "Oxygen_umol_kg",
]
VARIABLE_UNITS = ["PSS-78", "C", "mS/cm", "%", "ug/L", "mL/L", "umol/kg"]
VARIABLE_COLOURS = ["b", "r", "goldenrod", "grey", "g", "grey", "grey"]


# ------- step 1: .rsk file in EPdesktop structure - Export profile data to .csv files from .rsk files-------

# run function to export files
# sh - call function at bottom of function


# def EXPORT_FILES(dest_dir, file, year, cruise_number, event_start):
#     """
#     Read in a rsk file and output in csv format
#     Inputs:
#         - folder, file, year, cruise_number: rsk-format file containing raw RBR data
#         - event_start: taken from ctd log or cruise log
#     Outputs:
#         - csv file: csv files containing the profile data
#     """
#
#     # function to export .csv files from .rsk file.
#     filename = str(dest_dir) + str(file)  # full path and name of .rsk file
#     path_slash_type = '/' if '/' in filename else '\\'
#     rsk = pyrsktools.open(filename)  # load up an RSK
#
#     # check the number of profiles
#     n_profiles = len(list(rsk.profiles()))  # get the number of profiles recorded
#
#     ctd_data = pd.DataFrame()  # set an empty pandas dataframe to store data
#
#     # print(list(rsk.channels.items()))
#     # export .csv file for each profile
#     for i in range(0, n_profiles, 1):
#
#         downcast = list(itertools.islice(rsk.casts(pyrsktools.Region.CAST_DOWN), i, i + 1))[
#             0].npsamples()  # separate samples for each downcast file
#         downcast_dataframe = pd.DataFrame(data=downcast,
#                                           columns=downcast.dtype.names)  # convert data into pandas data frame
#
#         upcast = list(itertools.islice(rsk.casts(pyrsktools.Region.CAST_UP), i, i + 1))[0].npsamples()
#         upcast_dataframe = pd.DataFrame(data=upcast, columns=upcast.dtype.names)
#         # print(downcast[0])
#
#         column_names = list(downcast.dtype.names)  # read the column names
#         # print(list(rsk.channels))
#         # print(column_names)
#         # print(len(column_names))
#         col_range = len(column_names)  # added this to get them all
#         column_names[0] = 'Time(yyyy-mm-dd HH:MM:ss.FFF)'  # update time
#         for j in range(1, col_range, 1):  # update the column names
#             column_names[j] = column_names[j][0: -3] + "(" + list(rsk.channels.items())[j - 1][1][4] + ")"
#
#         downcast_dataframe.columns = column_names  # update column names in downcast data frame
#         downcast_dataframe["Cast_direction"] = "d"  # add a column for cast direction
#         downcast_dataframe["Event"] = i + event_start  # add a column for event number - get event_start from Logs
#         upcast_dataframe.columns = column_names
#         upcast_dataframe["Cast_direction"] = "u"
#         upcast_dataframe["Event"] = i + event_start
#         # downcast_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(i + event_start).zfill(
#         #     4) + "_DOWNCAST.csv"  # downcast file name
#         # upcast_name = filename.split("/")[-1][0:-4].upper() + "_profile" + str(i + event_start).zfill(
#         #     4) + "_UPCAST.csv"  # upcast file name
#         profile_name = filename.split(path_slash_type)[-1][0:-4].upper() + "_profile" + str(i + event_start).zfill(
#             4) + ".csv"  # profile name
#         # combine downcast and upcast into one profile
#         profile_data = pd.concat([downcast_dataframe, upcast_dataframe])
#         ctd_data = ctd_data.append(profile_data, ignore_index=True)  # combine profiles into one file
#
#         # downcast_dataframe.to_csv(folder + downcast_name)
#         # upcast_dataframe.to_csv(folder + upcast_name)
#         profile_data.to_csv(dest_dir + profile_name)  # export each profile
#
#     output_filename = year + "-" + cruise_number + "_CTD_DATA.csv"  # all data file name
#     ctd_data.to_csv(dest_dir + output_filename)  # export all data in one .csv file
#
#     return


# -------------------  Exploring if we can input mulitple .rsk files----------------------


def READ_RSK(
    dest_dir: str,
    year: str,
    cruise_number: str,
    skipcasts,
    rsk_time1=None,
    rsk_time2=None,
) -> None:
    """
    Formerly EXPORT_MULTIFILES()
    Replace all_last parameter to skipcasts: list-like or int For cases when there
    are data from previous cruise in the file, or initial test casts, or soaks

    Read in a directory of rsk files and output in csv format
    Inputs:
        - dest_dir
        - year
        - cruise_number:
        - skipcasts: number of casts to skip over in an rsk file when writing data to output format.
            Input format as either an integer or as a list-like object with one integer per excel file
            representing the number of initial casts to skip in each excel file.
        - rsk_time1, rsk_time2
    Outputs:
        - csv file: csv files containing the profile data
    """
    files = os.listdir(dest_dir)  # list all the files in dest_dir
    files = list(filter(lambda f: f.endswith(".rsk"), files))  # keep the rsk files only
    n_files = len(files)  # get the number of files

    # current_profile = event_start
    header_merge_file = os.path.join(
        dest_dir, f"{year}-{cruise_number}_header-merge.csv"
    )
    header_merge_df = pd.read_csv(header_merge_file)
    # Get event numbers from the header-merge.csv file
    header_event_no = header_merge_df.loc[:, "LOC:Event Number"].to_numpy()
    event_number_idx = (
        0  # initialize counter to go through the event numbers in header_merge_df
    )

    # Iterate through all the rsk files that were found
    for k in range(n_files):
        print(files[k])
        # Open the rsk file and read the data within it
        filename = os.path.join(
            str(dest_dir), str(files[k])
        )  # full path and name of .rsk file
        # readHiddenChannels=True does not reveal the derived variables
        rsk = pyrsktools.RSK(filename, readHiddenChannels=False)  # load up an RSK
        rsk.open()
        rsk.readdata(t1=rsk_time1, t2=rsk_time2)
        # rsk.data now returns data in a structured array format

        # Compute the derived channels
        rsk.derivesalinity()
        rsk.deriveseapressure()
        rsk.derivedepth()
        # rsk.deriveO2()  # Derive later

        # # Obtain the number of profiles in the rsk file
        # try:
        #     rsk_num_profiles = len(rsk.getprofilesindices(direction="down"))
        # except AttributeError:
        #     rsk.computeprofiles()  # This should fix the problem
        #     rsk_num_profiles = len(rsk.getprofilesindices(direction="down"))

        # # Add a check on the number of profiles in the rsk file
        # # compare this number to the number from the rbr/cruise log
        # if rsk_num_profiles != log_num_profiles and all_last == 'ALL':
        #     warnings.warn(f'Number of rsk profiles does not match input number of profiles, '
        #                   f'{rsk_num_profiles} != {log_num_profiles}...')
        # # If the number of profiles does not match the number of profiles in the log,
        # # compute the profiles (distinguish the different profiles) based on
        # # pressure and conductivity data
        # # Choose the best pressure and conductivity thresholds?
        # pressureThreshold = (
        #                             max(rsk.data['sea_pressure']) -
        #                             min(rsk.data['sea_pressure'])
        #                     ) * 1 / 4
        # rsk.computeprofiles(pressureThreshold=pressureThreshold,
        #                     conductivityThreshold=0.05)
        # if rsk_num_profiles != len(rsk.getprofilesindices(direction="down")):
        #     print(f'Recomputed number of rsk profiles does not match input '
        #           f'number of profiles, {rsk_num_profiles} != {log_num_profiles}. '
        #           f'Ending process!')
        #     return
        # profileIndices = rsk.getprofilesindices()  # Returns a list of lists
        try:
            downcastIndices = rsk.getprofilesindices(direction="down")
            upcastIndices = rsk.getprofilesindices(direction="up")
        except AttributeError:
            rsk.computeprofiles()  # This should fix the problem
            downcastIndices = rsk.getprofilesindices(direction="down")
            upcastIndices = rsk.getprofilesindices(direction="up")

        # check the number of profiles
        n_profiles = len(downcastIndices)  # get the number of profiles recorded

        if type(skipcasts) == int:
            profile_range = range(skipcasts, n_profiles)
        elif len(skipcasts) == len(header_event_no):
            # Check that for each event there is a corresponding number of casts to skip
            profile_range = range(skipcasts[k], n_profiles)
        else:
            print(
                f"Invalid skipcasts value, {skipcasts}. Must be int or same length as number of rsk files."
            )
            return

        if len(profile_range) == 0:
            print(
                "Warning: Number of casts to skip in file",
                files[k],
                "equals the number of casts in the file",
            )

        # Iterate through the selected profiles in one rsk file
        for i in profile_range:
            downcast_dataframe = pd.DataFrame(rsk.data).loc[downcastIndices[i]]
            upcast_dataframe = pd.DataFrame(rsk.data).loc[upcastIndices[i]]
            # downcast = list(itertools.islice(
            #     rsk.casts(pyrsktools.Region.CAST_DOWN), i, i + 1))[
            #     0].npsamples()  # separate samples for each downcast file
            # downcast_dataframe = pd.DataFrame(data=downcast,
            #                                   columns=downcast.dtype.names)  # convert data into pandas data frame
            #
            # upcast = list(itertools.islice(
            #     rsk.casts(pyrsktools.Region.CAST_UP), i, i + 1))[0].npsamples()
            # upcast_dataframe = pd.DataFrame(data=upcast, columns=upcast.dtype.names)

            column_names = list(downcast_dataframe.columns)  # read the column names
            channel_units = [chan.units for chan in rsk.channels]

            # update the column names to include units
            column_names[0] = "Time(yyyy-mm-dd HH:MM:ss.FFF)"  # update time first
            for j in range(1, len(column_names)):
                column_names[j] = column_names[j] + "(" + channel_units[j - 1] + ")"

            # update column names in downcast and upcast data frames
            downcast_dataframe.columns = column_names
            upcast_dataframe.columns = column_names

            # add a column for cast direction
            downcast_dataframe["Cast_direction"] = "d"
            downcast_dataframe["Event"] = header_event_no[event_number_idx]

            upcast_dataframe["Cast_direction"] = "u"
            upcast_dataframe["Event"] = header_event_no[event_number_idx]

            # combine downcast and upcast into one profile
            profile_data = pd.concat([downcast_dataframe, upcast_dataframe])

            profile_name = "{}_profile{}.csv".format(
                os.path.basename(filename)[:-4],
                str(header_event_no[event_number_idx]).zfill(4),
            )

            profile_data.to_csv(
                dest_dir + profile_name, index=False
            )  # export each profile

            # Update counter
            event_number_idx += 1
    return


def READ_EXCELrsk(dest_dir: str, year: str, cruise_number: str, skipcasts) -> None:
    """
    Function to read in an excel (.xlsx) file exported from RUSKIN software using
    the rbr .rsk file. An alternative to reading data directly from the .rsk files.
    inputs:
        - dest_dir: working & destination directory
        - year: year of input cruise
        - cruise_number: cruise number of input cruise
        - skipcasts: number of casts to skip over in an excel file when writing data to output format.
            Input format as either an integer or as a list-like object with one integer per excel file
            representing the number of initial casts to skip in each excel file.
    outputs:
        - one csv file per RBR cast, containing the data from that cast
    """

    files = os.listdir(dest_dir)  # list all the files in dest_dir
    files = list(
        filter(lambda f: f.endswith(".xlsx"), files)
    )  # keep the rsk xlsx files only (make sure no other xlsx)

    header_merge_file = os.path.join(
        dest_dir, f"{year}-{cruise_number}_header-merge.csv"
    )
    header_merge_df = pd.read_csv(header_merge_file)
    header_event_no = header_merge_df.loc[:, "LOC:Event Number"].to_numpy()
    event_number_idx = (
        0  # initialize index counter to go through the event numbers in header_merge_df
    )

    # Iterate through all available excel files
    for k in range(len(files)):
        filename = os.path.join(dest_dir, files[k])  # full path and name of .rsk file

        # extract a dataframe from the excel sheets
        # hourstonh 2023-01-25: xlrd package does not read xlsx files any more; See release notes
        # https://groups.google.com/g/python-excel/c/IRa8IWq_4zk/m/Af8-hrRnAgAJ?pli=1
        df1 = pd.read_excel(
            filename, sheet_name="Data", skiprows=[0], engine="openpyxl"
        )
        try:
            df2 = pd.read_excel(
                filename,
                sheet_name="Profile_annotation",
                skiprows=[0],
                engine="openpyxl",
            )
        except ValueError:
            print(
                'User must change "Profile annotation" sheet name to "Profile_annotation"'
            )
            return
        # df3 = pd.read_excel(filename, sheet_name='Metadata',
        #                     skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        #                     usecols=[2], engine='openpyxl')

        # Convert time data to pandas datetime format
        df1["Time"] = pd.to_datetime(df1["Time"])
        df2["Time 1"] = pd.to_datetime(df2["Time 1"])
        df2["Time 2"] = pd.to_datetime(df2["Time 2"])

        down_times = pd.DataFrame()
        up_times = pd.DataFrame()

        # find the start and end times for each profile
        down_times["Start"] = df2["Time 1"][1::3]
        down_times["End"] = df2["Time 2"][1::3]
        down_times.index = range(len(down_times))
        up_times["Start"] = df2["Time 1"][2::3]
        up_times["End"] = df2["Time 2"][2::3]
        up_times.index = range(len(up_times))

        total_casts = len(
            list(down_times["Start"])
        )  # Number of casts in the input file

        # Iterate through the profiles in a single file
        if type(skipcasts) == int:
            profile_range = range(skipcasts, total_casts)
        elif len(skipcasts) == len(header_event_no):
            profile_range = range(skipcasts[k], total_casts)
        else:
            print(
                f"Invalid skipcasts value, {skipcasts}. Must be int or same length as number of rsk files."
            )
            return

        if len(profile_range) == 0:
            print(
                "Warning: Number of casts to skip in file",
                files[k],
                "equals the number of casts in the file",
            )

        # Iterate through the selected profiles
        for i in profile_range:
            down_start_time = down_times["Start"][i]
            down_end_time = down_times["End"][i]
            up_start_time = up_times["Start"][i]
            up_end_time = up_times["End"][i]

            # extract data for each profile - using start and end times
            downcast = df1.loc[
                (df1["Time"] > down_start_time) & (df1["Time"] <= down_end_time)
            ]
            downcast["Cast_direction"] = "d"
            # Add event numbers
            downcast["Event"] = header_event_no[event_number_idx]  # current_profile
            upcast = df1.loc[
                (df1["Time"] > up_start_time) & (df1["Time"] <= up_end_time)
            ]
            upcast["Cast_direction"] = "u"
            upcast["Event"] = header_event_no[event_number_idx]  # current_profile
            # combine downcast and upcast into one profile
            profile_data = pd.concat([downcast, upcast])

            profile_name = "{}_profile{}.csv".format(
                os.path.basename(filename)[:-5],
                str(header_event_no[event_number_idx]).zfill(4),
            )
            print(profile_name)

            # Change column name for time
            profile_data.rename(
                {"Time": "Time(yyyy-mm-dd HH:MM:ss.FFF)"}, axis=1, inplace=True
            )

            profile_data.to_csv(os.path.join(dest_dir, profile_name), index=False)

            # Update event number index counter
            event_number_idx += 1

    return


def MERGE_FILES(dest_dir: str, year: str, cruise_number: str) -> None:
    """
    Read in multiple csv file and output one csv format
    Inputs:
        - dest_dir
        - year
        - cruise_number
    Outputs:
        - csv file: csv files containing the profile data
    """
    files = os.listdir(dest_dir)  # list all the files in dest_dir
    # keep the profile csv files only
    files = list(filter(lambda fname: "profile" in fname, files))
    # event_csv = pd.read_csv(dest_dir + year + '-'
    #                         + cruise_number + '_header-merge.csv')
    files.sort(
        key=lambda x: int(x[-8:-4])
    )  # reorder the files according to profile number

    ctd_data = pd.DataFrame()  # set an empty pandas dataframe to store data

    # Iterate through all profiles and concatenate them together to form one dataset
    for i in range(len(files)):
        input_filename = str(dest_dir) + str(files[i])
        # data = pd.read_csv(input_filename, sep=',', skiprows=range(0, 20),
        #                    encoding= 'unicode_escape') #original by Lu
        data = pd.read_csv(input_filename, sep=",", encoding="unicode_escape")
        # if event_from == 'filename':
        data.loc[:, "Event"] = np.repeat(int(files[i][-8:-4]), len(data))
        # elif event_from == 'header-merge':
        # Use this if you want non-sequential events, but won't work auto processing (cast i+1 etc...)
        # data.loc[:, 'Event'] = np.repeat(event_csv.loc[i, 'LOC:Event Number'], len(data))
        # data['Event'] = int(files[i][-8:-4])
        # data['Time'] = pd.to_datetime(data['Time'])        #might not be needed
        ctd_data = pd.concat([ctd_data, data], ignore_index=True)

    if ctd_data.columns.to_list()[0] == "//Time(yyyy-mm-dd HH:MM:ss.FFF)":
        ctd_data.rename(
            {"//Time(yyyy-mm-dd HH:MM:ss.FFF)": "Time(yyyy-mm-dd HH:MM:ss.FFF)"},
            axis=1,
            inplace=True,
        )

    # Output merged dataset in csv format
    output_filename = year + "-" + cruise_number + "_CTD_DATA.csv"
    ctd_data.to_csv(dest_dir + output_filename, index=False)
    return


def CREATE_META_DICT(
        dest_dir: str,
        rsk_file: str,
        year: str,
        cruise_number: str,
        rsk_time1=None,
        rsk_time2=None,
) -> dict:
    """
    Read in a csv file and output a metadata dictionary.
    Use a .rsk file to get channel details.
    Inputs:
        - folder, file, year, cruise_number:
        - rsk_file: rsk-format file containing raw RBR data & csv file containing metadata
    Outputs:
        - metadata dictionary
    """
    # Initialize metadata dictionary
    meta_dict = {}
    # Read in a .rsk file just to get the channel details
    rsk_full_path = os.path.join(
        str(dest_dir), str(rsk_file)
    )  # full path and name of a .rsk file
    rsk = pyrsktools.RSK(rsk_full_path, readHiddenChannels=False)  # load up an RSK
    rsk.open()
    rsk.readdata(t1=rsk_time1, t2=rsk_time2)

    # Compute the derived channels so as to make number of channels accurate
    rsk.derivesalinity()
    rsk.deriveseapressure()
    rsk.derivedepth()

    header_input_name = str(year) + "-" + str(cruise_number) + "_header-merge.csv"
    header_input_filename = dest_dir + header_input_name
    header = pd.read_csv(header_input_filename, header=0)

    # get the time interval for the IOS Header (sampling period)
    time_input_name = str(year) + "-" + str(cruise_number) + "_CTD_DATA.csv"
    time_input_filename = dest_dir + time_input_name
    time_input = pd.read_csv(time_input_filename)
    time_interval = pd.to_datetime(
        time_input["Time(yyyy-mm-dd HH:MM:ss.FFF)"][2]
    ) - pd.to_datetime(time_input["Time(yyyy-mm-dd HH:MM:ss.FFF)"][1])
    time_interval = str(time_interval)
    time_interval = time_interval[-8:-3]

    # Metadata file
    csv_input_name = str(year) + "-" + str(cruise_number) + "_METADATA.csv"
    csv_input_filename = dest_dir + csv_input_name
    meta_csv = pd.read_csv(csv_input_filename)

    # Fill in metadata values
    meta_dict["Processing_Start_time"] = datetime.now()
    meta_dict["Instrument_information"] = rsk.instrument
    meta_dict["Sampling_Interval"] = time_interval
    # print("time_interval", time_interval)
    # meta_dict['RSK_filename'] = rsk.name
    meta_dict["RSK_filename"] = meta_csv["Value"][
                                    meta_csv["Name"] == "RSK_filename"
                                    ].values[0:]  # if more than one (list??)
    # meta_dict['Channels'] = list(rsk.channels.keys())
    meta_dict["Channels"] = rsk.channelNames
    # meta_dict['Channel_details'] = list(rsk.channels.items())
    meta_dict["Channel_details"] = rsk.channels
    meta_dict["Number_of_channels"] = len(rsk.channels)
    meta_dict["Location"] = header

    for key in [
        "number_of_profiles",
        "Data_description",
        "Final_file_type",
        "Mission",
        "Agency",
        "Country",
        "Project",
        "Scientist",
        "Platform",
        "Instrument_Model",
        "Serial_number",
        "Instrument_type",
    ]:
        meta_dict[key] = meta_csv["Value"][meta_csv["Name"] == key].values[0]

    # Fill out the rest of the METADATA.csv file and export it to csv as a record for users
    meta_csv.loc[meta_csv["Name"] == "Instrument_information", "Value"] = str(
        meta_dict["Instrument_information"]
    )
    meta_csv.loc[meta_csv["Name"] == "Channels", "Value"] = str(meta_dict["Channels"])
    meta_csv.loc[meta_csv["Name"] == "Channel_details", "Value"] = str(
        meta_dict["Channel_details"]
    )
    meta_csv.loc[meta_csv["Name"] == "Number_of_channels", "Value"] = str(
        meta_dict["Number_of_channels"]
    )

    updated_meta_filename = csv_input_filename.replace(".csv", "_FILLED.csv")
    meta_csv.to_csv(updated_meta_filename, index=False)

    return meta_dict


def ADD_6LINEHEADER_2(dest_dir: str, year: str, cruise_number: str) -> None:
    """
    Read in a csv file and output in csv format for use in IOSShell
    Filter through the csv and remove un-needed columns.
    Inputs:
        - folder, file, year, cruise: csv-format file containing raw RBR CTD data
        exported from rsk file
    Outputs:
        - csv file: csv files containing 6-header line for IOSShell
    """
    # Add six-line header to the .csv file.
    # This file could be used for data processing via IOSShell
    input_name = str(year) + "-" + str(cruise_number) + "_CTD_DATA.csv"
    output_name = str(year) + "-" + str(cruise_number) + "_CTD_DATA-6linehdr.csv"
    input_filename = os.path.join(dest_dir, input_name)
    ctd_data = pd.read_csv(input_filename, header=0)

    ctd_data["Time(yyyy-mm-dd HH:MM:ss.FFF)"] = ctd_data[
        "Time(yyyy-mm-dd HH:MM:ss.FFF)"
    ].str[:19]
    ctd_data["Date"] = pd.to_datetime(
        ctd_data["Time(yyyy-mm-dd HH:MM:ss.FFF)"]
    )  # add new column of Date
    ctd_data["Date"] = [d.date() for d in ctd_data["Date"]]
    ctd_data["TIME:UTC"] = pd.to_datetime(
        ctd_data["Time(yyyy-mm-dd HH:MM:ss.FFF)"]
    )  # add new column of time
    ctd_data["TIME:UTC"] = [d.time() for d in ctd_data["TIME:UTC"]]
    ctd_data["Date"] = pd.to_datetime(ctd_data["Date"], format="%Y-%m-%d").dt.strftime(
        "%d/%m/%Y"
    )

    # make a list of columns to drop
    # drop_list = ['Temperature .1(Â°C)', 'temperature1(Ã‚Â°C)', 'Temperature .1', 'Temperature.1',
    #              'Speed of sound ', 'Speed of sound', 'speed_of_sound(m/s)', 'speedofsound(m/s)',
    #              'specificconductivity(ÂµS/cm)', 'specific_conductivity(ÂµS/cm)',
    #              'Specific conductivity ', 'Specific conductivity',
    #              'Density anomaly', 'Density anomaly ',
    #              'Dissolved Oâ,,Concentration', 'Dissolved OÃ¢Â‚Â‚ concentration ',
    #              'Dissolved OÃ¢Â‚Â‚ concentration', 'Dissolved O₂ concentration',
    #              'Dissolved O2 concentration', 'Dissolved Oâ\x82\x82 concentration ',
    #              'Dissolved Oâ\x82\x82 concentration',
    #              'Turbidity', 'turbidity(NTU)',
    #              'Time(yyyy-mm-dd HH:MM:ss.FFF)', 'Unnamed: 0']

    # Rules to replace list of columns to drop; may not be full-proof
    is_temperature_secondary = (
        lambda name: "temperature" in name.lower() and "1" in name
    )
    is_speedofsound = lambda name: "speed" in name.lower() and "sound" in name.lower()
    is_specificconductivity = (
        lambda name: "specific" in name.lower() and "conductivity" in name.lower()
    )
    is_densityanomaly = (
        lambda name: "density" in name.lower() and "anomaly" in name.lower()
    )
    is_dissolvedO2concentration = (
        lambda name: "dissolved" in name.lower()
        and "O" in name
        and "concentration" in name.lower()
    )
    is_turbidity = lambda name: "turbidity" in name.lower()

    drop_list = ["Time(yyyy-mm-dd HH:MM:ss.FFF)", "Unnamed: 0"]
    for col in ctd_data.columns:
        if any(
            [
                is_temperature_secondary(col),
                is_speedofsound(col),
                is_specificconductivity(col),
                is_densityanomaly(col),
                is_dissolvedO2concentration(col),
                is_turbidity(col),
            ]
        ):
            drop_list.append(col)

    # Drop first indexing row and Time(HH....) plus everything else in drop_list if present
    # Ignore KeyError if columns in drop_list do not exist in the dataframe
    ctd_data.drop(columns=drop_list, inplace=True, errors="ignore")
    # ctd_data.reset_index(drop=True, inplace=True)

    col_list = ctd_data.columns.tolist()
    print(col_list)

    # set empty column names
    # dict.fromkeys(ctd_data.columns, np.arange(len(ctd_data.columns)))
    column_names = {
        old: new for old, new in zip(col_list, np.arange(len(ctd_data.columns)))
    }
    ctd_data.rename(columns=column_names, inplace=True)  # remove column names

    # append header information into the empty lists
    channel_list = ["Y", "Y", "N", "Y", "Y", "Y", "Y", "Y", "N", "Y", "Y", "Y"]
    index_list = [
        "Conductivity",
        "Temperature",
        "Pressure_Air",
        "Fluorescence",
        "Oxygen:Dissolved:Saturation",
        "Pressure",
        "Depth",
        "Salinity:CTD",
        "Cast_direction",
        "Event_number",
        "Date",
        "TIME:UTC",
    ]
    unit_list = [
        "mS/cm",
        "deg C(ITS90)",
        "decibar",
        "mg/m^3",
        "%",
        "decibar",
        "meters",
        "PSS-78",
        "n/a",
        "n/a",
        "n/a",
        "n/a",
    ]
    input_format_list = [
        "R4",
        "R4",
        "R4",
        "R4",
        "R4",
        "R4",
        "R4",
        "R4",
        " ",
        "I4",
        "D:dd/mm/YYYY",
        "T:HH:MM:SS",
    ]
    output_format_list = [
        "R4:F11.4",
        "R4:F9.4",
        "R4:F7.1",
        "R4:F8.3",
        "R4:F11.4",
        "R4:F7.1",
        "R4:F7.1",
        "R4:F9.4",
        " ",
        "I:I4",
        "D:YYYY/mm/dd",
        "T:HH:MM:SS",
    ]
    na_value_list = [
        "-99",
        "-99",
        "-99",
        "-99",
        "-99",
        "-99",
        "-99",
        "-99",
        "",
        "-99",
        "",
        "",
    ]

    # # Need to add to this list of lists as more spellings come up from Ruskin
    # channel_spellings = [['conductivity(mS/cm)', 'Conductivity ', 'Conductivity'],
    #                      ['temperature(Â°C)', 'Temperature ', 'Temperature'],
    #                      ['pressure(dbar)', 'Pressure ', 'Pressure'],
    #                      ['chlorophyll(Âµg/l)', 'Chlorophyll a ', 'Chlorophyll a',
    #                       'chlorophyll_a(Âµg/l)'],
    #                      ['oxygensaturation(%)', 'Dissolved OÃ¢Â‚Â‚ saturation ',
    #                       'Dissolved Oâ\x82\x82 saturation ', 'Dissolved O2 saturation',
    #                       'Dissolved OÃ¢Â‚Â‚ saturation', 'Dissolved Oâ\x82\x82 saturation',
    #                       'dissolved_o2_saturation(%)'],
    #                      ['seapressure(dbar)', 'Sea pressure', 'Sea pressure ',
    #                       'sea_pressure(dbar)'],
    #                      ['depth(m)', 'Depth ', 'Depth'],
    #                      ['salinity(PSU)', 'Salinity ', 'Salinity'],
    #                      ['Cast_direction'],
    #                      ['Event'],
    #                      ['Date'],
    #                      ['TIME:UTC']]

    is_conductivity = lambda name: "conductivity" in name.lower()
    is_temperature_primary = lambda name: "temperature" in name.lower()
    is_pressure = lambda name: "pressure" in name.lower() and "sea" not in name.lower()
    is_chlorophyll = lambda name: "chlorophyll" in name.lower()
    is_O2saturation = (
        lambda name: "saturation" in name.lower()
        and "o" in name.replace("saturation", "").lower()
    )
    is_seapressure = lambda name: "pressure" in name.lower() and "sea" in name.lower()
    is_depth = lambda name: "depth" in name.lower()
    is_salinity = lambda name: "salinity" in name.lower()
    is_castdirection = lambda name: name == "Cast_direction"
    is_event = lambda name: name == "Event"
    is_date = lambda name: name == "Date"
    is_time = lambda name: name == "TIME:UTC"

    is_channel_functions = [
        is_conductivity,
        is_temperature_primary,
        is_pressure,
        is_chlorophyll,
        is_O2saturation,
        is_seapressure,
        is_depth,
        is_salinity,
        is_castdirection,
        is_event,
        is_date,
        is_time,
    ]

    # Create empty lists to add information of available channels to
    channel = []
    index = []
    unit = []
    input_format = []
    output_format = []
    na_value = []

    for col in col_list:
        for i in range(len(channel_list)):
            if is_channel_functions[i](col):
                channel.append(channel_list[i]), index.append(
                    index_list[i]
                ), unit.append(unit_list[i]), input_format.append(
                    input_format_list[i]
                ), output_format.append(
                    output_format_list[i]
                ), na_value.append(
                    na_value_list[i]
                )

    header = pd.DataFrame([channel, index, unit, input_format, output_format, na_value])
    # column_names_header = dict.fromkeys(header.columns, '')  # set empty column names
    # header = header.rename(columns=column_names_header)

    # If plan shapes not aligned error, then may need to add another variable to channel_spellings
    try:
        ctd_data_header = pd.concat((header, ctd_data))
        ctd_data_header.to_csv(dest_dir + output_name, index=False, header=False)
        # Do a check to make sure the right columns were captured or discarded
        if np.nan in ctd_data_header.columns:
            print("Warning: Not all the right columns have been captured or discarded")
            print("Output csv file data column names:", ctd_data_header.columns)
    except ValueError as e:
        print(
            "If plan shapes are not aligned, then may need to add new channel name "
            "spellings to channel_spellings list"
        )
        print(e.args)

    return


def plot_track_location(
    dest_dir: str,
    year: str,
    cruise_number: str,
    left_lon=None,
    right_lon=None,
    bot_lat=None,
    top_lat=None,
):
    """
    Read in a csv file and output a map
    Inputs:
        - dest_dir
        - year
        - cruise_number
        - left_lon, right_lon, bot_lat, top_lat: longitude and latitude range of map extent
    Outputs:
        - A map showing sampling locations
    """
    # Create a folder for figures if it doesn't already exist
    figure_dir = os.path.join(dest_dir, "FIG")
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    input_name = str(year) + "-" + str(cruise_number) + "_header-merge.csv"
    input_filename = dest_dir + input_name
    header = pd.read_csv(input_filename, header=0)
    header["lat_degree"] = header["LOC:LATITUDE"].str[:2].astype(int)
    header["lat_min"] = header["LOC:LATITUDE"].str[3:10].astype(float)
    header["lat"] = header["lat_degree"] + header["lat_min"] / 60
    header["lon_degree"] = header["LOC:LONGITUDE"].str[:3].astype(int)
    header["lon_min"] = header["LOC:LONGITUDE"].str[4:12].astype(float)
    header["lon"] = 0 - (header["lon_degree"] + header["lon_min"] / 60)
    # event = header['LOC:STATION'].astype(str)

    lon = header["lon"].tolist()
    lat = header["lat"].tolist()
    # event = event.tolist()

    # If map extent limits are not provided, then make some up based on the data
    coord_limits = [
        np.floor(np.min(header["lon"]) * 10) / 10,
        np.ceil(np.max(header["lon"]) * 10) / 10,
        np.floor(np.min(header["lat"]) * 10) / 10,
        np.ceil(np.max(header["lat"]) * 10) / 10,
    ]

    left_lon = coord_limits[0] if left_lon is None else left_lon
    right_lon = coord_limits[1] if right_lon is None else right_lon
    bot_lat = coord_limits[2] if bot_lat is None else bot_lat
    top_lat = coord_limits[3] if top_lat is None else top_lat

    m = Basemap(
        llcrnrlon=left_lon,
        llcrnrlat=bot_lat,
        urcrnrlon=right_lon,
        urcrnrlat=top_lat,
        projection="lcc",
        resolution="h",
        lat_0=0.5 * (bot_lat + top_lat),
        lon_0=0.5 * (left_lon + right_lon),
    )  # lat_0=53.4, lon_0=-129.0)

    x, y = m(lon, lat)

    fig = plt.figure(num=None, figsize=(8, 6), dpi=100)
    m.drawcoastlines(linewidth=0.2)
    m.drawmapboundary(fill_color="white")
    # m.fillcontinents(color='0.8')
    m.drawrivers()

    m.scatter(x, y, marker="D", color="m", s=5)
    # m.plot(x, y, marker='D', color='m', markersize=4)
    #   for event, xpt, ypt in zip(event, x, y):
    #       plt.text(xpt, ypt, event)

    parallels = np.arange(bot_lat, top_lat, 0.5)
    m.drawparallels(
        parallels, labels=[True, False, True, False]
    )  # draw parallel lat lines
    meridians = np.arange(left_lon, right_lon, 0.5)
    m.drawmeridians(meridians, labels=[False, False, False, True])
    plt.title(year + "-" + cruise_number)
    plt.tight_layout()
    plt.savefig(
        os.path.join(figure_dir, f"{year}-{cruise_number}_basemap_sampling_area.png")
    )  # 'Fig_1.png'
    plt.close(fig)

    # ----- create a second map w/ Cartopy just to double check - had an issue with Basemap once only

    Map = plt.axes(projection=ccrs.PlateCarree())
    Map.set_extent(
        [left_lon, right_lon, bot_lat, top_lat]
    )  # try left_lon, right_lon, bot_lat, top_lat
    x, y = (lon, lat)
    Map.coastlines()
    gl = Map.gridlines(
        crs=ccrs.PlateCarree(),
        linewidth=0.5,
        color="black",
        alpha=0.5,
        linestyle="--",
        draw_labels=True,
    )
    gl.top_labels = False
    gl.left_labels = True
    gl.bottom_labels = True
    gl.right_labels = False
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()

    gl.xlabel_style = {"color": "black", "weight": "bold", "size": 6}
    gl.ylabel_style = {"color": "black", "weight": "bold", "size": 6}

    cax = plt.scatter(x, y, transform=ccrs.PlateCarree(), marker=".", color="red", s=25)
    plt.title(year + "-" + cruise_number)
    plt.tight_layout()
    plt.savefig(
        os.path.join(figure_dir, f"{year}-{cruise_number}_cartopy_sampling_area.png")
    )  # 'Figure_1.png'
    plt.close()
    return


def PLOT_PRESSURE_DIFF(dest_dir: str, year: str, cruise_number: str, input_ext: str):
    """
    Read in a csv file and output a plot to check zero-order holds
    Inputs:
        - dest_dir
        - year
        - cruise_number
        - input_ext: '_CTD_DATA-6linehdr.csv' or '_CTD_DATA-6linehdr_corr_hold.csv'
    Outputs:
        - a plot showing the time derivative of raw pressure
    """

    # Create a folder for figures if it doesn't already exist
    figure_dir = os.path.join(dest_dir, "FIG")
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    # input_name will depend on the need for zero order holds
    input_name = str(year) + "-" + str(cruise_number) + input_ext
    input_filename = dest_dir + input_name
    ctd_data = pd.read_csv(input_filename, header=None, low_memory=False)
    ctd_data = ctd_data.rename(
        columns=ctd_data.iloc[1]
    )  # assign the second row as column names
    ctd_data = ctd_data.rename(
        columns={
            "Oxygen:Dissolved:Saturation": "Oxygen",
            "Salinity:CTD": "Salinity",
            "TIME:UTC": "TIME",
        }
    )
    ctd = ctd_data.iloc[6:]
    ctd.index = np.arange(0, len(ctd))
    # ctd = ctd[1000:4000] # to limit the number of records plotted -

    pressure = ctd["Pressure"].apply(pd.to_numeric)
    pressure_lag = pressure[1:]
    pressure_lag.index = np.arange(0, len(pressure_lag))
    pressure_diff = pressure_lag - pressure

    fig = plt.figure(num=None, figsize=(14, 6), dpi=100)
    plt.plot(pressure_diff, color="blue", linewidth=0.5, label="Pressure_diff")
    plt.ylabel("Pressure (decibar)")
    plt.xlabel("Scans")
    plt.grid()
    plt.legend()
    plt.title(year + "-" + cruise_number + " " + input_ext)
    plt.tight_layout()
    plt.savefig(os.path.join(figure_dir, "zero_order_holds_" + input_ext + ".png"))
    plt.close(fig)
    return


def CREATE_CAST_VARIABLES(
    year: str, cruise_number: str, dest_dir: str, input_ext: str
) -> tuple:
    """
    Read in a csv file and output data dictionaries to hold profile data
    Inputs:
        - dest_dir
        - year
        - cruise_number
        - input_ext: '_CTD_DATA-6linehdr.csv' or '_CTD_DATA-6linehdr_corr_hold.csv'
    Outputs:
        - three dictionaries containing casts, downcasts and upcasts
    """
    input_name = str(year) + "-" + str(cruise_number) + input_ext
    input_filename = dest_dir + input_name
    ctd_data = pd.read_csv(
        input_filename, header=None, low_memory=False
    )  # read data without header
    ctd_data = ctd_data.rename(
        columns=ctd_data.iloc[1]
    )  # assign the second row as column names
    ctd_data = ctd_data.rename(
        columns={
            "Oxygen:Dissolved:Saturation": "Oxygen",
            "Salinity:CTD": "Salinity",
            "TIME:UTC": "TIME",
        }
    )
    ctd = ctd_data.iloc[6:]
    # drop NaNs from Zero order holds correction
    # (not including O or F in case they aren't present - but will capture)
    if input_ext == "_CTD_DATA-6linehdr.csv":
        pass  # ctd = ctd
    elif input_ext == "_CTD_DATA-6linehdr_corr_hold.csv":
        ctd = ctd.dropna(
            axis=0,
            subset=[
                "Conductivity",
                "Temperature",
                "Pressure_Air",
                "Pressure",
                "Depth",
                "Salinity",
            ],
            how="all",
        )  # I don't think this does anything - NaNs now dropped in Correct_Hold stage
    ctd = ctd.copy()
    cols = ctd.columns[0:-4]
    # cols_names = ctd.columns.tolist()
    ctd[cols] = ctd[cols].apply(pd.to_numeric, errors="coerce", axis=1)
    ctd["Cast_direction"] = ctd["Cast_direction"].str.strip()

    # Fix the iteration to account for event numbers that aren't sequential or start at n>1
    unique_event_numbers = ctd["Event_number"].unique()

    var_holder = {}
    for i in unique_event_numbers:
        # Assign values of type DataFrame
        var_holder["cast" + str(i)] = ctd.loc[(ctd["Event_number"] == str(i))]
    # var_holder['Processing_history'] = ""

    # Downcast dictionary
    var_holder_d = {}
    # for i in range(1, n + 1):
    for i in unique_event_numbers:
        var_holder_d["cast" + str(i)] = ctd.loc[
            (ctd["Event_number"] == str(i)) & (ctd["Cast_direction"] == "d")
        ]
    # var_holder_d['Processing_history'] = ""

    # Upcast dictionary
    var_holder_u = {}
    # for i in range(1, n + 1, 1):
    for i in unique_event_numbers:
        var_holder_u["cast" + str(i)] = ctd.loc[
            (ctd["Event_number"] == str(i)) & (ctd["Cast_direction"] == "u")
        ]
    # var_holder_u['Processing_history'] = ""

    return var_holder, var_holder_d, var_holder_u


def format_processing_plot(
    ax: plt.Axes,
    x_var_name: str,
    x_var_units,
    y_var_name: str,
    y_var_units: str,
    plot_title: str,
    invert_yaxis: bool,
    add_legend: bool = False,
) -> None:
    """
    Format a plot that has already been initialized.
    inputs:
        - ax: from fig, ax = plt.subplots()
        - var_name: one of Temperature, Conductivity, Salinity,
                    Fluorescence, Oxygen, Oxygen_mL_L, Oxygen_umol_kg
        - var_units: the units corresponding to the selected var_name
        - plot_title: Should indicate which processing step the plots are at
        - add_legend: If True then add a legend to the plot, default False
    """
    if invert_yaxis:
        ax.invert_yaxis()
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position("top")

    # Add ticks to top and right sides of the plot, like IOS Shell does
    ax.tick_params(
        bottom=True,
        top=True,
        left=True,
        right=True,
        labelbottom=True,
        labeltop=True,
        labelleft=True,
        labelright=True,
    )

    # For Oxygen_mL_L and Oxygen_umol_kg, remove the units at the end of the var_name
    # since the units will go in brackets after
    x_var_name = x_var_name.split("_")[0]

    if x_var_units is not None:
        ax.set_xlabel(f"{x_var_name} ({x_var_units})")
    else:
        ax.set_xlabel(f"{x_var_name}")
    ax.set_ylabel(f"{y_var_name} ({y_var_units})")
    ax.set_title(plot_title, fontsize=5)
    if add_legend:
        ax.legend()
    plt.tight_layout()
    return


def do_ts_plot(
        figure_dir: str, plot_title: str, figure_filename: str, cast_d: dict, cast_u=None,
):
    """Do T-S plot with isopycnal lines
    Reference: https://github.com/larsonjl/earth_data_tools/
    """
    fig, ax = plt.subplots()

    cast_numbers = [cast_i for cast_i in cast_d.keys()]

    # Initialize the min / max values for plotting isopycnals
    t_min = cast_d[cast_numbers[0]].Temperature.min()
    t_max = cast_d[cast_numbers[0]].Temperature.max()
    s_min = cast_d[cast_numbers[0]].Salinity.min()
    s_max = cast_d[cast_numbers[0]].Salinity.max()

    # Plot each cast
    for cast_i in cast_numbers:
        ax.plot(cast_d[cast_i].Salinity, cast_d[cast_i].Temperature, color="b")

        # Update temperature and salinity ranges
        t_min = min(t_min, cast_d[cast_i].Temperature.min())
        t_max = max(t_max, cast_d[cast_i].Temperature.max())
        s_min = min(s_min, cast_d[cast_i].Salinity.min())
        s_max = max(s_max, cast_d[cast_i].Salinity.max())

        if cast_u is not None:
            ax.plot(cast_u[cast_i].Salinity, cast_u[cast_i].Temperature, color="b")

            # Update temperature and salinity ranges
            t_min = min(t_min, cast_u[cast_i].Temperature.min())
            t_max = max(t_max, cast_u[cast_i].Temperature.max())
            s_min = min(s_min, cast_u[cast_i].Salinity.min())
            s_max = max(s_max, cast_u[cast_i].Salinity.max())

    # Finalize the min / max values for plotting isopycnals
    t_min -= 1
    t_max += 1
    s_min -= 1
    s_max += 1

    # Calculate how many gridcells we need in the x and y dimensions
    xdim = int(np.ceil(s_max - s_min) / 0.1)
    ydim = int(np.ceil(t_max - t_min))
    dens = np.zeros((ydim, xdim))

    # Create temp and salt vectors of appropiate dimensions
    ti = np.linspace(0, ydim, ydim) + t_min
    si = np.linspace(1, xdim, xdim) * 0.1 + s_min

    # Loop to fill in grid with densities
    for j in range(0, int(ydim)):
        for i in range(0, int(xdim)):
            dens[j, i] = gsw.rho(si[i], ti[j], 0)

    # Subtract 1000 to convert to sigma-t
    sigmat = dens - 1000

    # Add the isopycnal contour lines to the plot
    CS = ax.contour(si, ti, sigmat, linestyles='dashed', colors='k')
    plt.clabel(CS, fontsize=12, inline=1, fmt='%.2f')  # Label every second level

    # Do final plot formatting
    format_processing_plot(
        ax,
        x_var_name="Salinity",
        x_var_units="PSS-78",
        y_var_name="Temperature",
        y_var_units="C",
        plot_title=plot_title,
        invert_yaxis=False,
    )
    plt.savefig(os.path.join(figure_dir, figure_filename))
    plt.close(fig)
    return


def plot_by_other(
        year: str, cruise_number: str, dest_dir: str, input_ext: str
) -> None:
    """
    DEPRECATED
    Originally part of first_plots() function
    Plot by group, by index, and by profile
    """
    # Input ext specifies whether to use 6lineheader file with or without
    # the zero-order hold removed
    cast, cast_d, cast_u = CREATE_CAST_VARIABLES(
        year, cruise_number, dest_dir, input_ext
    )

    # get variables; cast numbering might not start from 1!
    for cast_i in cast.keys():
        vars_available = list(dict.fromkeys(cast[cast_i]))
        break
    # vars_available = list(dict.fromkeys(cast['cast1']))

    # -----------------------  Plot profiles by group --------------------------
    plot_by_group = False
    if plot_by_group:
        # T, C, O, S, F for each profile in one plot
        fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, sharey=True)
        # Temperature
        ax1.plot(
            cast_d["cast1"].Temperature,
            cast_d["cast1"].Pressure,
            color="red",
            label="cast_down",
        )
        ax1.plot(
            cast_u["cast1"].Temperature,
            cast_u["cast1"].Pressure,
            "--",
            color="red",
            label="cast_up",
        )
        ax1.set_ylabel("Pressure(decibar)", fontsize=8)
        ax1.set_ylim(ax1.get_ylim()[::-1])
        ax1.set_xlabel("Temperature(C)", fontsize=8)
        ax1.xaxis.set_label_position("top")
        ax1.xaxis.set_ticks_position("top")
        ax1.set_title("Pre-Processing", fontsize=5)
        ax1.legend()

        # Conductivity
        ax2.plot(
            cast_d["cast1"].Conductivity,
            cast_d["cast1"].Pressure,
            color="goldenrod",
            label="cast_down",
        )
        ax2.plot(
            cast_u["cast1"].Conductivity,
            cast_u["cast1"].Pressure,
            "--",
            color="goldenrod",
            label="cast1_up",
        )
        ax2.set_ylabel("Pressure(decibar)", fontsize=8)
        ax2.set_ylim(ax1.get_ylim()[::-1])
        ax2.set_xlabel("Conductivity (S/cm)", fontsize=8)
        ax2.xaxis.set_label_position("top")
        ax2.xaxis.set_ticks_position("top")
        ax2.set_title("Pre-Processing", fontsize=5)
        ax2.legend()

        # Oxygen
        for var in vars_available:
            if var == "Oxygen":
                ax3.plot(
                    cast_d["cast1"].Oxygen,
                    cast_d["cast1"].Pressure,
                    color="black",
                    label="cast_down",
                )
                ax3.plot(
                    cast_u["cast1"].Oxygen,
                    cast_u["cast1"].Pressure,
                    "--",
                    color="black",
                    label="cast_up",
                )
                ax3.set_ylabel("Pressure(decibar)", fontsize=8)
                ax3.set_ylim(ax1.get_ylim()[::-1])
                ax3.set_xlabel("Oxygen Saturation (%)", fontsize=8)
                ax3.xaxis.set_label_position("top")
                ax3.xaxis.set_ticks_position("top")
                ax3.set_title("Pre-Processing", fontsize=5)
                ax3.legend()
            elif var == "Fluorescence":
                ax5.plot(
                    cast_d["cast1"].Fluorescence,
                    cast_d["cast1"].Pressure,
                    color="green",
                    label="cast_down",
                )
                ax5.plot(
                    cast_u["cast1"].Fluorescence,
                    cast_u["cast1"].Pressure,
                    "--",
                    color="green",
                    label="cast1_up",
                )
                ax5.set_ylabel("Pressure(decibar)", fontsize=8)
                ax5.set_ylim(ax1.get_ylim()[::-1])
                ax5.set_xlabel("Fluoresence(ug/L)", fontsize=8)
                ax5.xaxis.set_label_position("top")
                ax5.xaxis.set_ticks_position("top")
                ax5.set_title("Pre-Processing", fontsize=5)
                ax5.legend()

        # Salinity
        ax4.plot(
            cast_d["cast1"].Salinity,
            cast_d["cast1"].Pressure,
            color="blue",
            label="cast_down",
        )
        ax4.plot(
            cast_u["cast1"].Salinity,
            cast_u["cast1"].Pressure,
            "--",
            color="blue",
            label="cast_up",
        )
        ax4.set_ylabel("Pressure(decibar)", fontsize=8)
        ax4.set_ylim(ax1.get_ylim()[::-1])
        ax4.set_xlabel("Salinity", fontsize=8)
        ax4.xaxis.set_label_position("top")
        ax4.xaxis.set_ticks_position("top")
        ax4.set_title("Pre-Processing", fontsize=5)
        ax4.legend()

    # -------------------------   Plot by Index and by Profile-------------------------
    # separate plot for T, C, O, S, F of each profile
    plot_by_index_and_profile = False
    if plot_by_index_and_profile:
        # Temperature
        fig, ax = plt.subplots()
        ax.plot(
            cast_d["cast1"].Temperature,
            cast_d["cast1"].Pressure,
            color="red",
            label="cast_down",
        )
        ax.plot(
            cast_u["cast1"].Temperature,
            cast_u["cast1"].Pressure,
            "--",
            color="red",
            label="cast_up",
        )
        ax.invert_yaxis()
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        ax.set_xlabel("Temperature(C)")
        ax.set_ylabel("Pressure (decibar)")
        ax.set_title("Pre-Processing", fontsize=5)
        plt.tight_layout()
        plt.savefig(dest_dir + "Pre_Cast1_T")
        ax.legend()

        # Salinity
        fig, ax = plt.subplots()
        ax.plot(
            cast_d["cast1"].Salinity,
            cast_d["cast1"].Pressure,
            color="blue",
            label="cast1_d",
        )
        ax.plot(
            cast_u["cast1"].Salinity,
            cast_u["cast1"].Pressure,
            "--",
            color="blue",
            label="cast1_u",
        )
        ax.invert_yaxis()
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        ax.set_xlabel("Salinity")
        ax.set_ylabel("Pressure (decibar)")
        ax.set_title("Pre-Processing", fontsize=5)
        ax.legend()
        plt.tight_layout()
        plt.savefig(dest_dir + "Pre_Cast1_S")

        # Conductivity
        fig, ax = plt.subplots()
        ax.plot(
            cast_d["cast1"].Conductivity,
            cast_d["cast1"].Pressure,
            color="yellow",
            label="cast1_d",
        )
        ax.plot(
            cast_u["cast1"].Conductivity,
            cast_u["cast1"].Pressure,
            "--",
            color="yellow",
            label="cast1_u",
        )
        ax.invert_yaxis()
        ax.xaxis.set_label_position("top")
        ax.xaxis.set_ticks_position("top")
        ax.set_xlabel("Conductivity (S/cm)")
        ax.set_ylabel("Pressure (decibar)")
        ax.set_title("Pre-Processing", fontsize=5)
        ax.legend()
        plt.tight_layout()
        plt.savefig(dest_dir + "Pre_Cast1_C")

        # Oxygen
        for var in vars_available:
            if var == "Oxygen":
                fig, ax = plt.subplots()
                ax.plot(
                    cast_d["cast1"].Oxygen,
                    cast_d["cast1"].Pressure,
                    color="black",
                    label="cast1_d",
                )
                ax.plot(
                    cast_u["cast1"].Oxygen,
                    cast_u["cast1"].Pressure,
                    "--",
                    color="black",
                    label="cast1_u",
                )
                ax.invert_yaxis()
                ax.xaxis.set_label_position("top")
                ax.xaxis.set_ticks_position("top")
                ax.set_xlabel("Oxygen Saturation (%)")  # Check unit here
                ax.set_ylabel("Pressure (decibar)")
                ax.set_title("Pre-Processing", fontsize=5)
                ax.legend()
                plt.tight_layout()
                plt.savefig(dest_dir + "Pre_Cast1_O")

            elif var == "Fluorescence":
                fig, ax = plt.subplots()
                ax.plot(
                    cast_d["cast1"].Fluorescence,
                    cast_d["cast1"].Pressure,
                    color="green",
                    label="cast1_d",
                )
                ax.plot(
                    cast_u["cast1"].Fluorescence,
                    cast_u["cast1"].Pressure,
                    "--",
                    color="green",
                    label="cast1_u",
                )
                ax.invert_yaxis()
                ax.xaxis.set_label_position("top")
                ax.xaxis.set_ticks_position("top")
                ax.set_xlabel("Fluorescence (ug/L)")  # Check unit here
                ax.set_ylabel("Pressure (decibar)")
                ax.set_title("Pre-Processing", fontsize=5)
                ax.legend()
                plt.tight_layout()
                plt.savefig(dest_dir + "Pre_Cast1_F")

        fig, ax = plt.subplots()
        ax.plot(
            cast_d["cast1"].Salinity,
            cast_d["cast1"].Temperature,
            color="red",
            label="cast1_d",
        )
        ax.plot(
            cast_u["cast1"].Salinity,
            cast_u["cast1"].Temperature,
            "--",
            color="blue",
            label="cast1_u",
        )
        ax.set_xlabel("Salinity")
        ax.set_ylabel("Temperature (C)")
        ax.set_title("Pre-Processing T-S Plot")
        ax.legend()
        plt.tight_layout()
        plt.savefig(dest_dir + "Pre_Cast1_T-S.png")

        number_of_colors = len(cast)
        color = [
            "#" + "".join([random.choice("0123456789ABCDEF") for j in range(6)])
            for i in range(number_of_colors)
        ]

    return


def first_plots(year: str, cruise_number: str, dest_dir: str, input_ext: str) -> None:
    """
    Plot pre-processing and after Zero-order Holds if needed
    inputs:
        - year
        - cruise_number
        - dest_dir
        - input_ext: '_CTD_DATA-6linehdr.csv' or '_CTD_DATA-6linehdr_corr_hold.csv'
    outputs:
        - profile plots of all available variables plus a T-S plot are output to
        dest_dir/FIG/, but the function returns nothing
    """

    # Create a folder for figures if it doesn't already exist
    figure_dir = os.path.join(dest_dir, "FIG")
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    # Input ext specifies whether to use 6lineheader file with or without
    # the zero-order hold removed
    cast, cast_d, cast_u = CREATE_CAST_VARIABLES(
        year, cruise_number, dest_dir, input_ext
    )

    # number_of_colors = len(cast)
    # color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    #          for i in range(number_of_colors)]

    # get variables; cast numbering might not start from 1!
    for cast_i in cast.keys():
        vars_available = list(dict.fromkeys(cast[cast_i]))
        break
    # vars_available = list(dict.fromkeys(cast['cast1']))

    # Iterate through all the channels, plot data from all casts on one plot per channel
    for j, var in enumerate(VARIABLES_POSSIBLE):
        if var in vars_available:
            fig, ax = plt.subplots()
            for cast_i in cast.keys():
                ax.plot(
                    cast_d[cast_i].loc[:, var],
                    cast_d[cast_i].Pressure,
                    color=VARIABLE_COLOURS[j],
                )
                # label='cast' + str(i + 1))
                ax.plot(
                    cast_u[cast_i].loc[:, var],
                    cast_u[cast_i].Pressure,
                    color=VARIABLE_COLOURS[j],
                )  # '--'
                # label='cast' + str(i + 1))
                # ax.plot(cast_d['cast1'].Salinity, cast_d['cast1'].Pressure, color='blue', label='cast1')
                # ax.plot(cast_u['cast1'].Salinity, cast_u['cast1'].Pressure, '--', color='blue', label='cast1')
            format_processing_plot(
                ax,
                x_var_name=var,
                x_var_units=VARIABLE_UNITS[j],
                y_var_name="Pressure",
                y_var_units="dbar",
                plot_title="Pre-Processing",
                invert_yaxis=True,
            )
            plt.savefig(os.path.join(figure_dir, f"Pre_Processing_{var[0]}.png"))
            plt.close(fig)

    # TS Plot add labeled isopycnals to T-S plots as in IOS Shell
    do_ts_plot(
        figure_dir, "Pre-Processing T-S Plot", "Pre_Processing_T-S.png", cast_d, cast_u
    )
    plt.savefig(os.path.join(figure_dir, "Pre_Processing_T-S.png"))
    plt.close(fig)

    # pressure check
    fig, ax = plt.subplots()
    for cast_i in cast.keys():
        ax.plot(
            cast_d[cast_i].Conductivity[0:20],
            cast_d[cast_i].Pressure[0:20],
            color="goldenrod",
        )  # , label='cast' + str(i + 1))
        ax.plot(
            cast_u[cast_i].Conductivity[-20:-1],
            cast_u[cast_i].Pressure[-20:-1],
            color="goldenrod",
        )  # label='cast' + str(i + 1))
    # ax.plot(cast_d['cast1'].Conductivity[0:10], cast_d['cast1'].Pressure[0:10],
    #         color='blue', label='cast1')
    # ax.plot(cast_d['cast2'].Conductivity[0:10], cast_d['cast2'].Pressure[0:10],
    #         color='red', label='cast2')
    format_processing_plot(
        ax,
        x_var_name="Conductivity",
        x_var_units="mS/cm",
        y_var_name="Pressure",
        y_var_units="dbar",
        plot_title="Checking need for Pressure correction",
        invert_yaxis=True,
    )
    plt.savefig(os.path.join(figure_dir, "PressureCorrection_need_CvP.png"))
    plt.close(fig)

    fig, ax = plt.subplots()
    for cast_i in cast.keys():
        ax.plot(
            cast_d[cast_i].Conductivity[0:20],
            cast_d[cast_i].Depth[0:20],
            color="goldenrod",
        )
        # label='cast' + str(i + 1))
        ax.plot(
            cast_u[cast_i].Conductivity[-20:-1],
            cast_u[cast_i].Depth[-20:-1],
            color="goldenrod",
        )  # label='cast' + str(i + 1))
    # ax.plot(cast_d['cast1'].Conductivity[0:10], cast_d['cast1'].Depth[0:10],
    #         color='blue', label='cast1')
    # ax.plot(cast_d['cast2'].Conductivity[0:10], cast_d['cast2'].Depth[0:10],
    #         color='red', label='cast2')
    format_processing_plot(
        ax,
        x_var_name="Conductivity",
        x_var_units="mS/cm",
        y_var_name="Depth",
        y_var_units="m",
        plot_title="Checking need for Pressure correction",
        invert_yaxis=True,
    )
    plt.savefig(os.path.join(figure_dir, "PressureCorrection_need_CvD.png"))
    plt.close(fig)
    return


def check_for_zoh(
    dest_dir, year: str, cruise_number: str, sampling_interval: float
) -> bool:
    """
    Compute first order differences on pressure data to determine whether
    a correction for zero-order holds is needed.
    From DFO Technical Report 314:
    'The analog-to-digital (A2D) converter on RBR instruments must recalibrate once per
    minute.'
    inputs
        - dest_dir: destination/working directory
        - year
        - cruise_number
        - sampling_interval: amount of time in seconds between records
    outputs
        - boolean flag indicating whether zero-order holds are present in the data
    """

    input_name = str(year) + "-" + str(cruise_number) + "_CTD_DATA-6linehdr.csv"
    input_filename = dest_dir + input_name
    ctd_data = pd.read_csv(input_filename, header=None, low_memory=False)
    # assign the second row as column names
    ctd_data = ctd_data.rename(columns=ctd_data.iloc[1])
    ctd_data = ctd_data.rename(
        columns={
            "Oxygen:Dissolved:Saturation": "Oxygen",
            "Salinity:CTD": "Salinity",
            "TIME:UTC": "TIME",
        }
    )
    ctd = ctd_data.iloc[6:]
    ctd.index = np.arange(0, len(ctd))

    pressure = ctd["Pressure"].apply(pd.to_numeric)
    pressure_diffs = np.diff(pressure)

    print("Number of pressure records:", len(pressure))
    print("Sum of zero pressure differences:", sum(pressure_diffs == 0))
    print(
        "Intervals between zero pressure differences:",
        np.diff(np.where(pressure_diffs == 0)[0]),
        sep="\n",
    )

    sec2min = 1 / 60  # Convert seconds to minutes b/c sampling interval in seconds
    if sum(pressure_diffs == 0) >= np.floor(
        len(pressure) * sampling_interval * sec2min
    ):
        zoh_correction_needed = True
    else:
        zoh_correction_needed = False
    return zoh_correction_needed


def CORRECT_HOLD(
    dest_dir: str, year: str, cruise_number: str, metadata_dict: dict
) -> None:
    """
    Read 6linehdr.csv and correct for zero order holds.  Look for repeat values in
    Pressure and replace with NaN, then
    look for repeats in the other sensors at the same place and replace
    with NaN..

    Adapted from RSKtools function RSKcorrecthold: 'This function identifies zero-hold
    points by looking
    for where consecutive differences for each channel are equal to zero, and replaces
    them with Nan or an
    interpolated value."  This function uses Nan. SH

    Output a new csv with the corrected values.
    """
    input_name = str(year) + "-" + str(cruise_number) + "_CTD_DATA-6linehdr.csv"
    output_name = (
        str(year) + "-" + str(cruise_number) + "_CTD_DATA-6linehdr_corr_hold.csv"
    )
    input_filename = dest_dir + input_name
    ctd_data = pd.read_csv(input_filename, header=None, low_memory=False)
    ctd_data = ctd_data.rename(
        columns=ctd_data.iloc[1]
    )  # assign the second row as column names
    ctd_data = ctd_data.rename(
        columns={
            "Oxygen:Dissolved:Saturation": "Oxygen",
            "Salinity:CTD": "Salinity",
            "TIME:UTC": "TIME",
        }
    )

    have_fluor = True if "Fluorescence" in ctd_data.columns else False
    have_oxy = True if "Oxygen" in ctd_data.columns else False

    header = ctd_data.iloc[0:6]  # keep the header to use at the end
    ctd = ctd_data.iloc[6:]
    ctd = ctd.copy()
    # vars_available = list(dict.fromkeys(ctd))
    # print(vars_available)
    # cols = ctd.columns[0:-4]
    # ctd[cols] = ctd[cols].apply(pd.to_numeric, errors='coerce', axis=1)
    ctd.index = np.arange(0, len(ctd))
    pressure = ctd["Pressure"].apply(pd.to_numeric)
    pressure_lag = pressure[1:]
    pressure_lag.index = np.arange(0, len(pressure_lag))
    air = ctd["Pressure_Air"].apply(pd.to_numeric)
    air_lag = air[1:]
    air_lag.index = np.arange(0, len(air_lag))

    conductivity = ctd["Conductivity"].apply(pd.to_numeric)
    conductivity_lag = conductivity[1:]
    conductivity_lag.index = np.arange(0, len(conductivity_lag))

    temperature = ctd["Temperature"].apply(pd.to_numeric)
    temperature_lag = temperature[1:]
    temperature_lag.index = np.arange(0, len(temperature_lag))

    if have_fluor:
        fluorescence = ctd["Fluorescence"].apply(pd.to_numeric)
        fluorescence_lag = fluorescence[1:]
        fluorescence_lag.index = np.arange(0, len(fluorescence_lag))
    if have_oxy:
        oxygen = ctd["Oxygen"].apply(pd.to_numeric)
        oxygen_lag = oxygen[1:]
        oxygen_lag.index = np.arange(0, len(oxygen_lag))

    depth = ctd["Depth"].apply(pd.to_numeric)
    depth_lag = depth[1:]
    depth_lag.index = np.arange(0, len(depth_lag))

    salinity = ctd["Salinity"].apply(pd.to_numeric)
    salinity_lag = salinity[1:]
    salinity_lag.index = np.arange(0, len(salinity_lag))

    for i in range(len(ctd) - 1):
        if pressure[i] == pressure_lag[i]:
            # pressure.iloc[i + 1] = np.nan
            if conductivity[i] == conductivity_lag[i]:
                conductivity.iloc[i + 1] = np.nan
            if air[i] == air_lag[i]:
                air.iloc[i + 1] = np.nan
            if temperature[i] == temperature_lag[i]:
                temperature.iloc[i + 1] = np.nan
            if have_fluor:
                if fluorescence[i] == fluorescence_lag[i]:
                    fluorescence.iloc[i + 1] = np.nan
            if have_oxy:
                if oxygen[i] == oxygen_lag[i]:
                    oxygen.iloc[i + 1] = np.nan
            # if depth[i] == depth_lag[i]:
            # depth.iloc[i + 1] = np.nan
            if salinity[i] == salinity_lag[i]:
                salinity.iloc[i + 1] = np.nan

    # ctd['Pressure'] = pressure  # this worked when pressure was set to NaN
    ctd["Conductivity"] = conductivity
    ctd["Temperature"] = temperature
    ctd["Pressure_Air"] = air
    if have_fluor:
        ctd["Fluorescence"] = fluorescence
    if have_oxy:
        ctd["Oxygen"] = oxygen
    # ctd['Depth'] = depth
    ctd["Salinity"] = salinity

    # drop the NaNs before they get into the CSV for IOS Shell Processing
    ctd.dropna(
        axis=0,
        subset=["Conductivity", "Temperature", "Pressure_Air", "Salinity"],
        how="all",
        inplace=True,
    )  # 'Pressure', 'Depth',  used to be in this list  #sometimes 'any' is required
    # ctd = ctd.reset_index(drop=True)

    metadata_dict["Processing_history"] = (
        "-Zero-Order Holds Correction:|"
        " Correction type = Substitute with Nan"
        " Corrections applied:|"
        " All channels corrected where zero-order holds concur with Pressure Holds:|"
    )
    metadata_dict["ZEROORDER_Time"] = datetime.now()

    # SH: Improve metadata_dict entry here? Ask Lu

    # ctd = ctd.rename(columns={'Oxygen': 'Oxygen:Dissolved:Saturation', 'Salinity': 'Salinity:CTD',
    #                           'TIME': 'TIME:UTC'})
    column_names = dict.fromkeys(ctd.columns, "")  # set empty column names
    ctd = ctd.rename(columns=column_names)  # remove column names

    column_names_header = dict.fromkeys(header.columns, "")  # set empty column names
    header = header.rename(columns=column_names_header)

    ctd_header = pd.concat((header, ctd))  # header.append(ctd) deprecated
    ctd_header.to_csv(dest_dir + output_name, index=False, header=False)
    return


def CALIB(
    var: dict,
    var_downcast: dict,
    var_upcast: dict,
    metadata_dict: dict,
    zoh: bool,
    pd_correction_value=0,
) -> tuple:
    """
    Correct pressure and depth data
    Inputs:
        - cast, downcast, upcast and metadata dictionaries
        - correction_value: for pressure and depth data
        - metadata_dict
        - zoh: if "no", then zero-order hold correction was not applied, else "yes"
    Outputs:
        - cast, downcast, upcast and metadata dictionaries after pressure correction
    """
    var1 = deepcopy(var)
    var2 = deepcopy(var_downcast)
    var3 = deepcopy(var_upcast)

    for cast_i in var1.keys():
        var1[cast_i].Pressure = var1[cast_i].Pressure + pd_correction_value
        var1[cast_i].Depth = var1[cast_i].Depth + pd_correction_value
        var2[cast_i].Pressure = var2[cast_i].Pressure + pd_correction_value
        var2[cast_i].Depth = var2[cast_i].Depth + pd_correction_value
        var3[cast_i].Pressure = var3[cast_i].Pressure + pd_correction_value
        var3[cast_i].Depth = var3[cast_i].Depth + pd_correction_value

    # check if a correction was done - need to see if this is the first addition of "processing history'
    if not zoh:
        metadata_dict["Processing_history"] = ""

    metadata_dict["Processing_history"] += (
        "-CALIB parameters:|"
        " Calibration type = Correct|"
        " Calibrations applied:|"
        f" Pressure (decibar) = {pd_correction_value}|"
        f" Depth (meters) = {pd_correction_value}|"
    )

    metadata_dict["CALIB_Time"] = datetime.now()

    return var1, var2, var3


# Data despiking has not been implemented before clipping.


def CLIP_CAST(
    var: dict, metadata_dict: dict, limit_pressure_change: float, cast_direction: str
):
    """
    CLIP the unstable measurement from sea surface and bottom on the downcast OR upcast
     Inputs:
         - Upcast, metadata dictionary,
         - limit_pressure_change: limit drop for downcast, limit rise for upcast
     Outputs:
         - Upcast after removing records near surface and bottom
    direction: 'down' or 'up'
    """
    var_clip = deepcopy(var)
    for cast_i in var.keys():
        pressure = var_clip[cast_i].Pressure
        diff = var_clip[cast_i].Pressure.diff()
        index_start = pressure.index[0]
        # index_end = pressure.index[-1]
        if cast_direction == "down":
            limit_drop = limit_pressure_change
            diff_mask = diff > limit_drop
        elif cast_direction == "up":
            limit_rise = limit_pressure_change
            diff_mask = diff < limit_rise
        else:
            print(f"cast_direction {cast_direction} is invalid. Ending program")
            return
        diff_rise = diff.loc[diff_mask]
        for j in range(len(diff.loc[diff_mask])):
            index_1 = diff_rise.index[j]
            if (
                (diff_rise.index[j + 1] == index_1 + 1)
                and (diff_rise.index[j + 2] == index_1 + 2)
                and (diff_rise.index[j + 3] == index_1 + 3)
                and (diff_rise.index[j + 4] == index_1 + 4)
                and (diff_rise.index[j + 5] == index_1 + 5)
                and (diff_rise.index[j + 6] == index_1 + 6)
                and (diff_rise.index[j + 7] == index_1 + 7)
                and (diff_rise.index[j + 8] == index_1 + 8)
            ):
                index_end_1 = index_1 - 1
                break
        cut_start = index_end_1 - index_start

        for j in range(-1, -len(diff.loc[diff_mask]), -1):
            index_2 = diff_rise.index[j]
            if (
                (diff_rise.index[j - 1] == index_2 - 1)
                and (diff_rise.index[j - 2] == index_2 - 2)
                and (diff_rise.index[j - 3] == index_2 - 3)
                and (diff_rise.index[j - 4] == index_2 - 4)
                and (diff_rise.index[j - 5] == index_2 - 5)
            ):
                index_end_2 = index_2 + 1
                break
        cut_end = index_end_2 - index_start
        var_clip[cast_i] = var_clip[cast_i][cut_start:cut_end]

        metadata_dict["Processing_history"] += (
            "-CLIP_{}{}".format(cast_direction, cast_i)
            + ": First Record = {}".format(str(cut_start))
            + ", Last Record = {}".format(str(cut_end))
            + "|"
        )
        metadata_dict[
            "CLIP_{}_Time{}".format(cast_direction[0].upper(), cast_i.split("cast")[-1])
        ] = datetime.now()
    return var_clip


def plot_clip(cast_d_clip: dict, cast_d_pc: dict, dest_dir: str) -> None:
    """
    plot the clipped casts to check that the clipping worked
    inputs:
        - cast_d_clip: downcast data dictionary from after clip step
        - cast_d_pc: downcast data dictionary from before clip step
        - dest_dir: working/destination directory
    outputs:
        - plots of downcast pressure vs time for before and after the clip
        step are saved to dest_dir/FIG/, but nothing is returned by the function
    """

    # Create a folder for figures if it doesn't already exist
    figure_dir = os.path.join(dest_dir, "FIG")
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    fig, ax = plt.subplots()  # Before
    for cast_i in cast_d_pc.keys():
        ax.plot(cast_d_pc[cast_i].TIME, cast_d_pc[cast_i].Pressure, color="blue")
        # ax.plot(cast_u_pc[cast_i].TIME, cast_u_pc[cast_i].Pressure,
        #         color='blue')
    xticks = ax.get_xticks()
    ax.set_xticks(ticks=xticks[:: int(len(xticks) / 4)])  # Make ticks farther apart
    format_processing_plot(
        ax,
        x_var_name="Time",
        x_var_units=None,
        y_var_name="Pressure",
        y_var_units="dbar",
        plot_title="Downcasts before clip",
        invert_yaxis=True,
    )
    plt.savefig(os.path.join(figure_dir, "Before_Clip_P_vs_t.png"))
    plt.close(fig)

    fig, ax = plt.subplots()  # After
    for cast_i in cast_d_pc.keys():
        ax.plot(cast_d_clip[cast_i].TIME, cast_d_clip[cast_i].Pressure, color="blue")
        # ax.plot(cast_u_clip[cast_i].TIME, cast_u_clip[cast_i].Pressure,
        #         color='blue')
    xticks = ax.get_xticks()
    ax.set_xticks(ticks=xticks[:: int(len(xticks) / 4)])
    format_processing_plot(
        ax,
        x_var_name="Time",
        x_var_units=None,
        y_var_name="Pressure",
        y_var_units="dbar",
        plot_title="Downcasts after clip",
        invert_yaxis=True,
    )
    plt.savefig(os.path.join(figure_dir, "After_Clip_P_vs_t.png"))
    plt.close(fig)

    # # plot all cast together
    #
    # number_of_colors = len(cast)
    # color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    #          for i in range(number_of_colors)]
    #
    # fig, ax = plt.subplots()
    # for i in range(0, len(cast), 1):
    #     ax.plot(cast_d_clip['cast' + str(i + 1)].TIME,
    #             cast_d_clip['cast' + str(i + 1)].Pressure, color=color[i],
    #             label='cast' + str(i + 1))
    #     # ax.plot(cast_u['cast' + str(i+1)].Salinity, cast_u['cast' + str(i+1)].Pressure, '--', color=color[i],
    #     #         label= 'cast' + str(i+1))
    # # ax.plot(cast_d['cast1'].Salinity, cast_d['cast1'].Pressure, color='blue', label='cast1')
    # # ax.plot(cast_u['cast1'].Salinity, cast_u['cast1'].Pressure, '--', color='blue', label='cast1')
    # ax.invert_yaxis()
    # ax.xaxis.set_label_position('top')
    # ax.xaxis.set_ticks_position('top')
    # ax.set_xlabel('Time')
    # ax.set_ylabel('Pressure (decibar)')
    # ax.set_title('After Clip')
    # ax.legend()
    return


# apply a moving average FIR filter (a simple low pass )

# def filter(x, n):# n -  filter size, 9 suggest by RBR manual, choose the smallest one which can do the job
#    b = (np.ones(n))/n #numerator co-effs of filter transfer function
#    #b = repeat(1.0/n, n)
#    a = np.ones(1)  #denominator co-effs of filter transfer function
#    #y = signal.convolve(x,b) #filter output using convolution
#    #y = signal.lfilter(b, a, x) #filter output using lfilter function
#    y = signal.filtfilt(b, a, x)  # Apply a digital filter forward and backward to a signal.
#    return y


def FILTER(
    var_downcast: dict,
    var_upcast: dict,
    metadata_dict: dict,
    have_fluor: bool,
    window_width=6,
    sample_rate: int = 8,
    time_constant: float = 1 / 8,
    filter_type: int = 1,
):
    """
    Filter the temperature, conductivity, pressure, and fluorescence (if available)
    data using a low pass filter: moving average
    Inputs:
        - downcast and upcast data dictionaries
        - metadata_dict
        - have_fluor: boolean flag indicating whether fluorescence data are present;
           if they are then they are also filtered
        - window_width:
        - sample_rate:
        - time_constant:
        - filter_type: 0 or 1, corresponding to FIR or Moving Average filter
    Outputs:
        - two dictionaries containing downcast and upcast profiles after applying filter
    """

    # cast_number = len(var_downcast.keys())
    if filter_type == 0:
        Wn = (1.0 / time_constant) / (sample_rate * 2)
        # Numerator (b) and denominator (a) polynomials of the IIR filter
        b, a = signal.butter(2, Wn, "low")
        filter_name = "FIR"
    elif filter_type == 1:
        b = (
            np.ones(window_width)
        ) / window_width  # numerator co-effs of filter transfer function
        a = np.ones(1)  # denominator co-effs of filter transfer function
        filter_name = "Moving average filter"
    else:
        print("Invalid filter type:", filter_type)
        return

    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)

    # Filter select variables in each cast
    for cast_i in var1.keys():
        var1[cast_i].Temperature = signal.filtfilt(b, a, var1[cast_i].Temperature)
        var1[cast_i].Conductivity = signal.filtfilt(b, a, var1[cast_i].Conductivity)
        var1[cast_i].Pressure = signal.filtfilt(b, a, var1[cast_i].Pressure)
        var2[cast_i].Temperature = signal.filtfilt(b, a, var2[cast_i].Temperature)
        var2[cast_i].Conductivity = signal.filtfilt(b, a, var2[cast_i].Conductivity)
        var2[cast_i].Pressure = signal.filtfilt(b, a, var2[cast_i].Pressure)
        if have_fluor:
            var1[cast_i].Fluorescence = signal.filtfilt(b, a, var1[cast_i].Fluorescence)
            var2[cast_i].Fluorescence = signal.filtfilt(b, a, var2[cast_i].Fluorescence)

    metadata_dict["Processing_history"] += (
        "-FILTER parameters:|"
        f" {filter_name} was used.|"
        f" Filter width = {window_width}.|"
        "The following channel(s) were filtered.|"
        " Pressure|"
        " Temperature|"
        " Conductivity|"
    )
    if have_fluor:
        metadata_dict["Processing_history"] += " Fluorescence|"

    metadata_dict["FILTER_Time"] = datetime.now()

    return var1, var2


def plot_filter(
    cast_d_filtered: dict,
    cast_u_filtered: dict,
    cast_d_clip: dict,
    cast_u_clip: dict,
    dest_dir: str,
    have_fluor: bool,
) -> None:
    """
    Make plots to show channel values before and after filtering
    inputs:
        - cast_d_filtered, cast_u_filtered: downcast and upcast data dictionaries after filtering
        - cast_d_clip, cast_u_clip: downcast and upcast data dictionaries after clipping
        - dest_dir
        - have_fluor
    outputs:
        - plots of filtered channel vs pressure for before and after the filter
        step are saved to dest_dir/FIG/, but nothing is returned by the function
    """

    # Create a folder for figures if it doesn't already exist
    figure_dir = os.path.join(dest_dir, "FIG")
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    # vars_available = list(dict.fromkeys(cast_d_filtered['cast1']))
    # n_casts = len(cast_d_filtered)

    filtered_vars = ["Temperature", "Conductivity"]
    units = [VARIABLE_UNITS[1], VARIABLE_UNITS[2]]
    colours = [VARIABLE_COLOURS[1], VARIABLE_COLOURS[2]]
    if have_fluor:
        filtered_vars.append("Fluorescence")
        units.append(VARIABLE_UNITS[4])
        colours.append(VARIABLE_COLOURS[4])

    for j, var in enumerate(filtered_vars):  # Before
        fig, ax = plt.subplots()
        for cast_i in cast_d_clip.keys():
            ax.plot(
                cast_d_clip[cast_i].loc[:, var],
                cast_d_clip[cast_i].loc[:, "Pressure"],
                color=colours[j],
            )
            ax.plot(
                cast_u_clip[cast_i].loc[:, var],
                cast_u_clip[cast_i].loc[:, "Pressure"],
                color=colours[j],
            )
        format_processing_plot(
            ax,
            x_var_name=var,
            x_var_units=units[j],
            y_var_name="Pressure",
            y_var_units="dbar",
            plot_title="Pre-Filter",
            invert_yaxis=True,
        )
        plt.savefig(os.path.join(figure_dir, f"Pre_Filter_{var[0]}.png"))
        plt.close(fig)

    for j, var in enumerate(filtered_vars):  # After
        fig, ax = plt.subplots()
        for cast_i in cast_d_clip.keys():
            ax.plot(
                cast_d_filtered[cast_i].loc[:, var],
                cast_d_filtered[cast_i].loc[:, "Pressure"],
                color=colours[j],
            )
            ax.plot(
                cast_u_filtered[cast_i].loc[:, var],
                cast_u_filtered[cast_i].loc[:, "Pressure"],
                color=colours[j],
            )
        format_processing_plot(
            ax,
            x_var_name=var,
            x_var_units=units[j],
            y_var_name="Pressure",
            y_var_units="dbar",
            plot_title="Post-Filter",
            invert_yaxis=True,
        )
        plt.savefig(os.path.join(figure_dir, f"Post_Filter_{var[0]}.png"))
        plt.close(fig)

    return


def SHIFT_CONDUCTIVITY(
    var_downcast: dict, var_upcast: dict, metadata_dict: dict, shifted_scan_number=2
) -> tuple:
    """
    Delay the conductivity signal, and recalculate salinity
    Inputs:
        - downcast and upcast data dictionaries, metadata dictionary
        - shifted_scan_number: number of scans shifted. +: delay; -: advance
    Outputs:
        - two dictionaries containing downcast and upcast profiles
    """
    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)
    # Apply the shift to each cast
    for cast_i in var1.keys():
        index_1 = var1[cast_i].Conductivity.index[0]
        v1 = var1[cast_i].Conductivity[index_1]
        index_2 = var2[cast_i].Conductivity.index[0]
        v2 = var2[cast_i].Conductivity[index_2]
        # shift C for n scans
        var1[cast_i].Conductivity = var1[cast_i].Conductivity.shift(
            periods=shifted_scan_number, fill_value=v1
        )
        # calculates SP from C using the PSS-78 algorithm (2 < SP < 42)
        var1[cast_i].Salinity = gsw.SP_from_C(
            var1[cast_i].Conductivity, var1[cast_i].Temperature, var1[cast_i].Pressure
        )
        var2[cast_i].Conductivity = var2[cast_i].Conductivity.shift(
            periods=shifted_scan_number, fill_value=v2
        )
        var2[cast_i].Salinity = gsw.SP_from_C(
            var2[cast_i].Conductivity, var2[cast_i].Temperature, var2[cast_i].Pressure
        )

    metadata_dict["Processing_history"] += (
        "-SHIFT parameters:|"
        " Shift Channel: Conductivity|"
        " # of Records to Delay (-ve for Advance):|"
        f" Shift = {shifted_scan_number}|"
        " Salinity was recalculated after shift|"
    )
    metadata_dict["SHIFT_Conductivity_Time"] = datetime.now()

    return var1, var2


def plot_shift_c(
    cast_d_shift_c: dict,
    cast_u_shift_c: dict,
    cast_d_filtered: dict,
    cast_u_filtered: dict,
    dest_dir: str,
) -> None:
    """
    Plot Salinity and T-S to check the index after shift
    inputs:
        - cast_d_shift_c, cast_u_shift_c: downcast and upcast data dictionaries
        - cast_d_filtered, cast_u_filtered: downcast and upcast data dictionaries
        - dest_dir:
    outputs:
        - profile plots of salinity and T-S plots before and after shifting conductivity and
        recalculating salinity
    """

    # Create a folder for figures if it doesn't already exist
    figure_dir = os.path.join(dest_dir, "FIG")
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    fig, ax = plt.subplots()  # Before
    for cast_i in cast_d_filtered.keys():
        ax.plot(
            cast_d_filtered[cast_i].Salinity,
            cast_d_filtered[cast_i].Pressure,
            color="blue",
        )
        ax.plot(
            cast_u_filtered[cast_i].Salinity,
            cast_u_filtered[cast_i].Pressure,
            color="blue",
        )
    format_processing_plot(
        ax,
        x_var_name="Salinity",
        x_var_units="PSS-78",
        y_var_name="Pressure",
        y_var_units="dbar",
        plot_title="Before Shift Conductivity",
        invert_yaxis=True,
    )
    plt.savefig(os.path.join(figure_dir, "Before_Shift_Conductivity_S.png"))
    plt.close(fig)

    fig, ax = plt.subplots()  # After
    for cast_i in cast_d_filtered.keys():
        ax.plot(
            cast_d_shift_c[cast_i].Salinity,
            cast_d_shift_c[cast_i].Pressure,
            color="blue",
        )
        ax.plot(
            cast_u_shift_c[cast_i].Salinity,
            cast_u_shift_c[cast_i].Pressure,
            color="blue",
        )
    format_processing_plot(
        ax,
        x_var_name="Salinity",
        x_var_units="PSS-78",
        y_var_name="Pressure",
        y_var_units="dbar",
        plot_title="After Shift Conductivity",
        invert_yaxis=True,
    )
    plt.savefig(os.path.join(figure_dir, "After_Shift_Conductivity_S.png"))
    plt.close(fig)

    # TS Plot before
    do_ts_plot(
        figure_dir, "Before Shift Conductivity T-S Plot",
        "Before_Shift_Conductivity_T-S.png", cast_d_filtered, cast_u_filtered,
    )

    # T-S plot After
    do_ts_plot(
        figure_dir, "After Shift Conductivity T-S Plot",
        "After_Shift_Conductivity_T-S.png", cast_d_shift_c, cast_u_shift_c,
    )

    return


def SHIFT_OXYGEN(
    var_downcast: dict, var_upcast: dict, metadata_dict: dict, shifted_scan_number=-11
) -> tuple:
    """
    Advance oxygen data by 2-3s (12-18 scans for 6Hz)
    Inputs:
        - downcast and upcast data dictionaries, metadata dictionary
        - - shifted_scan_number: number of scans shifted. +: delay; -: advance
    Outputs:
        - two dictionaries containing downcast and upcast profiles
    """
    # cast_number = len(var_downcast.keys())
    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)
    for cast_i in var1.keys():
        index_1 = var1[cast_i].Oxygen.index[-1]
        v1 = var1[cast_i].Oxygen[index_1]
        index_2 = var2[cast_i].Oxygen.index[-1]
        v2 = var2[cast_i].Oxygen[index_2]
        # shift C for n scans
        var1[cast_i].Oxygen = var1[cast_i].Oxygen.shift(
            periods=shifted_scan_number, fill_value=v1
        )
        var2[cast_i].Oxygen = var2[cast_i].Oxygen.shift(
            periods=shifted_scan_number, fill_value=v2
        )

    metadata_dict["Processing_history"] += (
        "-SHIFT parameters:|"
        " Shift Channel: Oxygen:Dissolved:Saturation|"
        " # of Records to Delay (-ve for Advance):|"
        f" Shift = {shifted_scan_number}|"
    )
    metadata_dict["SHIFT_Oxygen_Time"] = datetime.now()

    return var1, var2


def plot_shift_o(
    cast_d_shift_o: dict,
    cast_u_shift_o: dict,
    cast_d_shift_c: dict,
    cast_u_shift_c: dict,
    dest_dir: str,
) -> None:
    """
    Check Oxy plots after shift. Plot temperature vs oxygen saturation
    before and after the shift to confirm that the alignment has improved
    inputs:
        - cast_d_shift_o, cast_u_shift_o: downcast and upcast data dictionaries (after o shift)
        - cast_d_shift_c, cast_u_shift_c: downcast and upcast data dictionaries (before o shift)
        - dest_dir
    outputs:
        - Profile plots of oxygen before and after shifting oxygen, but nothing is returned
        by the function
    """

    # Create a folder for figures if it doesn't already exist
    figure_dir = os.path.join(dest_dir, "FIG")
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    # num_casts = len(cast_d_shift_o)

    # Before shift
    fig, ax = plt.subplots()
    for cast_i in cast_d_shift_c.keys():
        ax.plot(
            cast_d_shift_c[cast_i].Oxygen,
            cast_d_shift_c[cast_i].Temperature,
            color="blue",
        )
        ax.plot(
            cast_u_shift_c[cast_i].Oxygen,
            cast_u_shift_c[cast_i].Temperature,
            color="blue",
        )

    format_processing_plot(
        ax,
        x_var_name="Oxygen",
        x_var_units="%",
        y_var_name="Temperature",
        y_var_units="C",
        plot_title="Before Shift Oxygen T-O Plot",
        invert_yaxis=False,
    )
    plt.savefig(os.path.join(figure_dir, "Before_Shift_Oxygen_T-O.png"))
    plt.close(fig)

    # After shift
    fig, ax = plt.subplots()
    for cast_i in cast_d_shift_c.keys():
        ax.plot(
            cast_d_shift_o[cast_i].Oxygen,
            cast_d_shift_o[cast_i].Temperature,
            color="blue",
        )
        ax.plot(
            cast_u_shift_o[cast_i].Oxygen,
            cast_u_shift_o[cast_i].Temperature,
            color="blue",
        )

    format_processing_plot(
        ax,
        x_var_name="Oxygen",
        x_var_units="%",
        y_var_name="Temperature",
        y_var_units="C",
        plot_title="After Shift Oxygen T-O Plot",
        invert_yaxis=False,
    )
    plt.savefig(os.path.join(figure_dir, "After_Shift_Oxygen_T-O.png"))
    plt.close(fig)
    return


def DERIVE_OXYGEN_CONCENTRATION(
    var_downcast: dict, var_upcast: dict, metadata_dict: dict
) -> tuple:
    """
    Derive oxygen concentration in umol/kg and mL/L from oxygen percent saturation
    using SCOR WG 142 (DOI:10.13155/45915)
    inputs:
        - var_downcast, var_upcast: downcast and upcast data dictionaries
        - metadata_dict
    outputs:
        - downcast and upcast data dictionaries with derived oxygen concentration
        variables added
    """
    umol_L_to_mL_L = 1 / 44.6596
    m3_to_L = 1e3
    # cast_number = len(var_downcast.keys())
    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)

    # O_sat_num_decimal_places = None

    for cast_i in var1.keys():
        for var in [var1, var2]:
            T = var[cast_i].Temperature.to_numpy()
            S = var[cast_i].Salinity.to_numpy()
            O_sat = var[cast_i].Oxygen.to_numpy()
            P = var[cast_i].Pressure.to_numpy()

            # Convert oxygen saturation to molar oxygen concentration
            # Use default P=0dbar and p_atm=1013.25 mbar, where
            # P: hydrostatic pressure in dbar, and
            # p_atm: atmospheric (air) pressure in mbar
            O_umol_L = O2stoO2c(O_sat, T, S)
            # Convert to mL/L
            O_mL_L = O_umol_L * umol_L_to_mL_L
            # Convert to umol/kg using potential density of seawater (kg/L) from
            # Fofonoff and Millard (1983) and Millero et al. (1980).
            # rho: potential density of seawater referenced to a hydrostatic pressure
            # of 0 dbar and using practical salinity.
            # Since TEOS-10 is based off ABSOLUTE salinity, it can't be used here!
            # One option is seawater.eos80.pden() potential density
            # Returns potential density relative to the ref. pressure [kg m :sup:3]
            rho_kg_m3 = eos80.pden(S, T, P)
            rho_kg_L = rho_kg_m3 / m3_to_L
            O_umol_kg = O_umol_L / rho_kg_L
            var[cast_i]["Oxygen_mL_L"] = O_mL_L
            var[cast_i]["Oxygen_umol_kg"] = O_umol_kg

    metadata_dict["Processing_history"] += (
        "-Oxygen concentration was calculated from oxygen "
        "saturation using SCOR WG 142|"
    )
    metadata_dict["DERIVE_OXYGEN_CONCENTRATION_Time"] = datetime.now()
    return var1, var2


# def pH2O_Weiss_Price(T, S):
#     """
#     Compute the vapour pressure of water from temperature and salinity following
#     Weiss and Price (1980; doi:10.1016/0304-4203(80)90024-9)
#     inputs:
#         - T: temperature in degrees Celsius
#         - S: salinity (dimensionless, Practical Salinity Scale 1978)
#     outputs:
#         - pH2O: vapour pressure of water with units of atm
#     """
#     D0 = 24.4543
#     D1 = -67.4509
#     D2 = -4.8489
#     D3 = -5.44e-4
#     T_abs = T + 273.15
#     return np.exp(D0 + D1 * (100/T_abs) + D2 * np.log(T_abs/100) + D3 * S)


# def T_corr_Garcia_Gordon(T):
#     """
#     Correction term; the temperature-dependent part of seawater O2 solubility.
#     From Garcia and Gordon 1992, Benson and Krause refit
#     """
#     A0 = 2.00907
#     A1 = 3.22014
#     A2 = 4.05010
#     A3 = 4.94457
#     A4 = -2.56847e-1
#     A5 = 3.88767  # (Garcia and Gordon 1992, Benson and Krause refit)
#
#     T_s = np.log((298.15 - T) / (273.15 + T))
#
#     return 44.6596 * np.exp(A0 + A1*T_s + A2*T_s**2 + A3*T_s**3 + A4*T_s**4 + A5*T_s**5)


# def S_corr_Garcia_Gordon(T, S):
#     """
#     Correction term; the salinity-dependent part of seawater O2 solubility.
#     From Garcia and Gordon 1992, Benson and Krause refit
#     """
#     B0 = -6.24523e-3
#     B1 = -7.37614e-3
#     B2 = -1.03410e-2
#     B3 = -8.17083e-3
#     C0 = -4.88682e-7
#     S_preset = 0
#
#     T_s = np.log((298.15 - T) / (273.15 + T))
#
#     return np.exp((S - S_preset) *
#                   (B0 + B1*T_s + B2*T_s**2 + B3 * T_s**3) +
#                   C0 * (S**2 - S_preset**2))


# def DERIVE_OXYGEN_ML_L(var_downcast, var_upcast, metadata_dict):
#     """
#     Convert oxygen percent saturation to concentration with mL/L units
#     using SCOR WG 142 (DOI:10.13155/45915) equation used by ARGO program.
#     Hakai uses this reference to derive oxygen concentration.
#     Inputs:
#          - downcast and upcast data dictionaries, metadata dictionary
#     Outputs:
#          - two dictionaries containing downcast and upcast profiles
#     """
#     # test this function
#
#     umol_L_to_mL_L = 1/44.6596
#
#     cast_number = len(var_downcast.keys())
#     var1 = deepcopy(var_downcast)
#     var2 = deepcopy(var_upcast)
#
#     for i in range(1, cast_number + 1, 1):
#         for var in [var1, var2]:
#             T = var['cast' + str(i)].Temperature.values
#             S = var['cast' + str(i)].Salinity.values
#             O_sat = var['cast' + str(i)].Oxygen.values
#
#             T_corr = T_corr_Garcia_Gordon(T)
#             S_corr = S_corr_Garcia_Gordon(T, S)
#
#             cO2_star = T_corr * S_corr  # Oxygen solubility
#             cO2_umol_L = cO2_star * O_sat  # Oxygen concentration in what units?
#             cO2_mL_L = cO2_umol_L * umol_L_to_mL_L  # Convert to the desired units
#             var['cast' + str(i)]['Oxygen_concentration'] = cO2_mL_L  # Add to data dictionary
#
#     if 'Processing_history' not in metadata_dict:  # For testing purposes
#         metadata_dict['Processing_history'] = ''
#
#     metadata_dict['Processing_history'] += 'Converted Dissolved O2 % to mL/L units by using SCOR WG 142 ' \
#                                            '(DOI:10.13155/45915) saturation concentrations equation with ' \
#                                            'CTD and DO data 4.1667s smoothed|'
#     metadata_dict['DERIVE_OXYGEN_ML_L_Time'] = datetime.now()
#
#     return var1, var2


def DELETE_PRESSURE_REVERSAL(
    var_downcast: dict, var_upcast: dict, metadata_dict: dict
) -> tuple:
    """
    Detect and delete pressure reversal (swells/slow drop),
    correct for the wake effect.
    Inputs:
        - downcast and upcast data dictionaries
        - metadata dictionary
    Outputs:
        - two dictionaries containing downcast and upcast profiles
    """
    # cast_number = len(var_downcast.keys())
    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)
    for cast_i in var1.keys():
        press = var1[cast_i].Pressure.values
        ref = press[0]
        inversions = np.diff(np.r_[press, press[-1]]) < 0  # a mask
        mask = np.zeros_like(inversions)
        for k, p in enumerate(inversions):
            if p:
                ref = press[k]
                cut = press[k + 1 :] < ref
                mask[k + 1 :][cut] = True
        var1[cast_i][mask] = np.NaN

    for cast_i in var2.keys():
        press = var2[cast_i].Pressure.values
        ref = press[0]
        inversions = np.diff(np.r_[press, press[-1]]) > 0
        mask = np.zeros_like(inversions)
        for k, p in enumerate(inversions):
            if p:
                ref = press[k]
                cut = press[k + 1 :] > ref
                mask[k + 1 :][cut] = True
        var2[cast_i][mask] = np.NaN
    metadata_dict["Processing_history"] += (
        "-DELETE_PRESSURE_REVERSAL parameters:|" " Remove pressure reversals|"
    )
    metadata_dict["DELETE_PRESSURE_REVERSAL_Time"] = datetime.now()

    return var1, var2


def plot_delete(cast_d_wakeeffect: dict, cast_d_shift_o: dict, dest_dir: str) -> None:
    """
    Plot downcast profiles before and after removing pressure reversals (delete step)
    Plot all variables including derived oxygen concentration if available
    inputs:
        - cast_d_wakeeffect: dictionary containing downcast data after the delete step
        - cast_d_shift_o: dictionary containing downcast data before the delete step
        - dest_dir: working and output directory
    outputs:
        - profile plots of all available variables plus T-S plots comparing before and
        after delete are saved to dest_dir/FIG/, but nothing is returned by the function
    """
    # Create a folder for figures if it doesn't already exist
    figure_dir = os.path.join(dest_dir, "FIG")
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    # Get a list of the available variables
    for cast_i in cast_d_wakeeffect.keys():
        vars_available = list(dict.fromkeys(cast_d_wakeeffect[cast_i]))
        break

    # Plot BEFORE: Iterate through the casts and variables
    for j, var in enumerate(VARIABLES_POSSIBLE):
        if var in vars_available:
            fig, ax = plt.subplots()
            for cast_i in cast_d_shift_o.keys():
                ax.plot(
                    cast_d_shift_o[cast_i].loc[:, var],
                    cast_d_shift_o[cast_i].Pressure,
                    color=VARIABLE_COLOURS[j],
                )
            format_processing_plot(
                ax,
                x_var_name=var,
                x_var_units=VARIABLE_UNITS[j],
                y_var_name="Pressure",
                y_var_units="dbar",
                plot_title="Before Delete",
                invert_yaxis=True,
            )

            if len(var.split("_")) > 1:
                var_abbrev = var[0] + "_" + var.split("_")[1] + "_" + var.split("_")[2]
            else:
                var_abbrev = var[0]
            plt.savefig(os.path.join(figure_dir, f"Before_Delete_{var_abbrev}.png"))
            plt.close(fig)

    # Iterate through the casts and variables
    for j, var in enumerate(VARIABLES_POSSIBLE):
        if var in vars_available:
            fig, ax = plt.subplots()
            for cast_i in cast_d_wakeeffect.keys():
                ax.plot(
                    cast_d_wakeeffect[cast_i].loc[:, var],
                    cast_d_wakeeffect[cast_i].Pressure,
                    color=VARIABLE_COLOURS[j],
                )
            format_processing_plot(
                ax,
                x_var_name=var,
                x_var_units=VARIABLE_UNITS[j],
                y_var_name="Pressure",
                y_var_units="dbar",
                plot_title="After Delete",
                invert_yaxis=True,
            )
            if len(var.split("_")) > 1:
                # Capture derived oxygen concentration without overwriting O sat plot
                var_abbrev = var[0] + "_" + var.split("_")[1] + "_" + var.split("_")[2]
            else:
                var_abbrev = var[0]
            plt.savefig(os.path.join(figure_dir, f"After_Delete_{var_abbrev}.png"))
            plt.close(fig)

    # TS Plots
    do_ts_plot(
        figure_dir, plot_title="T-S Plot (before delete pressure reversal)",
        figure_filename="Before_Delete_T-S.png", cast_d=cast_d_shift_o
    )
    do_ts_plot(
        figure_dir, plot_title="T-S Plot (after delete pressure reversal)",
        figure_filename="After_Delete_T-S.png", cast_d=cast_d_wakeeffect
    )
    plt.savefig(os.path.join(figure_dir, "Before_Delete_T-S.png"))
    plt.close(fig)

    fig, ax = plt.subplots()  # After
    for cast_i in cast_d_wakeeffect.keys():
        ax.plot(
            cast_d_wakeeffect[cast_i].Salinity,
            cast_d_wakeeffect[cast_i].Temperature,
            color="blue",
        )
    format_processing_plot(
        ax,
        x_var_name="Salinity",
        x_var_units="PSS-78",
        y_var_name="Temperature",
        y_var_units="C",
        plot_title="T-S Plot (after delete pressure reversal)",
        invert_yaxis=False,
    )
    plt.savefig(os.path.join(figure_dir, "After_Delete_T-S.png"))
    plt.close(fig)
    return


def BINAVE(
    var_downcast: dict, var_upcast: dict, metadata_dict: dict, interval=1
) -> tuple:
    """
    Bin average the profiles
    Note: Bin width and spacing are both universally chosen to be 1m in coastal waters
    Inputs:
        - downcast and upcast data dictionaries, metadata dictionary
        - default bin average interval set to 1 dbar
    Outputs:
        - two dictionaries containing downcast and upcast profiles
    """
    # cast_number = len(var_downcast.keys())
    var1 = deepcopy(var_downcast)
    var2 = deepcopy(var_upcast)
    # Iterate through all the casts
    for cast_i in var1.keys():
        start_d = np.floor(np.nanmin(var1[cast_i].Pressure.values))
        # start_d = np.round(start_d)
        stop_d = np.ceil(np.nanmax(var1[cast_i].Pressure.values))
        # stop_d = np.round(stop_d)
        new_press_d = np.arange(start_d - 0.5, stop_d + 1.5, interval)
        binned_d = pd.cut(var1[cast_i].Pressure, bins=new_press_d)
        obs_count_d = var1[cast_i].groupby(binned_d).size()

        var1[cast_i] = var1[cast_i].groupby(binned_d).mean()
        var1[cast_i]["Observation_counts"] = obs_count_d
        # Potential for whole row Nan values at top and bottom of output files
        var1[cast_i] = var1[cast_i].dropna(
            axis=0, how="any"
        )  # drop the nans - ask if this is OK?
        var1[cast_i].reset_index(drop=True, inplace=True)

        start_u = np.ceil(np.nanmax(var2[cast_i].Pressure.values))
        stop_u = np.floor(np.nanmin(var2[cast_i].Pressure.values))
        new_press_u = np.arange(start_u + 0.5, stop_u - 1.5, -interval)
        binned_u = pd.cut(var2[cast_i].Pressure, bins=new_press_u[::-1])
        obs_count_u = var2[cast_i].groupby(binned_u).size()

        var2[cast_i] = var2[cast_i].groupby(binned_u).mean()
        var2[cast_i] = var2[cast_i].sort_values("Depth", ascending=False)
        var2[cast_i]["Observation_counts"] = obs_count_u
        var2[cast_i] = var2[cast_i].dropna(axis=0, how="any")
        var2[cast_i].reset_index(drop=True, inplace=True)

    metadata_dict["Processing_history"] += (
        "-BINAVE parameters:"
        " Bin channel = Pressure|"
        " Averaging interval = 1.00|"
        " Minimum bin value = 0.000|"
        " Average value were used|"
        " Interpolated values were NOT used for empty bins|"
        " Channel NUMBER_OF_BIN_RECORDS was added to file|"
    )
    metadata_dict["BINAVE_Time"] = datetime.now()

    return var1, var2


def FINAL_EDIT(
    var_cast: dict, have_oxy: bool, have_fluor: bool, metadata_dict: dict
) -> dict:
    """
    Final editing the profiles: edit header information, correct the unit of conductivity
    Inputs:
        - downcast and upcast data dictionaries, metadata dictionary
    Outputs:
        - two dictionaries containing downcast and upcast profiles
    """

    # vars = list(dict.fromkeys(var_cast['cast1']))
    # cast_number = len(var_cast.keys())
    var = deepcopy(var_cast)

    if have_oxy and have_fluor:
        col_list = [
            "Pressure",
            "Depth",
            "Temperature",
            "Salinity",
            "Fluorescence",
            "Oxygen",
            "Oxygen_mL_L",
            "Oxygen_umol_kg",
            "Conductivity",
            "Observation_counts",
        ]
    elif have_oxy and not have_fluor:
        col_list = [
            "Pressure",
            "Depth",
            "Temperature",
            "Salinity",
            "Oxygen",
            "Oxygen_mL_L",
            "Oxygen_umol_kg",
            "Conductivity",
            "Observation_counts",
        ]
    elif have_fluor and not have_oxy:
        col_list = [
            "Pressure",
            "Depth",
            "Temperature",
            "Salinity",
            "Fluorescence",
            "Conductivity",
            "Observation_counts",
        ]
    else:
        col_list = [
            "Pressure",
            "Depth",
            "Temperature",
            "Salinity",
            "Conductivity",
            "Observation_counts",
        ]

    # Do channel format corrections, unit conversions for each cast
    for cast_i in var.keys():
        var[cast_i] = var[cast_i].reset_index(drop=True)  # drop index column
        var[cast_i] = var[cast_i][col_list]  # select columns
        var[cast_i].Conductivity = (
            var[cast_i].Conductivity * 0.1
        )  # convert Conductivity to S/m

        var[cast_i].Pressure = var[cast_i].Pressure.apply("{:,.1f}".format)
        var[cast_i].Depth = var[cast_i].Depth.apply("{:,.1f}".format)
        var[cast_i].Temperature = var[cast_i].Temperature.apply("{:,.4f}".format)
        var[cast_i].Salinity = var[cast_i].Salinity.apply("{:,.4f}".format)
        if have_fluor:
            var[cast_i].Fluorescence = var[cast_i].Fluorescence.apply("{:,.3f}".format)
        if have_oxy:
            var[cast_i].Oxygen = var[cast_i].Oxygen.apply("{:,.2f}".format)
            var[cast_i].Oxygen_mL_L = var[cast_i].Oxygen_mL_L.apply("{:,.2f}".format)
            var[cast_i].Oxygen_umol_kg = var[cast_i].Oxygen_umol_kg.apply(
                "{:,.2f}".format
            )
        var[cast_i].Conductivity = var[cast_i].Conductivity.apply("{:,.6f}".format)

        # HH DO NOT rename columns because would have repeating oxygen concentration names
        # # Rename the columns in each cast taking into account the existing order of columns
        # if have_oxy and have_fluor:
        #     var['cast' + str(i)].columns = ['Pressure', 'Depth', 'Temperature', 'Salinity',
        #                                     'Fluorescence:URU', 'Oxygen:Dissolved:Saturation:RBR',
        #                                     'Conductivity', 'Number_of_bin_records']
        # elif have_oxy and not have_fluor:
        #     var['cast' + str(i)].columns = ['Pressure', 'Depth', 'Temperature', 'Salinity',
        #                                     'Oxygen:Dissolved:Saturation:RBR',
        #                                     'Conductivity', 'Number_of_bin_records']
        # elif have_fluor and not have_oxy:
        #     var['cast' + str(i)].columns = ['Pressure', 'Depth', 'Temperature', 'Salinity',
        #                                     'Fluorescence:URU',
        #                                     'Conductivity', 'Number_of_bin_records']
        # else:
        #     var['cast' + str(i)].columns = ['Pressure', 'Depth', 'Temperature', 'Salinity',
        #                                     'Conductivity', 'Number_of_bin_records']

    # Update the Processing history
    metadata_dict["Processing_history"] += (
        "-Remove Channels:|"
        " The following CHANNEL(S) were removed:|"
        " Date|"
        " TIME:UTC|"
        "-CALIB parameters:|"
        " Calobration type = Correct|"
        " Calibration applied:|"
        " Conductivity (S/m) = 0.1* Conductivity (mS/cm)|"
    )
    metadata_dict["FINALEDIT_Time"] = datetime.now()

    return var


def plot_processed(cast_final: dict, dest_dir: str) -> None:
    """Final plots: Plot the processed casts after bin averaging
    inputs:
        - cast_final: output from FINAL_EDIT() function
        - dest_dir: working and output directory
    outputs:
        - profile plots of all available channels plus a T-S plot
    """
    # Create a folder for figures if it doesn't already exist
    figure_dir = os.path.join(dest_dir, "FIG")
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    for cast_i in cast_final.keys():
        vars_available = list(dict.fromkeys(cast_final[cast_i]))
        break

    # -------------------------  Plot processed profiles ------------------------------
    # number_of_colors = len(cast)
    # color = ["#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    #          for i in range(number_of_colors)]

    # Change conductivity units from mS/cm to agree with FINAL_EDIT()
    variable_units_final = VARIABLE_UNITS
    variable_units_final[2] = "S/m"

    # Plot the first cast before and after processing
    for j, var in enumerate(VARIABLES_POSSIBLE):
        if var in vars_available:
            fig, ax = plt.subplots()
            for cast_i in cast_final.keys():
                # Need to convert to float because FINAL_EDIT() converted values to strings
                ax.plot(
                    cast_final[cast_i].loc[:, var].astype(float),
                    cast_final[cast_i].Pressure.astype(float),
                    color=VARIABLE_COLOURS[j],
                )
            format_processing_plot(
                ax,
                x_var_name=var,
                x_var_units=variable_units_final[j],
                y_var_name="Pressure",
                y_var_units="dbar",
                plot_title="Post-Processing",
                invert_yaxis=True,
            )
            if len(var.split("_")) > 1:
                # Capture derived oxygen concentration without overwriting O sat plot: O_ml_l and O_umol_kg
                var_abbrev = var[0] + "_" + var.split("_")[1] + "_" + var.split("_")[2]
            else:
                var_abbrev = var[0]
            plt.savefig(os.path.join(figure_dir, f"Post_Processing_{var_abbrev}.png"))
            plt.close(fig)

    # T-S plot
    do_ts_plot(
        figure_dir, "Post-Processing T-S Plot", "Post_Processing_T-S.png", cast_final
    )

    return


def write_file(
        cast_number,
        cast_original: dict,
        cast_final: dict,
        metadata_dict: dict,
) -> None:
    """
    Write file section of IOS header file
    Inputs:
        - cast_number, cast_original = cast,
          cast_final = cast_d_final,
          metadata_dict = metadata
        - have_fluor, have_oxy: boolean flags; True if fluorescence and oxygen data
        are present in the dataset, respectively
    Outputs:
        - File information added to IOS header file that is open, but nothing returned
        by the function
    """

    # vars = list(dict.fromkeys(cast_original['cast1']))

    start_time = pd.to_datetime(
        cast_original["cast" + str(cast_number)].Date.values[0]
        + " "
        + cast_original["cast" + str(cast_number)].TIME.values[0]
    ).strftime("%Y/%m/%d %H:%M:%S.%f")[0:-3]
    end_time = pd.to_datetime(
        cast_original["cast" + str(cast_number)].Date.values[-1]
        + " "
        + cast_original["cast" + str(cast_number)].TIME.values[-1]
    ).strftime("%Y/%m/%d %H:%M:%S.%f")[0:-3]

    sample_interval = metadata_dict["Sampling_Interval"]
    time_increment = "0 0 0 " + sample_interval + " 0  ! (day hr min sec ms)"

    # Number of ensembles
    number_of_records = str(cast_final["cast" + str(cast_number)].shape[0])
    data_description = metadata_dict["Data_description"]
    number_of_channels = str(cast_final["cast" + str(cast_number)].shape[1])
    nan = -99
    file_type = "ASCII"

    print("*FILE")
    print("    " + "{:20}".format("START TIME") + ": UTC " + start_time)
    print("    " + "{:20}".format("END TIME") + ": UTC " + end_time)
    print("    " + "{:20}".format("TIME INCREMENT") + ": " + time_increment)
    print("    " + "{:20}".format("NUMBER OF RECORDS") + ": " + number_of_records)
    print("    " + "{:20}".format("DATA DESCRIPTION") + ": " + data_description)
    print("    " + "{:20}".format("FILE TYPE") + ": " + file_type)
    print("    " + "{:20}".format("NUMBER OF CHANNELS") + ": " + number_of_channels)
    print()
    print("{:>20}".format("$TABLE: CHANNELS"))
    print(
        "    "
        + "! No Name                             Units          Minimum        Maximum"
    )
    print(
        "    "
        + "!--- -------------------------------- -------------- -------------- --------------"
    )

    # Create a dictionary of channel information
    # Don't use the channel names as the dict keys because Oxygen:Dissolved:RBR repeats
    # format: {df_column: (ios_channel_name, channel_unit, format_channel_max, channel_pad,
    #                   channel_width, channel_format, channel_type, channel_decimal_places)}

    # Column names no longer changed in FINAL_EDIT()
    channel_dict = {
        "Pressure": ("Pressure", "decibar", False, nan, 7, "F", "R4", 1),
        "Depth": ("Depth", "metres", False, nan, 7, "F", "R4", 1),
        "Temperature": ("Temperature", "'deg C (ITS90)'", False, nan, 9, "F", "R4", 4),
        "Salinity": ("Salinity", "PSS-78", "%.04f", nan, 9, "F", "R4", 4),
        "Fluorescence": ("Fluorescence:URU", "mg/m^3", "%.03f", nan, 8, "F", "R4", 3),
        "Oxygen": (
            "Oxygen:Dissolved:Saturation:RBR",
            "%",
            "%.04f",
            nan,
            8,
            "F",
            "R4",
            2,
        ),
        "Oxygen_mL_L": ("Oxygen:Dissolved:RBR", "mL/L", "%.04f", nan, 8, "F", "R4", 2),
        "Oxygen_umol_kg": (
            "Oxygen:Dissolved:RBR",
            "umol/kg",
            "%.04f",
            nan,
            8,
            "F",
            "R4",
            2,
        ),
        "Conductivity": ("Conductivity", "S/m", "%.05f", nan, 10, "F", "R4", 6),
        "Observation_counts": (
            "Number_of_bin_records",
            "n/a",
            False,
            "' '",
            5,
            "I",
            "I",
            0,
        ),
    }

    # Remove unavailable channels
    for channel in VARIABLES_POSSIBLE:
        if channel not in cast_final[f"cast{cast_number}"]:
            _ = channel_dict.pop(channel)

    current_chan_no = 1  # Update current channel number as iteration progresses
    for i, key in enumerate(channel_dict):
        df_name = key
        # Unpack the first three elements of the tuple to more descriptive names
        ios_name, unit, format_max = channel_dict[key][:3]
        # Add this variable since number of bin records is the only int type channel, the rest are float
        dtype_minmax = float if channel_dict[key][5] == "F" else int

        if format_max:
            # The min is not formatted but the max is
            # nanmin and nanmax will still include -99 pad values in computation????
            print(
                "{:>8}".format(str(current_chan_no))
                + " "
                + "{:33}".format(ios_name)
                + "{:15}".format(unit)
                + "{:15}".format(
                    str(
                        np.nanmin(
                            cast_final[f"cast{cast_number}"]
                            .loc[:, df_name]
                            .astype(dtype_minmax)
                        )
                    )
                )
                + "{:14}".format(
                    str(
                        float(
                            format_max
                            % np.nanmax(
                                cast_final[f"cast{cast_number}"]
                                .loc[:, df_name]
                                .astype(dtype_minmax)
                            )
                        )
                    )
                )
            )
        else:
            print(
                "{:>8}".format(str(current_chan_no))
                + " "
                + "{:33}".format(ios_name)
                + "{:15}".format(unit)
                + "{:15}".format(
                    str(
                        np.nanmin(
                            cast_final[f"cast{cast_number}"]
                            .loc[:, df_name]
                            .astype(dtype_minmax)
                        )
                    )
                )
                + "{:14}".format(
                    str(
                        np.nanmax(
                            cast_final[f"cast{cast_number}"]
                            .loc[:, df_name]
                            .astype(dtype_minmax)
                        )
                    )
                )
            )
        # Update channel counter
        current_chan_no += 1

    # Add in table of Channel summary
    print("{:>8}".format("$END"))
    print()
    print("{:>26}".format("$TABLE: CHANNEL DETAIL"))
    print("    " + "! No  Pad   Start  Width  Format  Type  Decimal_Places")
    print("    " + "!---  ----  -----  -----  ------  ----  --------------")
    # print('{:>8}'.format('1') + "  " + '{:15}'.format("' '") + '{:7}'.format(' ') +
    # '{:7}'.format("' '") + '{:22}'.format('YYYY-MM-DDThh:mm:ssZ') +
    # '{:6}'.format('D, T') + '{:14}'.format("' '"))

    current_chan_no = 1
    for key, value in channel_dict.items():
        # Unpack the tuple to more descriptive names
        pad, width, chan_format, chan_type, decimals = value[3:]
        print(
            "{:>8}".format(str(current_chan_no))
            + "  "
            + "{:6}".format(str(pad))
            + "{:7}".format("' '")
            + "{:7}".format(str(width))
            + "{:8}".format(chan_format)
            + "{:6}".format(chan_type)
            + "{:3}".format(decimals)
        )
        # Update channel counter
        current_chan_no += 1

    # Add in table of Channel detail summary
    print("{:>8}".format("$END"))
    print()
    return


def write_admin(metadata_dict: dict) -> None:
    """
    function to write administation section of IOS header file
    inputs:
        - metadata_dict: dictionary containing metadata for the RBR CTD
    outputs:
        - Administration section added to open IOS header file, but nothing returned
        by the function
    """
    mission = metadata_dict["Mission"]
    agency = metadata_dict["Agency"]
    country = metadata_dict["Country"]
    project = metadata_dict["Project"]
    scientist = metadata_dict["Scientist"]
    platform = metadata_dict["Platform"]
    print("*ADMINISTRATION")
    print("    " + "{:20}".format("MISSION") + ": " + mission)
    print("    " + "{:20}".format("AGENCY") + ": " + agency)
    print("    " + "{:20}".format("COUNTRY") + ": " + country)
    print("    " + "{:20}".format("PROJECT") + ": " + project)
    print("    " + "{:20}".format("SCIENTIST") + ": " + scientist)
    print("    " + "{:20}".format("PLATFORM ") + ": " + platform)
    print()
    return


def write_location(cast_number: int, metadata_dict: dict) -> None:
    """
    write location part in IOS header file
    Inputs:
        - cast_number
        - metadata_dict: dictionary containing metadata for the RBR CTD
    Outputs:
        - location section added to open IOS header file, but nothing returned
        by the function
    """
    station = metadata_dict["Location"]["LOC:STATION"].to_numpy()
    event_number = metadata_dict["Location"]["LOC:Event Number"].to_numpy()
    lon = metadata_dict["Location"]["LOC:LONGITUDE"].to_numpy()
    lat = metadata_dict["Location"]["LOC:LATITUDE"].to_numpy()
    water_depth = metadata_dict["Location"]["LOC:Water Depth"].to_numpy()

    event_mask = event_number == cast_number

    station = str(station[event_mask][0])
    event_number = str(event_number[event_mask][0])
    lon = lon[event_mask][0].split(" ")
    lat = lat[event_mask][0].split(" ")
    water_depth = str(water_depth[event_mask][0])

    # Correct lat and lon formatting
    # Fill lat and lon degree placement to 3 characters with spaces if needed
    lat[0] = f"  {lat[0]}"[-3:]
    lon[0] = f"  {lon[0]}"[-3:]
    # Remove filler zero from minutes part of coordinates
    if lat[1][0] == "0":
        lat[1] = lat[1][1:]
    if lon[1][0] == "0":
        lon[1] = lon[1][1:]

    print("*LOCATION")
    # print("    " + '{:20}'.format('STATION') + ": " + str(station_number[cast_number - 1]))
    # print("    " + '{:20}'.format('EVENT NUMBER') + ": " + str(event_number[cast_number - 1]))
    # print("    " + '{:20}'.format('LATITUDE') + ":  " + lat[cast_number - 1][0:10] +
    #       "0 " + lat[cast_number - 1][-14:-1] + ")")
    # print("    " + '{:20}'.format('LONGITUDE') + ": " + lon[cast_number - 1])
    print("    " + "{:20}".format("STATION") + ": " + station)
    print("    " + "{:20}".format("EVENT NUMBER") + ": " + str(event_number))
    # print("    " + '{:20}'.format('LATITUDE') + ":  " + lat[0:10] + "0" + lat[-14:])
    # print("    " + '{:20}'.format('LONGITUDE') + ": " + lon[0:11] + "0" + lon[-14:])
    print(
        "    "
        + "{:20}".format("LATITUDE")
        + ": "
        + lat[0]
        + " "
        + lat[1]
        + "0 "
        + lat[2]
        + "  ! (deg min)"
    )
    print(
        "    "
        + "{:20}".format("LONGITUDE")
        + ": "
        + lon[0]
        + " "
        + lon[1]
        + "0 "
        + lon[2]
        + "  ! (deg min)"
    )
    print("    " + "{:20}".format("WATER DEPTH") + ": " + water_depth)
    print()
    return


def write_instrument(metadata_dict: dict) -> None:
    """function to write instrument info
    inputs:
        - metadata_dict: dictionary containing metadata for the RBR CTD
    outputs:
        - Instrument section added to open IOS header file, but nothing returned
        by the function
    """
    model = metadata_dict["Instrument_Model"]
    if int(metadata_dict["Serial_number"]) < 1000:
        serial_number = f"{0:0}" + metadata_dict["Serial_number"]
    else:
        serial_number = metadata_dict["Serial_number"]
    data_description = metadata_dict["Data_description"]
    instrument_type = metadata_dict["Instrument_type"]
    print("*INSTRUMENT")
    print("    MODEL               : " + model)
    print("    SERIAL NUMBER       : " + serial_number)
    print(
        "    DATA DESCRIPTION    : "
        + data_description
        + "                               ! custom item"
    )
    print(
        "    INSTRUMENT TYPE     : "
        + instrument_type
        + "                           ! custom item"
    )
    print()
    return


def write_history(
        cast_original: dict,
        cast_correct_time: dict,
        cast_calib: dict,
        cast_clip: dict,
        cast_filtered: dict,
        cast_shift_c: dict,
        cast_shift_o: dict,
        cast_d_o_conc: dict,
        cast_wakeeffect: dict,
        cast_binned: dict,
        cast_dropvars: dict,
        cast_final: dict,
        cast_number: int,
        metadata_dict: dict,
):
    """
    function to write raw info
    inputs:
        - have_fluor: boolean flag to indicate if fluorescence channel is available
        - have_oxy: boolean flag to indicate if oxygen channels are available
        (saturation + derived concentration)
        - cast_original:
        - cast_clip:
        - cast_filtered:
        - cast_shift_c:
        - cast_shift_o: dict or None if no oxygen sensor
        - cast_wakeeffect:
        - cast_binned:
        - cast_final:
        - cast_number:
        - metadata_dict: dictionary containing metadata for the RBR CTD
    outputs:
        - history section added to open IOS header file, but nothing returned
        by the function
    """

    time_format = "%Y/%m/%d %H:%M:%S.%f"
    print("*HISTORY")
    print()
    print("    $TABLE: PROGRAMS")
    print("    !   Name     Vers   Date       Time     Recs In   Recs Out")
    print("    !   -------- ------ ---------- -------- --------- ---------")
    print(
        "        Z ORDER  "
        + "{:7}".format(str(1.0))
        + "{:11}".format(
            metadata_dict["ZEROORDER_Time"].strftime(time_format)[0:-7].split(" ")[0]
        )
        + "{:9}".format(
            metadata_dict["ZEROORDER_Time"].strftime(time_format)[0:-7].split(" ")[1]
        )
        + "{:>9}".format(str(cast_original["cast" + str(cast_number)].shape[0]))
        + "{:>10}".format(str(cast_original["cast" + str(cast_number)].shape[0]))
    )
    if "CORRECT_TIME_OFFSET_Time" in metadata_dict.keys():
        print(
            "        CORRECT_T"
            + "{:7}".format(str(1.0))
            + "{:11}".format(
                metadata_dict["CORRECT_TIME_OFFSET_Time"].strftime(time_format)[0:-7].split(" ")[0]
            )
            + "{:9}".format(
                metadata_dict["CORRECT_TIME_OFFSET_Time"].strftime(time_format)[0:-7].split(" ")[1]
            )
            + "{:>9}".format(str(cast_original["cast" + str(cast_number)].shape[0]))
            + "{:>10}".format(str(cast_correct_time["cast" + str(cast_number)].shape[0]))
        )
    if "CALIB_Time" in metadata_dict.keys():
        print(
            "        CALIB    "
            + "{:7}".format(str(1.0))
            + "{:11}".format(
                metadata_dict["CALIB_Time"].strftime(time_format)[0:-7].split(" ")[0]
            )
            + "{:9}".format(
                metadata_dict["CALIB_Time"].strftime(time_format)[0:-7].split(" ")[1]
            )
            + "{:>9}".format(str(cast_correct_time["cast" + str(cast_number)].shape[0]))
            + "{:>10}".format(str(cast_calib["cast" + str(cast_number)].shape[0]))
        )
    print(
        "        CLIP     "
        + "{:7}".format(str(1.0))
        + "{:11}".format(
            metadata_dict["CLIP_D_Time" + str(cast_number)]
                .strftime(time_format)[0:-7]
                .split(" ")[0]
        )
        + "{:9}".format(
            metadata_dict["CLIP_D_Time" + str(cast_number)]
                .strftime(time_format)[0:-7]
                .split(" ")[1]
        )
        + "{:>9}".format(str(cast_calib["cast" + str(cast_number)].shape[0]))
        + "{:>10}".format(str(cast_clip["cast" + str(cast_number)].shape[0]))
    )
    print(
        "        FILTER   "
        + "{:7}".format(str(1.0))
        + "{:11}".format(
            metadata_dict["FILTER_Time"].strftime(time_format)[0:-7].split(" ")[0]
        )
        + "{:9}".format(
            metadata_dict["FILTER_Time"].strftime(time_format)[0:-7].split(" ")[1]
        )
        + "{:>9}".format(str(cast_clip["cast" + str(cast_number)].shape[0]))
        + "{:>10}".format(str(cast_filtered["cast" + str(cast_number)].shape[0]))
    )
    print(
        "        SHIFT    "
        + "{:7}".format(str(1.0))
        + "{:11}".format(
            metadata_dict["SHIFT_Conductivity_Time"]
                .strftime(time_format)[0:-7]
                .split(" ")[0]
        )
        + "{:9}".format(
            metadata_dict["SHIFT_Conductivity_Time"]
                .strftime(time_format)[0:-7]
                .split(" ")[1]
        )
        + "{:>9}".format(str(cast_filtered["cast" + str(cast_number)].shape[0]))
        + "{:>10}".format(str(cast_shift_c["cast" + str(cast_number)].shape[0]))
    )
    if "SHIFT_Oxygen_Time" in metadata_dict.keys():
        # Add entries for oxygen saturation shift and oxygen concentration derivation
        print(
            "        SHIFT    "
            + "{:7}".format(str(1.0))
            + "{:11}".format(
                metadata_dict["SHIFT_Oxygen_Time"]
                    .strftime(time_format)[0:-7]
                    .split(" ")[0]
            )
            + "{:9}".format(
                metadata_dict["SHIFT_Oxygen_Time"]
                    .strftime(time_format)[0:-7]
                    .split(" ")[1]
            )
            + "{:>9}".format(str(cast_shift_c["cast" + str(cast_number)].shape[0]))
            + "{:>10}".format(str(cast_shift_o["cast" + str(cast_number)].shape[0]))
        )
        print(
            "        DERIVE   "
            + "{:7}".format(str(1.0))
            + "{:11}".format(
                metadata_dict["DERIVE_OXYGEN_CONCENTRATION_Time"]
                    .strftime(time_format)[0:-7]
                    .split(" ")[0]
            )
            + "{:9}".format(
                metadata_dict["DERIVE_OXYGEN_CONCENTRATION_Time"]
                    .strftime(time_format)[0:-7]
                    .split(" ")[1]
            )
            + "{:>9}".format(str(cast_shift_o["cast" + str(cast_number)].shape[0]))
            + "{:>10}".format(str(cast_d_o_conc["cast" + str(cast_number)].shape[0]))
        )
    # applies to all cases
    print(
        "        DELETE   "
        + "{:7}".format(str(1.0))
        + "{:11}".format(
            metadata_dict["DELETE_PRESSURE_REVERSAL_Time"]
                .strftime(time_format)[0:-7]
                .split(" ")[0]
        )
        + "{:9}".format(
            metadata_dict["DELETE_PRESSURE_REVERSAL_Time"]
                .strftime(time_format)[0:-7]
                .split(" ")[1]
        )
        + "{:>9}".format(str(cast_shift_o["cast" + str(cast_number)].shape[0]))
        + "{:>10}".format(
            str(
                cast_wakeeffect["cast" + str(cast_number)].shape[0]
                - list(cast_wakeeffect["cast" + str(cast_number)].isna().sum())[0]
            )
        )
    )
    print(
        "        BINAVE   "
        + "{:7}".format(str(1.0))
        + "{:11}".format(
            metadata_dict["BINAVE_Time"].strftime(time_format)[0:-7].split(" ")[0]
        )
        + "{:9}".format(
            metadata_dict["BINAVE_Time"].strftime(time_format)[0:-7].split(" ")[1]
        )
        + "{:>9}".format(
            str(
                cast_wakeeffect["cast" + str(cast_number)].shape[0]
                - list(cast_wakeeffect["cast" + str(cast_number)].isna().sum())[0]
            )
        )
        + "{:>10}".format(str(cast_binned["cast" + str(cast_number)].shape[0]))
    )
    if "DROP_SELECT_VARS_Time" in metadata_dict.keys():
        print(
            "        DROP_SLCT"
            + "{:7}".format(str(1.0))
            + "{:11}".format(
                metadata_dict["DROP_SELECT_VARS_Time"]
                    .strftime(time_format)[0:-7]
                    .split(" ")[0]
            )
            + "{:9}".format(
                metadata_dict["DROP_SELECT_VARS_Time"]
                    .strftime(time_format)[0:-7]
                    .split(" ")[1]
            )
            + "{:>9}".format(str(cast_binned["cast" + str(cast_number)].shape[0]))
            + "{:>10}".format(
                str(
                    cast_dropvars["cast" + str(cast_number)].shape[0]
                    - list(cast_dropvars["cast" + str(cast_number)].isna().sum())[0]
                )
            )
        )
    print(
        "        EDIT     "
        + "{:7}".format(str(1.0))
        + "{:11}".format(
            metadata_dict["FINALEDIT_Time"].strftime(time_format)[0:-7].split(" ")[0]
        )
        + "{:9}".format(
            metadata_dict["FINALEDIT_Time"].strftime(time_format)[0:-7].split(" ")[1]
        )
        + "{:>9}".format(str(cast_dropvars["cast" + str(cast_number)].shape[0]))
        + "{:>10}".format(str(cast_final["cast" + str(cast_number)].shape[0]))
    )

    print("    $END")
    print(" $REMARKS")

    list_number = len(metadata_dict["Processing_history"].split("|"))
    for i in range(list_number):
        print("     " + metadata_dict["Processing_history"].split("|")[i])
    print("$END")
    print()
    return


def write_comments(
        processing_report_name: str, channel_names   # have_fluor: bool, have_oxy: bool,
) -> None:
    """Write comments section in the IOS header file
    inputs:
        - have_fluor: boolean flag, True if fluorescence channel is available
        - have_oxy: boolean flag, True if oxygen channels are available
        - processing_report_name: name of the RBR processing report for the selected cruise
    outputs:
        - Comments section added to open IOS header file, but nothing returned by
        the function
    """
    # cruise_ID = metadata_dict["Mission"]
    print("*COMMENTS")
    print("    " + "-" * 85)
    print()
    print("    Data Processing Notes:")
    print("    " + "-" * 22)
    print("       " + "No calibration sampling was available.")
    print()
    # print("       " + "For details on the processing see document: " + cruise_ID + "RBR_Processing_Report.doc.")
    print(
        "       "
        + "For details on the processing see document: "
        + processing_report_name
        + "."
    )
    # add name of processing report as an input parameter as some say "RBR" and have docx suffix
    print()
    print("    " + "-" * 85)

    # ----------------------
    # 7-line table header
    print_list = ["!--1--- --2--- ",
                  "!Pressu Depth  ",
                  "!re            ",
                  "!              ",
                  "!              ",
                  "!              ",
                  "!------ ------ "]

    print_dict = {"Temperature": ["---*---- ",
                                  "Temperat ",
                                  "ure      ",
                                  "         ",
                                  "         ",
                                  "         ",
                                  "-------- "],
                  "Salinity": ["---*---- ",
                               "Salinity ",
                               "         ",
                               "         ",
                               "         ",
                               "         ",
                               "-------- "],
                  "Fluorescence": ["---*--- ",
                                   "Fluores ",
                                   "cence:  ",
                                   "URU     ",
                                   "        ",
                                   "        ",
                                   "------- "],
                  "Oxygen": ["---*--- "
                             "Oxygen: ",
                             "Dissolv ",
                             "ed:     ",
                             "Saturat ",
                             "ion:RBR ",
                             "------- "],
                  "Oxygen_mL_L": ["---*--- ",
                                  "Oxygen: ",
                                  "Dissolv ",
                                  "ed: RBR ",
                                  "        ",
                                  "        ",
                                  "------- "],
                  "Oxygen_umol_kg": ["---*--- ",
                                     "Oxygen: ",
                                     "Dissolv ",
                                     "ed: RBR ",
                                     "        ",
                                     "        ",
                                     "------- "],
                  "Conductivity": ["----*---- ",
                                   "Conductiv ",
                                   "ity       ",
                                   "          ",
                                   "          ",
                                   "          ",
                                   "--------- "],
                  "Observation_counts": ["-*--",
                                         "Numb",
                                         "er_o",
                                         "~bin",
                                         "_rec",
                                         "ords",
                                         "----"]}

    channel_counter = 3
    # Skip Pressure and Depth
    for cn in channel_names[2:]:
        # Add columns to print list
        if cn == "Observation_counts":
            # account for single vs double digit channel numbers observation counts
            # in terms of maintaining constant channel width
            print_list[0] += print_dict[cn][0].replace("*", str(channel_counter)[:4])
        else:
            print_list[0] += print_dict[cn][0].replace("*", str(channel_counter))
        for i in range(1, len(print_list)):
            print_list[i] += print_dict[cn][i]
        # Update channel counter
        channel_counter += 1

    print_list.append("*END OF HEADER")
    # Now print the statements
    for line in print_list:
        print(line)
    # ----------------------

    # Make more flexible if channels other than oxygen or fluorescence were dropped
    # if have_fluor and have_oxy:
    #     print("!--1--- --2--- ---3---- ---4---- ---5--- ---6--- ---7--- ---8--- ----9---- -10-")
    #     print("!Pressu Depth  Temperat Salinity Fluores Oxygen: Oxygen: Oxygen: Conductiv Numb")
    #     print("!re            ure               cence:  Dissolv Dissolv Dissolv ity       er_o")
    #     print("!                                URU     ed:     ed: RBR ed: RBR           ~bin")
    #     print("!                                        Saturat                           _rec")
    #     print("!                                        ion:RBR                           ords")
    #     print("!------ ------ -------- -------- ------- ------- ------- ------- --------- ----")
    #     print("*END OF HEADER")
    # elif have_fluor:
    #     print("!--1--- --2--- ---3---- ---4---- ---5--- ----7---- -8--")
    #     print("!Pressu Depth  Temperat Salinity Fluores Conductiv Numb")
    #     print("!re            ure               cence:  ity       er_o")
    #     print("!                                URU               ~bin")
    #     print("!                                                  _rec")
    #     print("!                                                  ords")
    #     print("!------ ------ -------- -------- ------- --------- ----")
    #     print("*END OF HEADER")
    # elif have_oxy:
    #     print("!--1--- --2--- ---3---- ---4---- ---6--- ---7--- ---8--- ----7---- -8--")
    #     print("!Pressu Depth  Temperat Salinity Oxygen: Oxygen: Oxygen: Conductiv Numb")
    #     print("!re            ure               Dissolv Dissolv Dissolv ity       er_o")
    #     print("!                                ed:     ed: RBR ed: RBR           ~bin")
    #     print("!                                Saturat                           _rec")
    #     print("!                                ion:RBR                           ords")
    #     print("!------ ------ -------- -------- ------- ------- ------- --------- ----")
    #     print("*END OF HEADER")
    # else:
    #     print("!--1--- --2--- ---3---- ---4---- ----5---- -6--")
    #     print("!Pressu Depth  Temperat Salinity Conductiv Numb")
    #     print("!re            ure               ity       er_o")
    #     print("!                                          ~bin")
    #     print("!                                          _rec")
    #     print("!                                          ords")
    #     print("!------ ------ -------- -------- --------- ----")
    #     print("*END OF HEADER")

    return


def write_data(
        cast_data: dict, cast_number: int, channel_names
) -> None:  # , have_fluor: bool, have_oxy: bool, cast_d):
    """
    Write data to header file, taking into account if fluorescence and oxygen data are there
    inputs:
        - have_fluor: boolean flag, True if fluorescence channel is available
        - have_oxy: boolean flag, True if oxygen channels are available
        - cast_data: dictionary containing the processed data for the selected cast
        - cast_number: the event number for the cast
    outputs:
        - Data values printed to open IOS header file, but nothing returned by the function
    """
    # --------------
    channel_widths_all = {"Pressure": "{:>7}",
                          "Depth": "{:>6}",
                          "Temperature": "{:>8}",
                          "Salinity": "{:>8}",
                          "Fluorescence": "{:>7}",
                          "Oxygen": "{:>7}",
                          "Oxygen_mL_L": "{:>7}",
                          "Oxygen_umol_kg": "{:>7}",
                          "Conductivity": "{:>9}",
                          "Observation_counts": "{:>4}"}

    channel_widths_available = {}
    for cn in channel_names:
        channel_widths_available[cn] = channel_widths_all[cn]

    # for cn in channel_widths.keys():  RuntimeError: dictionary changed size during iteration
    #     if cn not in channel_names:
    #         _ = channel_widths.pop(cn)

    for i in range(len(cast_data["cast" + str(cast_number)])):
        print_line = ""
        for cn, wd in channel_widths_available.items():
            if cn != "Observation_counts":
                print_line += (
                    wd.format(cast_data["cast" + str(cast_number)].loc[i, cn])
                    + " "
                )
            else:
                print_line += (
                        wd.format(cast_data["cast" + str(cast_number)].loc[i, cn])
                )  # Omit space after the last column, which is Observation_counts
        print(print_line)
    # ---------------
    # if have_fluor and have_oxy:
    #     for i in range(len(cast_data["cast" + str(cast_number)])):
    #         # print(cast_data['cast' + str(cast_number)]['Pressure'][i] +
    #         # cast_data['cast' + str(cast_number)]['Depth'][i] + "  ")
    #         print(
    #             "{:>7}".format(cast_data["cast" + str(cast_number)].loc[i, "Pressure"])
    #             + " "
    #             + "{:>6}".format(cast_data["cast" + str(cast_number)].loc[i, "Depth"])
    #             + " "
    #             + "{:>8}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Temperature"]
    #             )
    #             + " "
    #             + "{:>8}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Salinity"]
    #             )
    #             + " "
    #             + "{:>7}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Fluorescence"]
    #             )
    #             + " "
    #             + "{:>7}".format(cast_data["cast" + str(cast_number)].loc[i, "Oxygen"])
    #             + " "
    #             + "{:>7}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Oxygen_mL_L"]
    #             )
    #             + " "
    #             + "{:>7}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Oxygen_umol_kg"]
    #             )
    #             + " "
    #             + "{:>9}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Conductivity"]
    #             )
    #             + " "
    #             + "{:>4}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Observation_counts"]
    #             )
    #             + " "
    #         )
    # elif have_fluor:
    #     for i in range(len(cast_data["cast" + str(cast_number)])):
    #         print(
    #             "{:>7}".format(cast_data["cast" + str(cast_number)].loc[i, "Pressure"])
    #             + " "
    #             + "{:>6}".format(cast_data["cast" + str(cast_number)].loc[i, "Depth"])
    #             + " "
    #             + "{:>8}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Temperature"]
    #             )
    #             + " "
    #             + "{:>8}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Salinity"]
    #             )
    #             + " "
    #             + "{:>7}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Fluorescence"]
    #             )
    #             + " "
    #             + "{:>9}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Conductivity"]
    #             )
    #             + " "
    #             + "{:>4}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Observation_counts"]
    #             )
    #             + " "
    #         )
    # elif have_oxy:
    #     for i in range(len(cast_data["cast" + str(cast_number)])):
    #         print(
    #             "{:>7}".format(cast_data["cast" + str(cast_number)].loc[i, "Pressure"])
    #             + " "
    #             + "{:>6}".format(cast_data["cast" + str(cast_number)].loc[i, "Depth"])
    #             + " "
    #             + "{:>8}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Temperature"]
    #             )
    #             + " "
    #             + "{:>8}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Salinity"]
    #             )
    #             + " "
    #             + "{:>7}".format(cast_data["cast" + str(cast_number)].loc[i, "Oxygen"])
    #             + " "
    #             + "{:>7}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Oxygen_mL_L"]
    #             )
    #             + " "
    #             + "{:>7}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Oxygen_umol_kg"]
    #             )
    #             + " "
    #             + "{:>9}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Conductivity"]
    #             )
    #             + " "
    #             + "{:>4}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Observation_counts"]
    #             )
    #             + " "
    #         )
    # else:
    #     for i in range(len(cast_data["cast" + str(cast_number)])):
    #         # print(cast_data['cast' + str(cast_number)]['Pressure'][i] +
    #         # cast_data['cast' + str(cast_number)]['Depth'][i] + "  ")
    #         print(
    #             "{:>7}".format(cast_data["cast" + str(cast_number)].loc[i, "Pressure"])
    #             + " "
    #             + "{:>6}".format(cast_data["cast" + str(cast_number)].loc[i, "Depth"])
    #             + " "
    #             + "{:>8}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Temperature"]
    #             )
    #             + " "
    #             + "{:>8}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Salinity"]
    #             )
    #             + " "
    #             + "{:>9}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Conductivity"]
    #             )
    #             + " "
    #             + "{:>4}".format(
    #                 cast_data["cast" + str(cast_number)].loc[i, "Observation_counts"]
    #             )
    #             + " "
    #         )

    return


def main_header(
        dest_dir: str,
        n_cast: int,
        metadata_dict: dict,
        cast: dict,
        cast_d: dict,
        cast_d_correct_t: dict,
        cast_d_calib: dict,
        cast_d_clip: dict,
        cast_d_filtered: dict,
        cast_d_shift_c: dict,
        cast_d_shift_o: dict,
        cast_d_o_conc: dict,
        cast_d_wakeeffect: dict,
        cast_d_binned: dict,
        cast_d_dropvars: dict,
        cast_d_final: dict,
        channel_names,
        processing_report_name: str,
) -> str:
    """
    Main function for creating an IOS header file containing final processed RBR CTD data
    inputs:
        - dest_dir: working directory that output files are saved to
        - n_cast: cast/event number
        - metadata_dict: dictionary containing metadata for the selected RBR cruise
        - cast: dictionary containing the original raw data from both the upcast and
        downcast
        - cast_d: dictionary containing the original raw data from the downcast only
        - cast_d_clip: dictionary containing the clipped downcast data
        - cast_d_filtered: dictionary containing the low-pass filtered downcast data
        - cast_d_shift_c: dictionary containing the downcast data that had conductivity
        shifted
        - cast_d_shift_o: dictionary containing the downcast data that had oxygen shifted
        - cast_d_wakeeffect: dictionary containing the downcast data that had pressure
        reversals deleted
        - cast_d_binned: dictionary containing the downcast data that was binned into
        1dbar vertical bins
        - cast_d_final: dictionary containing the final processed data, with conductivity
        in units of S/m
        - have_fluor: boolean flag, True if fluorescence channel is available
        - have_oxy: boolean flag, True if oxygen channels are available
    outputs:
        - absolute path of the output data file
    """
    path_slash_type = "/" if "/" in dest_dir else "\\"
    f_name = dest_dir.split(path_slash_type)[-2]
    f_output = f_name.split("_")[0] + "-" + f"{n_cast:04}" + ".CTD"
    new_dir = os.path.join(dest_dir, f"CTD{path_slash_type}")
    output = new_dir + f_output
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    # Start
    # datetime object containing current date and time
    now = datetime.now()

    # dd/mm/YY H:M:S
    dt_string = now.strftime("%Y/%m/%d %H:%M:%S.%f")[0:-4]

    IOS_string = "*IOS HEADER VERSION 2.0      2020/03/01 2020/04/15 PYTHON"

    orig_stdout = sys.stdout
    file_handle = open(output, "wt")
    try:
        sys.stdout = file_handle
        print("*" + dt_string)
        print(IOS_string)
        print()  # print("\n") pring("\n" * 40)
        write_file(n_cast, cast, cast_d_final, metadata_dict)
        write_admin(metadata_dict=metadata_dict)
        write_location(cast_number=n_cast, metadata_dict=metadata_dict)
        write_instrument(metadata_dict=metadata_dict)
        write_history(
            cast_d,
            cast_d_correct_t,
            cast_d_calib,
            cast_d_clip,
            cast_d_filtered,
            cast_d_shift_c,
            cast_d_shift_o,
            cast_d_o_conc,
            cast_d_wakeeffect,
            cast_d_binned,
            cast_d_dropvars,
            cast_d_final,
            cast_number=n_cast,
            metadata_dict=metadata_dict,
        )
        write_comments(
            processing_report_name, channel_names
        )  # , metadata_dict=meta_data, cast_d=cast_d) have_fluor, have_oxy,
        write_data(
            cast_d_final, cast_number=n_cast, channel_names=channel_names
        )  # , cast_d=cast_d) have_fluor, have_oxy,
        sys.stdout.flush()  # Recommended by Tom
    finally:
        sys.stdout = orig_stdout

    return os.path.abspath(output)


def get_started(dest_dir: str):
    """Start by opening the RSK files, find out how many channels and profiles there are

    Compare this to the cast list given by chief scientist

    prep metadata .csv

    have the header-merge.csv ready
    """

    files = os.listdir(dest_dir)  # list all the files in dest_dir
    files = list(filter(lambda f: f.endswith(".rsk"), files))  # keep the rsk files only
    n_files = len(files)  # get the number of files
    print(n_files)

    for k in range(n_files):
        filename = str(dest_dir) + str(files[k])  # full path and name of .rsk file
        print(filename)
        rsk = pyrsktools.RSK(filename)  # load up an RSK
        rsk.open()
        rsk.readdata()

        # check the number of profiles
        try:
            downcastIndices = rsk.getprofilesindices(direction="down")
        except AttributeError:
            rsk.computeprofiles()  # This should fix the problem
            downcastIndices = rsk.getprofilesindices(direction="down")

        n_profiles = len(downcastIndices)  # get the number of profiles recorded

        # last_profile = n_profiles[:-1] # get the last profile only
        # last_profile = list(rsk.profiles())[:-1]
        print("Number of profiles:", n_profiles)
        print("Samples:\n", rsk.samples)
        print("Channel items:\n", list(rsk.channels.items()))

    return


def first_step(
    dest_dir,
    year: str,
    cruise_number: str,
    data_file_type: str,
    skipcasts,
    rsk_time1=None,
    rsk_time2=None,
    left_lon=None,
    right_lon=None,
    bot_lat=None,
    top_lat=None,
) -> None:
    """
    Choose how to export the csv files from the rsk files, in preparation for processing
    Plot cruise, plot pre-processing plots, determine need for zero-order holds correction

     inputs:
        - dest_dir, year, cruise_number
        - data_file_type: 'rsk' for *.rsk files or 'excel' for *.xlsx type files
        - skipcasts: number of casts to skip over in each data file (rsk or excel)
            when writing data to output format. Input format as either an integer
            or as a list-like object with one integer per excel file
            representing the number of initial casts to skip in each excel file.
        - data_file_type: "rsk" for single or multiple rsk files, or "excel"
        (use on .xlsx files exported from Ruskin)
        - left_lon, right_lon, bot_lat, top_lat: map extent for plotting cast locations
     Outputs:
        - None in terms of python objects, but data files and plots are saved

     ***IMPORTANT: DON"T FORGET TO CHANGE THE SHEETNAME TO Profile_annotation or
     openpyxl won't read it.
    """

    if data_file_type == "rsk":
        READ_RSK(dest_dir, year, cruise_number, skipcasts, rsk_time1, rsk_time2)
    elif data_file_type == "excel":
        # input file =  # not needed, filtered to keep only .xls
        READ_EXCELrsk(dest_dir, year, cruise_number, skipcasts)

    # Merge all data from a cruise into one csv file
    MERGE_FILES(dest_dir, year, cruise_number)
    print("files merged")
    ADD_6LINEHEADER_2(dest_dir, year, cruise_number)

    # Make preliminary plots
    plot_track_location(
        dest_dir, year, cruise_number, left_lon, right_lon, bot_lat, top_lat
    )
    first_plots(year, cruise_number, dest_dir, input_ext="_CTD_DATA-6linehdr.csv")
    PLOT_PRESSURE_DIFF(
        dest_dir, year, cruise_number, input_ext="_CTD_DATA-6linehdr.csv"
    )
    return


def second_step(
    dest_dir: str,
    year: str,
    cruise_number: str,
    processing_report_name: str,
    rsk_file,
    rsk_time1=None,
    rsk_time2=None,
    pd_correction_value=0,
    window_width=6,  # sample_rate=8, time_constant=1 / 8,
    filter_type=1,
    shift_recs_conductivity=2,
    shift_recs_oxygen=-11,
    verbose: bool = False,
):
    """
    Run the processing steps for RBR CTD data
    inputs:
        - dest_dir
        - year
        - cruise_number
        - correction_value: value to correct pressure and depth data by using CALIB
        - input_ext: '_CTD_DATA-6linehdr.csv' or '_CTD_DATA-6linehdr_corr_hold.csv'
    outputs:
        - IOS header-format files containing processed RBR CTD data are saved to
        dest_dir/CTD/. This function itself does not return anything
    """

    # Create metadata dict
    metadata_dict = CREATE_META_DICT(
        dest_dir=dest_dir,
        rsk_file=rsk_file,
        year=year,
        cruise_number=cruise_number,
        rsk_time1=rsk_time1,
        rsk_time2=rsk_time2,
    )

    # # initialize casts for if statement
    # cast, cast_d, cast_u = 0, 0, 0
    # cast_pc, cast_d_pc, cast_u_pc = 0, 0, 0

    if verbose:
        print("Checking need for zero-order hold correction...")
    # Check pressure channel for zero order holds
    zoh = check_for_zoh(
        dest_dir, year, cruise_number, float(metadata_dict["Sampling_Interval"])
    )

    if zoh:
        # Correct the zero-order-holds
        CORRECT_HOLD(dest_dir, year, cruise_number, metadata_dict)

        input_ext = "_CTD_DATA-6linehdr_corr_hold.csv"

        if verbose:
            print("Using zero-order holds corrected variables")

        # check the plot then save it as Fig_3
        PLOT_PRESSURE_DIFF(dest_dir, year, cruise_number, input_ext)
    else:
        input_ext = "_CTD_DATA-6linehdr.csv"

        if verbose:
            print("using original variables")

    cast, cast_d, cast_u = CREATE_CAST_VARIABLES(
        year, cruise_number, dest_dir, input_ext
    )

    for cast_i in cast_d.keys():
        have_oxy = True if "Oxygen" in cast_d[cast_i].columns else False
        have_fluor = True if "Fluorescence" in cast_d[cast_i].columns else False
        if verbose:
            print(f"have_oxy: {have_oxy}, have_fluor: {have_fluor}")
        break

    # Calibrate pressure and depth
    cast_pc, cast_d_pc, cast_u_pc = CALIB(
        cast, cast_d, cast_u, metadata_dict, zoh, pd_correction_value
    )  # 0 if no neg pressures
    if verbose:
        print(
            "The following correction value has been applied to Pressure and Depth:",
            pd_correction_value,
            sep="\n",
        )

    # Commented this out because first plots are made earlier and in the past haven't been
    # remade after the zero hold correction. Plus there are two possible input_ext now.
    # first_plots(year, cruise_number, dest_dir, input_ext='_CTD_DATA-6linehdr_corr_hold.csv')
    # print('finished plotting first plots')

    # clip the casts
    # cast_d_clip = CLIP_DOWNCAST(cast_d_pc, metadata_dict, limit_drop=0.02)
    # cast_u_clip = CLIP_UPCAST(cast_u_pc, metadata_dict, limit_rise=-0.02)
    cast_d_clip = CLIP_CAST(
        cast_d_pc, metadata_dict, limit_pressure_change=0.02, cast_direction="down"
    )
    cast_u_clip = CLIP_CAST(
        cast_u_pc, metadata_dict, limit_pressure_change=-0.02, cast_direction="up"
    )

    plot_clip(cast_d_clip, cast_d_pc, dest_dir)

    if verbose:
        print("Casts clipped")

    # Apply a low-pass filter
    # metadata_dict['Sampling_Interval']: time in seconds between records
    # (0.125-0.167s for an 8-6Hz instrument)
    sample_rate = int(np.round(1 / float(metadata_dict["Sampling_Interval"])))
    cast_d_filtered, cast_u_filtered = FILTER(
        cast_d_clip,
        cast_u_clip,
        metadata_dict,
        have_fluor,
        window_width,
        sample_rate=sample_rate,
        time_constant=float(metadata_dict["Sampling_Interval"]),
        filter_type=filter_type,
    )  # n = 5 should be good.

    # Plot the filtered data
    plot_filter(
        cast_d_filtered, cast_u_filtered, cast_d_clip, cast_u_clip, dest_dir, have_fluor
    )

    if verbose:
        print(
            f"Casts filtered, assuming a sample rate of {sample_rate} records per second"
        )

    # Shift the conductivity channel and recalculate salinity after
    # Default conductivity shift is a delay by 2 scans
    cast_d_shift_c, cast_u_shift_c = SHIFT_CONDUCTIVITY(
        cast_d_filtered,
        cast_u_filtered,
        metadata_dict=metadata_dict,
        shifted_scan_number=shift_recs_conductivity,
    )

    plot_shift_c(
        cast_d_shift_c, cast_u_shift_c, cast_d_filtered, cast_u_filtered, dest_dir
    )

    if verbose:
        print(f"Conductivity shifted {shift_recs_conductivity} scans")

    # Shift oxygen channel if available, then calculate oxygen concentration
    if have_oxy:
        cast_d_shift_o, cast_u_shift_o = SHIFT_OXYGEN(
            cast_d_shift_c,
            cast_u_shift_c,
            metadata_dict=metadata_dict,
            shifted_scan_number=shift_recs_oxygen,
        )

        plot_shift_o(
            cast_d_shift_o, cast_u_shift_o, cast_d_shift_c, cast_u_shift_c, dest_dir
        )

        if verbose:
            print(f"Oxygen shifted {shift_recs_oxygen} scans")

        # Add cast variables of oxygen conc. in ml/l and umol/kg
        cast_d_o_conc, cast_u_o_conc = DERIVE_OXYGEN_CONCENTRATION(
            cast_d_shift_o, cast_u_shift_o, metadata_dict
        )

        if verbose:
            print("Oxygen concentration derived from oxygen saturation")
    else:
        # Rename the variables
        cast_d_shift_o, cast_u_shift_o = cast_d_shift_c, cast_u_shift_c
        cast_d_o_conc, cast_u_o_conc = cast_d_shift_c, cast_u_shift_c

    # Delete the upcast and all other pressure change reversals
    cast_d_wakeeffect, cast_u_wakeeffect = DELETE_PRESSURE_REVERSAL(
        cast_d_o_conc, cast_u_o_conc, metadata_dict=metadata_dict
    )

    if verbose:
        print("Deleted pressure change reversals")

    # Plot before and after comparisons of the delete step
    plot_delete(cast_d_wakeeffect, cast_d_o_conc, dest_dir)

    # Average the data into 1-dbar bins
    cast_d_binned, cast_u_binned = BINAVE(
        cast_d_wakeeffect, cast_u_wakeeffect, metadata_dict=metadata_dict
    )

    if verbose:
        print("Records averaged into equal-width pressure bins")

    # Final edits: change conductivity units
    cast_d_final = FINAL_EDIT(cast_d_binned, have_oxy, have_fluor, metadata_dict)

    if verbose:
        print("Final edit completed")

    # Final plots after bin averaging
    # print(cast_d_final)
    plot_processed(cast_d_final, dest_dir)  # cast_d_o_conc, cast_d, cast,

    # # Produce the IOS Shell header files containing the final data
    # main_header(dest_dir, n_cast,
    #             metadata_dict, cast, cast_d, cast_d_clip, cast_d_filtered,
    #             cast_d_shift_c, cast_d_shift_o,
    #             cast_d_wakeeffect, cast_d_binned, cast_d_final, have_oxy)

    # Retrieve the cast numbers from the data dictionary
    # The cast numbers may not start at 1 and monotonously increase by 1
    cast_numbers = list(map(lambda x: int(x.split("cast")[-1]), cast_d_final.keys()))
    for i in cast_numbers:
        # # Update have_fluor and have_oxy flags, since it could be different for different casts
        # have_fluor = True if "Fluorescence" in cast_d_final[f"cast{i}"].columns else False
        # have_oxy = True if "Oxygen" in cast_d_final[f"cast{i}"].columns else False
        channel_names = cast_d_final[f"cast{i}"].columns
        # Call the main header creation function
        main_header(
            dest_dir,
            i,
            metadata_dict,
            cast,
            cast_d,
            cast_d_correct_t,
            cast_d_pc,
            cast_d_clip,
            cast_d_filtered,
            cast_d_shift_c,
            cast_d_shift_o,
            cast_d_o_conc,
            cast_d_wakeeffect,
            cast_d_binned,
            cast_d_dropvars,
            cast_d_final,
            channel_names,
            processing_report_name,
        )

    if verbose:
        print("Header files produced")

    return


def PROCESS_RBR(
    dest_dir,
    year: str,
    cruise_number: str,
    processing_report_name: str,
    rsk_file: str,
    data_file_type: str,
    window_width: int,
    skipcasts,
    filter_type: int = 1,
    # sample_rate: int, time_constant: float,
    shift_recs_conductivity: int = -2,
    shift_recs_oxygen=None,
    pd_correction_value=0,
    rsk_time1=None,
    rsk_time2=None,
    left_lon=None,
    right_lon=None,
    bot_lat=None,
    top_lat=None,
    verbose=False,
):
    """
    Main function to run the suite of processing functions on raw RBR data
    inputs:
        - dest_dir
        - year
        - cruise_number: 3-character string, e.g., '015', '101'
        - processing_report_name
        - rsk_file: relative file name of an rsk file in the dest_dir
        - skipcasts: number of casts to skip over in each data file (rsk or excel)
            when writing data to output format. Input format as either an integer
            or as a list-like object with one integer per excel file
            representing the number of initial casts to skip in each excel file.
        - data_file_type: 'rsk' for ruskin files or "excel"
        (use on .xlsx files exported from Ruskin)
        - nprof_per_rsk: int or list-like if multiple input rsk data files
        - left_lon, right_lon, bot_lat, top_lat: map extent for plotting cast locations
    """

    first_step(
        dest_dir,
        year,
        cruise_number,
        data_file_type,
        skipcasts,
        rsk_time1,
        rsk_time2,
        left_lon,
        right_lon,
        bot_lat,
        top_lat,
    )

    second_step(
        dest_dir,
        year,
        cruise_number,
        processing_report_name,
        rsk_file,
        rsk_time1,
        rsk_time2,
        pd_correction_value,
        window_width,  # sample_rate, time_constant,
        filter_type,
        shift_recs_conductivity,
        shift_recs_oxygen,
        verbose,
    )
    return


def test_process():
    # test_year = '2022'
    # test_cruise_num = '025'
    # test_file = '201172_20220925_1059.rsk'
    # # num_profiles = 6
    # skipcasts = 0

    # test_year = '2023'
    # test_cruise_num = '015'
    # test_file = '208765_20230121_2113_newDOcoefficients.rsk'
    # # num_profiles = 44
    # skipcasts = 54  # First 53 casts from an earlier cruise + 1 cast with error

    test_year = "2022"  # Issue with rsk.open() on files from this cruise
    test_cruise_num = "031"
    test_file = "204848_20220304_2322-Haro Strait.rsk"
    num_profiles = 1
    skipcasts = [1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0]

    test_dir = (
        "C:\\Users\\HourstonH\\Documents\\ctd_processing\\RBR\\"
        "python_testing\\{}-{}\\".format(test_year, test_cruise_num)
    )
    # test_event_start = 1
    processing_report_name = f"{test_year}-{test_cruise_num}_RBR_Processing_Report.docx"

    # EXPORT_MULTIFILES(test_dir, num_profiles, test_event_start, 'ALL')
    # MERGE_FILES(test_dir, test_year, test_cruise_num, event_from='header-merge')

    # first_step(test_dir, test_year, test_cruise_num, test_event_start, 'ALL', 'rsk',
    #            num_profiles)

    # second_step(test_dir, test_year, test_cruise_num, processing_report_name,
    #             test_file, window_width=3, filter_type=1, verbose=True)

    PROCESS_RBR(
        test_dir,
        test_year,
        test_cruise_num,
        processing_report_name=processing_report_name,
        rsk_file=test_file,
        data_file_type="rsk",
        skipcasts=skipcasts,
        window_width=3,
        shift_recs_conductivity=2,  # sample_rate=8, time_constant=1 / 8,
        shift_recs_oxygen=-11,
        verbose=True,
    )

    return
