#run e.g. via PyCharm, python3.11 generate_averages_and_kymograph_data.py
import pandas as pd
import numpy as np
import os
import pickle

import scipy
#import ffmpeg
import scipy.interpolate
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

# set which simulation run to analyze
name = "simpleformat_27_1"
# name="simpleformat_26_1_M1"


# name="simpleformat_28_1"
# name="simpleformat_29_3"
# name="simpleformals
# t_29_4_regenerated"
# name="simpleformat_30_1_M0_GaussPoints"
# name="simpleformat_30_1_M1"
# name="simpleformat_30_1_M1"

# select which data to run
#"simpleformat_26_1_M1" #with mass movement, varying inflow. 4 different D, 72h
#"simpleformat_26_1_M0" #without movement, varying inflow. 4 different D, 72h
#"simpleformat_27_1" #with mass movement, constant inflow and constant nutrient concentration in inflow. 4 different D, 72h

#simpleformat_26_1_SCFA #with mass movement, varying inflow. 4 different D, 72h - simulates SCFA secretion and uptake

#"simpleformat_28_1" #no mass movement, constant inflow and constant nutrient concentration in inflow. 4 different D, 72h

#"simpleformat_29_3" #WITH mass movement, varying inflow with varying nutrient concentrations. NO CECUM. 4 different D, 72h

#"simpleformat_29_4" #NO mass movement, varying inflow with varying nutrient concentrations. NO CECUM. 4 different D, 72h

#"simpleformat_30_1_M0" no mass movement, constant inflow constant nutrient concentration in inflow. NO CECUM. 4 different D, 72h
#"simpleformat_30_1_M1" mass movement constant inflow constant nutrient concentration in inflow. NO CECUM. 4 different D, 72h


#"30_1_M0","30_1_M1","27_1","28_1","29_3","29_4",


namelist = ["26_5_SCFA_Normal"]
#namelist = ["26_3_SCFA_FinerRes"]
includelambda=False #read in different lambda
#namelist = ["29_4"]
#already run 26_1_M0, 26_1_M1
for name in namelist:

    varnames=["c2","c3","u","w"]
    varnames_long=["nutrients (mM)","bacteria (OD)","$v_{r}$","$v_{z}$"]
    if "SCFA" in name:

        varnames = ["c2", "c3", "u", "w", "c4"]
        varnames_long = ["nutrients (mM)", "bacteria (OD)", "$v_{r}$", "$v_{z}$","SCFA (mM)"]

    #name="simpleformat_"+\
    #name

    foldername="data_simpleformat" #where data is stored


    #which variables to save as kymographs

    #set how many points to store along z direction
    zpoints=100 #how many points alon z

    #start preparing Kymographs
    data=pd.read_csv(os.path.join(foldername,name+".txt"))
    #display(data)

    #get unique values:
    t_unique=data["t"].unique()
    t_unique.sort()
    r_unique=data["r"].unique()
    r_unique.sort()
    z_unique=data["z"].unique()
    z_unique.sort()
    D_unique=data["D"].unique()
    D_unique.sort()
    M1_unique=data["M1"].unique()
    M1_unique.sort()



    #display unique values
    print("t values")
    print(t_unique)
    print("r values")
    print(len(r_unique))
    print("z values")
    print(len(z_unique))
    print("D values")
    print(D_unique)
    print("M1 values")
    print(M1_unique)

    print(max(z_unique))
    D_unique = sorted(data["D"].unique())
    print("D unique")
    print(D_unique)

    if includelambda:
        if len(D_unique)>1:
            error_several_diffusionvaluesaswell
        D_unique = data["Lambda"].unique()
        D_unique.sort()
        print("Lambda")
        print(D_unique)



    #set settings for z resolution when averaging
    x=data["z"].max()
    zrange=np.linspace(data["z"].min(),data["z"].max(),zpoints)
    delta_z=zrange[-1]/float(zpoints)

    #calculate index to plot over last 24 hours
    indext_last24=list(t_unique).index(48)

    #prepare dataframe to store average across 24h
    #dcolnames=[]
    #for Dindex in range(len(D_unique)):
    #    for var in varnames:
    #        dcolnames.append("D"+str(Dindex+1)+var)
    #data_avtime=pd.DataFrame(columns=["z"]+dcolnames)
    #data_avtime["z"]=zrange



    #go through every diffusion coefficient
    for Dindex in range(0,len(D_unique)):

        print("D index: "+str(Dindex))
        #Dindex=0
        if includelambda:
            select = data.loc[(data["Lambda"] == D_unique[Dindex])]
        else:
            select = data.loc[(data["D"] == D_unique[Dindex])]

        #prepare foldername to store plots
        foldername=name+"_D"+str(Dindex)
        if not os.path.exists(foldername):
           os.makedirs(foldername)

        #find for all variables the min and max values (to manaint same ranges across plots)
        var_list_min_global=[]
        var_list_max_global=[]
        for var in varnames:
            var_list=select[var].values
            var_list_min_global.append(var_list.min())
            var_list_max_global.append(var_list.max())

        #if you want to set plotting range by hand for full profile plots
        var_list_min_manual=[0,-10,-0.00001,-0.00002]
        var_list_max_manual=[30,200,0.00001,0.0002]

        ##############################
        #prepare data for kymograph
        #############################

        ### discussion: how to add kymograph???
        #one possibility: calculate average of values, weighting with radial distance

        dataout_lists={}  #for radial average
        dataout_lists["time"]=t_unique
        dataout_lists["z"] = zrange

        dataout_lists["D_values"]=D_unique
        dataout_lists["M1_values"]=M1_unique

        for iV in range(0,len(varnames)):
            print("prep kymograph for "+varnames[iV])

            #reserve arrays to generate kymographs
            kymo_radav=np.zeros([t_unique.shape[0],zpoints])
            kymo_centerline=np.zeros([t_unique.shape[0],zpoints])
            maxvalue=np.zeros([t_unique.shape[0],zpoints])
            minvalue = np.zeros([t_unique.shape[0], zpoints])
            systemradius = np.zeros([t_unique.shape[0]])

            av_value = np.zeros([t_unique.shape[0]])
            #sum_radial=np.zeros([t_unique.shape[0],zpoints])
            #sum_total=np.zeros([t_unique.shape[0]])

            #go through every time-pioint
            tcount=-1
            for t_cur in t_unique:
                select2=select.loc[ (select["t"]==t_cur)]

                #print(tcount)
                tcount=tcount+1

                av_value_cur=0
                for z_count in range(0,zpoints):

                    select3=select2.loc[(select2["z"]>z_count*delta_z) & (select2["z"]<=(z_count+1)*delta_z) ]
                    if z_count==0:
                        systemradius[tcount]=select3["r"].max()

                    #average over cross-section (for kymograph plots)
                    av=0
                    rad=0

                    # kymograph - only points at centerline
                    av2 = 0
                    rad2 = 0

                    #get minimal values
                    maxvalue[tcount,z_count]=select3[varnames[iV]].max()
                    minv = select3[varnames[iV]].min()
                    if minv<0:
                        minv=0
                    minvalue[tcount, z_count] = minv


                    for index, row in select3.iterrows():
                        if row["r"]>=0:
                            av=av+row[varnames[iV]]*row["r"]
                            rad=rad+row["r"]

                            if row["r"] < 0.001:  # radial average only of points close to center (average not weighted by radius)
                                av2 = av2 + row[varnames[iV]]
                                rad2 = rad2 + 1


                    av=av/rad #normalize by radial values

                    av_value_cur=av_value_cur+av




                    kymo_radav[tcount,z_count]=av  #add value to array
                    if rad2 == 0: #division by zero should not occur
                        currs=select3["r"].unique()
                        currs.sort()
                        print(currs)
                        error0
                    av2 = av2 / rad2  # normalize by radial values
                    kymo_centerline[tcount, z_count] = av2  # add value to array
                    ##total amount across radial section
                    #sum = 0
                    #for index, row in select3.iterrows():
                    #    if row["r"] > 0:
                    #        sum = av + np.pi*row[varnames[iV]] * row["r"]
                    #sum_radial[tcount, z_count] = sum  # add value to array
                    #sum_total[tcount] = sum_total[tcount] + sum*delta_z







                    #elif 3>2: #
                    #    tobeimplemented
                    #    data.groupby(pd.cut(df["B"], np.arange(0, 1.0+0.155, 0.155))).sum()
                    #else:
                    #    select=data.loc[ (data["t"]==t_cur) & (data["D"]==D_unique[Dindex]) & (data["z"]>z_count*delta_z) & (data["z"]<=(z_count+1)*delta_z) ]
                av_value_cur = av_value_cur / (1. * zpoints)
                av_value[tcount] = av_value_cur
            #take average of last 24 hours

            #data_avtime["D"+str(Dindex+1)+var]=kymograph_lists[-1][indext_last24:,:].mean(axis=0)
            #data_avtime["D"+str(Dindex+1)+var+"_centerline"]=kymograph_lists2[-1][indext_last24:,:].mean(axis=0)

            dataout_lists["kymo_radav_"+varnames[iV]]=kymo_radav
            dataout_lists["average_" + varnames[iV]] = av_value
            dataout_lists["kymo_centerline_" + varnames[iV]] = kymo_centerline
            dataout_lists["minradial_" + varnames[iV]] = minvalue
            dataout_lists["maxradial_" + varnames[iV]] = maxvalue
            dataout_lists["systemradius"] = systemradius


        with open(os.path.join("data_kymographs",foldername+'calculateddata.pickle'), 'wb') as handle:
            pickle.dump(dataout_lists, handle, protocol=pickle.HIGHEST_PROTOCOL)
        #np.save(os.path.join("data_kymographs",foldername+'kymographdata.npy'), dataout_lists, allow_pickle=True)
        print("saved sucessfully"+foldername)
    #data_avtime.to_csv(os.path.join("data_kymographs",name+"_24haverage.csv"))