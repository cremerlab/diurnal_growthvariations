import pandas as pd
import numpy as np
import os
import pickle
import scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from mpl_toolkits.axes_grid1.colorbar import colorbar
#import ffmpeg
import scipy.interpolate
import matplotlib.pyplot as plt
import matplotlib.colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import scipy.integrate as spi
#from scipy.integrate import odeint #this is the module to solve ODEs
#%matplotlib inline

#mode="linear" #"linear, log"
#"simpleformat_26_M0"


foldername_data="data_simpleformat"
foldername_outbase="videos"
foldername_kymograph="data_kymographs"
if not os.path.exists(foldername_outbase):
               os.makedirs(foldername_outbase)


# select which data to run
#"simpleformat_26_1_M1" #with mass movement, varying inflow. 4 different D, 72h
#"simpleformat_26_1_M0" #without movement, varying inflow. 4 different D, 72h
#"simpleformat_29_4" #NO mass movement, varying inflow with varying nutrient concentrations. NO CECUM. 4 different D, 72h

#"simpleformat_30_1_M0" no mass movement, constant inflow constant nutrient concentration in inflow. NO CECUM. 4 different D, 72h
#"simpleformat_30_1_M1" mass movement constant inflow constant nutrient concentration in inflow. NO CECUM. 4 different D, 72h

#run it...
#simpleformat_26_1_M1

#run ["26_1_M1"]



####################
##set settings
####################
#Done
namelist = ["26_1_M0"]
no_massflow_list=[True] #set if massflow is not included (this is important to set a different aspect ratio of plots when radius remains smaller
namelist = ["30_1_M0"]
no_massflow_list=[True]

namelist = ["26_1_M1"]
no_massflow_list=[False]

namelist = ["30_1_M1"]
no_massflow_list=[False]

#"simpleformat_27_1" #with mass movement, constant inflow and constant nutrient concentration in inflow. 4 different D, 72h
#simpleformat_26_1_SCFA #with mass movement, varying inflow. 4 different D, 72h - simulates SCFA secretion and uptake
#"simpleformat_28_1" #no mass movement, constant inflow and constant nutrient concentration in inflow. 4 different D, 72h
#"simpleformat_29_3" #WITH mass movement, varying inflow with varying nutrient concentrations. NO CECUM. 4 different D, 72h





namelist = ["29_4"]
no_massflow_list=[True]

namelist = ["26_4_SCFA"]
no_massflow_list=[False]
#options:
boundary_option="fixed"# "minmax"  #which ranges to plot (current "minmax", or fixed). minmax might lead to areas without data in LogNorm plot.
repeatrun=False #If false only plots are generated for which no file exist.


modelist=["log"] #lin or log
Dindexrange=[0] #2,3,1,
#Dindexrange=[2]
deltaT=48


conversionv=3600 #set unit of velocity (m/h instead of m/s)
###################



nc=-1
for name in namelist:
    nc=nc+1
    no_massflow=no_massflow_list[nc]
    if "SCFA" in name:
        varnames = ["c2", "c3","c4", "u", "w"]
        varnames_long = ["nutrients $n$ (mM)", 'bacteria $ \\rho $ (OD)', "SCFA (mmol)","velocity $v_{r}$ (m/h)", "velocity $v_{z}$ (m/h)"]
    else:
        varnames = ["c2", "c3", "u", "w"]
        varnames_long = ["nutrients $n$ (mM)", 'bacterial density $ \\rho$ (OD)', "velocity $v_{r}$ (m/h)", "velocity $v_{z}$ (m/h)"]
    print(varnames)
    data = pd.read_csv(os.path.join(foldername_data, name + ".txt"))
    #select only times larger than certain value
    data = data.loc[data["t"] >= deltaT]

    shortname = name[-5:]
    for mode in modelist:
        #read in file (simple format)
        




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


        for Dindex in Dindexrange: #add D again
            #decide which parameter to use (u)
            select=data.loc[(data["D"]==D_unique[Dindex])]
            runname=shortname+"D="+str(D_unique[Dindex])
            #which variables to plot

            #prepare foldername to store plots
            foldername=os.path.join(foldername_outbase,name+"_D"+str(Dindex))
            kymoname=name+"_D"+str(Dindex)+'kymographdata.npy'
            if mode=="log":
                foldernameout=foldername+"log"
            elif mode=="linear":
                foldernameout=foldername
            else:
                errormodenotfound
            if not os.path.exists(foldernameout):
               os.makedirs(foldernameout)

            #find for all variables the min and max values (to manaint same ranges across plots)
            var_list_min_global=[]
            var_list_max_global=[]
            for var in varnames:
                var_list=select[var].values
                var_list_min_global.append(np.nanmin(var_list))
                var_list_max_global.append(np.nanmax(var_list))




            #if you want to set plotting range by hand for full profile plots
            var_list_min_manual=[0,-10,-0.00001,-0.00002]
            var_list_max_manual=[30,200,0.00001,0.0002]

            # set range manually
            #vmin = 0
            #vmaxlistmanual = [180, 15, np.power(10., -5.), np.power(10., -5.)]  # for nutrients, OD, etc
            #vmax = vmaxlistmanual[iV]

            #load kymograph data

            picklename = name + "_D" + str(Dindex) + 'calculateddata.pickle'
            with open(os.path.join("data_kymographs", picklename), 'rb') as handle:
                data_calculated = pickle.load(handle)

            #look at r min and max

            max_r=[]
            tcount=0
            for t_cur in t_unique: #t_unique[-1]

                    tcount=tcount+1
                    select2=select.loc[ (select["t"]==t_cur)]
                    r=select2["r"].values #get coordinates as array
                    z=select2["z"].values #get coordinates as array
                    t=select2["t"]
                    max_r.append(r.max())

            plt.plot(t_unique,100*np.array(max_r),ls='-',marker='o')
            plt.xlabel("time (h)")
            plt.ylabel("radius (cm)")
            plt.ylim(0,3)
            #plt.show()


            tcount=0
            for t_cur in t_unique: #t_unique[-1]
                tcount = tcount + 1
                print(t_cur)
                timefilename = str(tcount).zfill(3) + ".png"
                if os.path.exists(os.path.join(foldernameout, timefilename))==False or repeatrun:

                        #plot one certain timepoint and parameter combination
                        select2=select.loc[ (select["t"]==t_cur)]
                        #display(select.shape)
                        #display(select)

                        ##############################
                        #interpolote for 2d plot: Set up a regular grid of interpolation points
                        #############################

                        r=select2["r"].values #get coordinates as array
                        z=select2["z"].values #get coordinates as array
                        t=select2["t"]

                        zi, ri = np.linspace(z.min(), z.max(), 100),np.linspace(r.min(), r.max(), 100) #prepare grid for interpolation
                        zi, ri = np.meshgrid(zi, ri) #prepare grid for interpolation

                        #prepare lists to store results
                        var_list_interpolated=[] #interpolated grid
                        var_list_interpolated2=[] #to storry mirrored version
                        var_list_min=[] #min and max values (to adjust plotting range)
                        var_list_max=[] #min and max values (to adjust plotting range
                        var_list=[]

                        for var in varnames:
                            var_list.append(select2[var].values)
                            var_list_min.append(var_list[-1].min())
                            var_list_max.append(var_list[-1].max())
                            #interpolate. Use resacling to get it to work, see https://stackoverflow.com/questions/17577587/matplotlib-2d-graph-with-interpolation
                            var_rescaled = (var_list[-1] - var_list[-1].min()) / var_list[-1].ptp()
                            rbf = scipy.interpolate.Rbf(z, r,var_rescaled , function='linear')
                            var_rescaled=rbf(zi,ri)
                            var_list_interpolated.append(var_list[-1].ptp()*var_rescaled + var_list[-1].min())
                            # mirror data along r=0:
                            mirroredcur = np.append(var_list_interpolated[-1][::-1, :], var_list_interpolated[-1], axis=0)
                            # flit x and y axis
                            var_list_interpolated2.append(np.transpose(mirroredcur))  #
                            # print(var_list_min_global)
                            # print(valr_list_max_global)

                        #plot everything in one plot

                        fig, ax = plt.subplots(len(varnames),2,figsize=(8,len(varnames)*5)) # its for the plots: 1 row and 3 panels in 1 row (figure size in cm)
                        fig.subplots_adjust(wspace=-.35)
                            #wighsize width and height
                        divider=[]
                        cax=[]

                        #plot for every variable
                        for iV in range(0,len(varnames)):
                            #decide which range for z values to use
                            if 3>4: #local min and max
                                vmin=var_list_min[iV]
                                vmax=var_list_max[iV]
                            elif 3>2: #global min and max
                                vmin=var_list_min_global[iV]
                                vmax=var_list_max_global[iV]
                            else: #manual
                                vmin=var_list_min_manual[iV]
                                vmax=var_list_max_manual[iV]

                            #choose colormap
                            #chose colormap here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
                            if iV<2:
                                cmap="viridis"
                            elif iV==2:
                                cmap="plasma"
                            else:
                                cmap="plasma_r"
                                #cmap="cividis"

                            ###################
                            #plot kymographs

                            extent = [t_unique.min()-deltaT, t_unique.max()-deltaT, z.min() * 100, z.max() * 100]
                            kymoselectionT=np.where(data_calculated["time"]>=deltaT)

                            #print(data_calculated["kymo_radav_" + varnames[iV]].shape)
                            #433 is num timesteps

                            #curkymodata = np.transpose(data_calculated["kymo_radav_" + varnames[iV]][kymoselectionT,:])
                            #print(curkymodata.shape)

                            aspect=2
                            if mode=="linear":
                                if varnames[iV] in ["u","w"]:



                                    imK = ax[iV,1].imshow(curkymodata, vmin=vmin*0.0001, vmax=vmax*0.1, origin='lower',
                                        extent=extent, aspect=aspect,cmap=cmap)
                                else:
                                    imK = ax[iV,1].imshow(curkymodata, vmin=vmin, vmax=vmax, origin='lower',
                                        extent=extent, aspect=aspect,cmap=cmap)
                            elif mode=="log":
                                #varnames = ["c2", "c3", "c4", "u", "w"]
                                #varnames_long = ["nutrients (mM)", "bacteria (OD)", "SCFA (mmol)", "$v_{r}$", "$v_{z}$"]

                                if varnames[iV]=="c2": #nutrients
                                    lower=0.001
                                    upper=250
                                    curarray=data_calculated["kymo_radav_"+varnames[iV]][kymoselectionT,:].squeeze()

                                    curarray[curarray <lower] = lower

                                    curkymodata = np.transpose(curarray)
                                    imK = ax[iV,1].imshow(curkymodata, norm=matplotlib.colors.LogNorm(vmin=lower, vmax=upper), origin='lower',
                                        extent=extent, aspect=aspect,cmap=cmap)
                                elif varnames[iV]=="c3": #bacteria
                                    lower=0.001
                                    upper=15
                                    curarray=data_calculated["kymo_radav_"+varnames[iV]][kymoselectionT,:].squeeze()
                                    curarray[curarray <lower] = lower
                                    curkymodata = np.transpose(curarray)

                                    imK = ax[iV,1].imshow(curkymodata, norm=matplotlib.colors.LogNorm(vmin=lower, vmax=upper), origin='lower',
                                        extent=extent, aspect=aspect,cmap=cmap)
                                elif varnames[iV] == "c4":  # scfa
                                    lower = 0.001
                                    upper = 500

                                    curarray = data_calculated["kymo_radav_" + varnames[iV]][kymoselectionT,:].squeeze()
                                    curarray[curarray < lower] = lower
                                    curkymodata = np.transpose(curarray)
                                    imK = ax[iV, 1].imshow(curkymodata, norm=matplotlib.colors.LogNorm(vmin=lower, vmax=upper),
                                                           origin='lower',
                                                           extent=extent, aspect=aspect, cmap=cmap)

                                elif  varnames[iV] in ["u","w"]:
                                    if no_massflow:
                                        lower=vmin
                                        upper=vmax
                                    else:
                                        lower = vmin * 0.0001
                                        upper = vmax * 0.1

                                    if boundary_option == "minmax":
                                        lower = var_list_interpolated2[iV].min()
                                        upper = var_list_interpolated2[iV].max()

                                    curarray = data_calculated["kymo_radav_" + varnames[iV]][kymoselectionT,:].squeeze()
                                    curkymodata = np.transpose(curarray)

                                    imK = ax[iV,1].imshow(curkymodata*conversionv, vmin=lower*conversionv, vmax=upper*conversionv, origin='lower',
                                        extent=extent, aspect=aspect,cmap=cmap)


                                    #or try symlog norm
                                    #imK = ax[0,iV].imshow(data_calculated[+varnames[iV]],  norm=matplotlib.colors.SymLogNorm(vmax*0.01,vmin=vmin, vmax=vmax), origin='lower',
                                    #    extent=[z.min()*100, z.max()*100,t_unique.min(), t_unique.max()], aspect=0.2,cmap=cmap)
                                else:
                                    error_Variablenotknown


                            ax[iV,1].set_xlabel("time (h)")
                            ax[iV,1].set_ylabel("z position (cm)")
                            ax[iV,1].axvline(t_cur-deltaT,ls='--',color='w')
                            ####################


                            #################
                            #plot k full spatial resolution
                            ##################

                            extent = [-1 * r.max() * 100, r.max() * 100, z.min() * 100, z.max() * 100]


                            aspect = .3 / (0.35 * r.max() * 100) #for wider plot smaller value
                            if no_massflow:
                                aspect=.3 / (0.35 * r.max() * 100) / 4.5




                            if mode=="linear":
                                if varnames[iV] in ["u", "w"]:
                                    im = ax[iV,0].imshow(var_list_interpolated2[iV], vmin=vmin * 0.0001, vmax=vmax * 0.1,
                                                               origin='lower',
                                                               extent=extent, aspect=aspect, cmap=cmap)
                                else:
                                    im = ax[iV,0].imshow(var_list_interpolated2[iV], vmin=vmin, vmax=vmax, origin='lower',
                                                               extent=extent, aspect=aspect, cmap=cmap)
                            elif mode=="log":
                                if iV==0: #nutrients
                                    lower=0.001
                                    upper=250
                                    if boundary_option=="minmax":
                                        lower = var_list_interpolated2[iV].min()
                                        upper = var_list_interpolated2[iV].max()
                                    curarray=var_list_interpolated2[iV]
                                    #print(curarray.min())
                                    curarray[curarray <lower] = lower
                                    #print(curarray.min())
                                    #print(lower)
                                    #lll
                                    im = ax[iV,0].imshow(curarray, norm=matplotlib.colors.LogNorm(vmin=lower*0.99, vmax=upper), origin='lower',
                                       extent=extent, aspect = aspect,cmap=cmap)

                                elif iV==1: #bacteria
                                    lower=0.01
                                    upper=15
                                    if boundary_option == "minmax":
                                        lower = var_list_interpolated2[iV].min()
                                        upper = var_list_interpolated2[iV].max()
                                    curarray=var_list_interpolated2[iV]
                                    curarray[curarray <lower] = lower
                                    im = ax[iV,0].imshow(curarray, norm=matplotlib.colors.LogNorm(vmin=lower, vmax=upper), origin='lower',
                                       extent=extent, aspect = aspect,cmap=cmap)
                                elif varnames[iV] =="c4":  # bacteria
                                    lower = 0.001
                                    upper = 500
                                    if boundary_option == "minmax":
                                        lower = var_list_interpolated2[iV].min()
                                        upper = var_list_interpolated2[iV].max()
                                    curarray = var_list_interpolated2[iV]
                                    curarray[curarray < lower] = lower
                                    im = ax[iV, 0].imshow(curarray,
                                                          norm=matplotlib.colors.LogNorm(vmin=lower, vmax=upper),
                                                          origin='lower',
                                                          extent=extent, aspect=aspect, cmap=cmap)


                                else: #velocities
                                    if no_massflow:
                                        lower=vmin
                                        upper=vmax
                                    else:
                                        lower = vmin * 0.0001
                                        upper = vmax * 0.1

                                    if boundary_option == "minmax":
                                        lower = var_list_interpolated2[iV].min()
                                        upper = var_list_interpolated2[iV].max()


                                    #curarray = var_list_interpolated2[iV]
                                    #curarray[curarray < lower] = lower

                                    im = ax[iV,0].imshow(var_list_interpolated2[iV]*conversionv, vmin=lower*conversionv, vmax=vmax*0.1*conversionv, origin='lower',
                                       extent=extent, aspect = aspect,cmap=cmap)

                                    #or try SymLogNorm
                                    #im = ax[1,iV].imshow(var_list_interpolated[iV], norm=matplotlib.colors.SymLogNorm(vmax*0.01,vmin=vmin, vmax=vmax), origin='lower',
                                    #   extent=[z.min()*100, z.max()*100,r.min()*1000, r.max()*1000], aspect = 0.35*r.max()*100,cmap=cmap)

                            ax[iV,0].axvline(0,ls='--',color='white')


                            #add scatter plot in addition (points should not be visible)
                            #ax[1,iV].scatter(r, z, marker="x",c=var_list[iV],alpha=0.2,vmin=vmin, vmax=vmax,cmap=cmap)
                            #ax[1,iV].scatter(r*1000, z*100, marker="x",c='k',alpha=0.2) #c=var
                            #plot scatter - with color code intensity

                            #add titles
                            #ax[iV,0].set_title(varnames_long[iV])
                            #if iV==0:
                            #ax[iV,0].set_title("time: "+str(round(t_cur-deltaT,2))+" h")

                            fig.suptitle("\n\n\n\n time: "+str(round(t_cur-deltaT,2))+" h", fontsize=20)

                            ax[iV,0].set_xlabel("r position (cm)")
                            ax[iV,0].set_ylabel("z position (cm)")

                            #colorbar
                            #try with divider
                            #fig.colorbar(im, ax=ax[1,iV],orientation="horizontal")
                            #divider.append(make_axes_locatable(ax[0,iV]))
                            #cax.append(divider[-1].append_axes("top", 0.4, pad="5%"))

                            #cax = divider.new_vertical(size = '2%', pad = 1)
                            #fig.add_axes(cax[-1])
                            #auto location (no top option)
                            #fig.colorbar(im, cax = cax[-1], orientation = 'horizontal')
                            #cb2 = colorbar(im, cax=cax[-1], orientation="horizontal")

                            #Give the colorbar its own axis to avoid resizing the parent axis:

                            if "c4" in varnames:
                                vertical_position = 1. - 0.15 * (iV + 1) - 0.088
                                height = 0.12
                                width = 0.025
                            else:
                                vertical_position = 1.-0.2*(iV+1)-0.088
                                height = 0.17
                                width = 0.025
                            horizontal_position = 0.82
                            cax.append(plt.axes([horizontal_position, vertical_position, width, height])) #the new axis for first colorbar
                            cb2 = fig.colorbar(imK, cax=cax[-1], orientation="vertical")
                            cb2.set_label(varnames_long[iV],fontsize=20,labelpad=15)
                            cax[-1].set_label(varnames_long[iV])
                            #cax[-1].set_yticks([])
                            #cax[-1].minorticks_off()

                            #adjustment of tick locations for SymLogNorm
                            #SymLogNorm needs adjustment of axes
                                    #https://stackoverflow.com/questions/11138706/colorbar-for-imshow-centered-on-0-and-with-symlog-scale
                            #need to fix
                            #if 3>5 and mode=="log" and iV in [2,3]:
                            #    maxlog=int(np.ceil( np.log10(vmax) ))
                            #    minlog=int(np.ceil( np.log10(-vmin) ))
                            #    #generate logarithmic ticks
                            #    tick_locations=([-(10**x) for x in xrange(minlog,-logthresh-1,-1)]
                            #        +[0.0]
                            #        +[(10**x) for x in xrange(-logthresh,maxlog+1)] )
                            #    cax[-1].set_xticks(tick_locations)

                            #plt.colorbar(im,cax=axColor,orientation='horizontal')

                            ax[iV,0].set_ylim(0,25)
                            #axColor.yaxis.set_label_position('left')
                            #axColor.yaxis.set_ticks_position('left')


                        #plt.show()


                        #fig.colorbar(im, ax=ax[1,iV])
                        #fig.tight_layout()
                        #plt.tight_layout(h_pad=2)


                        fig.savefig(os.path.join(foldernameout,timefilename), dpi=150)

                        #plt.show()x