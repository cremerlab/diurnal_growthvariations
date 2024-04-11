#
import os
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
#set filename

filename="26_3_SCFA.txt"
D=np.power(10.,-8.)
M1=1

#read in header and extract information
folder_comsoldata="data_fromcomsol/full_variables_spatiotemporal_resolution"
folder_output="data_simpleformat"

print(filename)
fp = open(os.path.join(folder_comsoldata,filename))

for i, line in enumerate(fp):
    if i == 8:
        headstr=line
    elif i > 8:
        break
fp.close()

headstr=headstr.replace('% r;', '% R;', 1)



columnname=headstr[2:].split(";")


timel=[]

colname=[]
varlist=[]

for il in range(0,len(columnname)):
    curstr=columnname[il]
    varname=curstr.split(" ")[0]
    varlist.append(varname)
    timel.append(np.nan)
    if il>1:
        #print(curstr)
        parinfo=curstr.split("@ ")[1].split(", ")
        for item in parinfo:
            itemsplit=item.split("=")
            if itemsplit[0]=="t":
                timel[-1]=itemsplit[1]

    #    lll
    #        
    #except:
    #    pass
        colname.append(varname+"t"+str(timel[-1])) #+"D"+str(Dlist[-1])
    else:
        colname.append(varname)

#read in data
data=pd.read_csv(os.path.join(folder_comsoldata,filename),skiprows=9,delimiter=";",names=colname)
print("File opend with shape:")
print(data.shape)


#generate dataframe with variables in different columns
#unique variables:
varunq = list(set(varlist[2:]))
tunq = list(set(timel[2:]))
print(len(tunq))

columnnames=data.columns
firstrun=True
#do this seperately for all D and M values to make sure concat to merge dataframes does not take too long


print("Going through different parameter combinations....")

for t in tunq:
    #print("D: "+str(D)+"M1: "+str(M1)+str(t))
    firstvar=True
    #for unique combination of t,D,M1, merge all variables into one table
    for var in varunq:
        cc=-1
        for column in columnnames:
            cc=cc+1
            if cc>1:
                if varlist[cc]==var and timel[cc] == t:
                    curdata=data[[columnnames[0],columnnames[1],column]]
                    curdata.columns = ['R', 'z', varlist[cc]]
                    curdata['M1']= M1
                    curdata['D']= D
                    curdata['t']= timel[cc]
                    if firstvar==True: #when running combination of D and M1 for first time, use concat. After that merge (to add other variables)
                        datanew2=curdata.copy()
                        
                        firstvar=False
                    else:
                        datanew2=pd.merge(datanew2, curdata,  how='outer', on=['R','z',"t","D","M1"])
                        #print("added")
    if firstrun:
        dataout=datanew2.copy()
        firstrun=False
    else:
        dataout = pd.concat([dataout,datanew2], ignore_index=True,axis=0, join='outer')

        #DataFrame.append(other, ignore_index=False, verify_integrity=False, sort=None)[source]¶
if "r" in dataout.columns:
    pass
else:
    dataout.rename(columns={"R": "r"},inplace=True)
        
print("Saving dataframe to file...(might take a few min)")
filepath=os.path.join(folder_output,filename)
dataout.to_csv(filepath,index=False)
print("output saved to: "+filepath)
print("output datafraome: ")
