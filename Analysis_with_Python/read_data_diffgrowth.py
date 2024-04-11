#
import os
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
#set filename


filename="26_11.txt"

#########
#start main script - no adjustment should be required
########

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
Dlist=[]
M1list=[]
colname=[]
varlist=[]
Lambdalist=[]
for il in range(0,len(columnname)):
    curstr=columnname[il]
    varname=curstr.split(" ")[0]
    varlist.append(varname)
    timel.append(np.nan)
    Dlist.append(np.nan)
    Lambdalist.append(np.nan)
    M1list.append("no_entry")
    if il>1:
        #print(curstr)
        parinfo=curstr.split("@ ")[1].split(", ")
        for item in parinfo:
            itemsplit=item.split("=")
            if itemsplit[0]=="t":
                timel[-1]=itemsplit[1]
            if itemsplit[0]=="D_AC":
                Dlist[-1]=itemsplit[1]
            if itemsplit[0]=="Lambda":
                Lambdalist[-1]=itemsplit[1]
            if itemsplit[0]=="M1":
                M1list[-1]=itemsplit[1].strip()

    #    lll
    #        
    #except:
    #    pass
        colname.append(varname+"t"+str(timel[-1])+"Lambda"+str(Lambdalist[-1])+"D"+str(Dlist[-1]))
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
#varunq = (list(varunq))
Dunq = list(set(Dlist[2:]))
Lambdaunq = list(set(Lambdalist[2:]))
M1unq = list(set(M1list[2:]))

#r=data.iloc[:, 0]
#z=data.iloc[:, 1]
columnnames=data.columns
firstrun=True
#dataout=pd.DataFrame(columns=["r","z","t","D","M1"]+varunq)
#do this seperately for all D and M values to make sure concat to merge dataframes does not take too long


print("Going through different parameter combinations....")
for Lambda in Lambdaunq:
    print("Lambeda: "+str(Lambda))
    for D in Dunq:
        print("D: "+str(D))
        for M1 in M1unq:
            print("M1: "+str(M1))
            for t in tunq:
                #print("D: "+str(D)+"M1: "+str(M1)+str(t))
                firstvar=True
                #for unique combination of t,D,M1, merge all variables into one table
                for var in varunq:
                    cc=-1
                    for column in columnnames:
                        cc=cc+1
                        if cc>1:
                            if M1list[cc]==M1 and Lambdalist[cc]==Lambda and Dlist[cc]==D and varlist[cc]==var and timel[cc] == t:
                                curdata=data[[columnnames[0],columnnames[1],column]]
                                curdata.columns = ['R', 'z', varlist[cc]]
                                curdata['M1']= M1list[cc]
                                curdata['D']= Dlist[cc]
                                curdata['Lambda'] = Lambdalist[cc]
                                curdata['t']= timel[cc]
                                if firstvar==True: #when running combination of D and M1 for first time, use concat. After that merge (to add other variables)
                                    datanew2=curdata.copy()
                                    firstvar=False
                                else:
                                    datanew2=pd.merge(datanew2, curdata,  how='outer', on=['R','z',"t","D","Lambda","M1"])
                                    #print("added")
                if firstrun:
                    dataout=datanew2.copy()
                    firstrun=False
                else:
                    #dataout=dataout.append(datanew2,ignore_index=True)
                    dataout = pd.concat([dataout,datanew2], ignore_index=True,axis=0, join='outer')

                    #display(dataout)
            #print(data[column])

        #DataFrame.append(other, ignore_index=False, verify_integrity=False, sort=None)[source]¶
if "r" in dataout.columns:
    pass
else:
    dataout.rename(columns={"R": "r"},inplace=True)
        
print("Saving dataframe to file...(might take a few min)")
filepath=os.path.join(folder_output,"simpleformat_"+filename)
dataout.to_csv(filepath,index=False)
print("output saved to: "+filepath)
print("output datafraome: ")
