import pandas as pd
import datetime as dt
import math
from pyproj import Geod
from obspy import UTCDateTime
from obspy.core import Stream
from obspy.clients.fdsn import Client
client = Client("http://192.168.35.10:8080")

def getChanString(chan):
    if len(chan)==3:
        return chan[0]+chan[2]
    else:
        return chan

def right_justify(s, numSpace=3):
    return f"%{numSpace}s" % s

def toFixed(numObj, digits=0):
    return f"{numObj:.{digits}f}"

def CreateNordicFile(eventid, df):
    global stations
    global client
    delt_time=300
    part_df = df.loc[df['eventid']==eventid]
    firstElem = part_df.iloc[0]
    g=Geod(ellps="WGS84")
    lat = toFixed(firstElem.lat,3)
    lon = toFixed(firstElem.lon,3)
    depth = toFixed(firstElem.depth,2)
    dh = toFixed(firstElem.dh,2)
    ml = toFixed(firstElem.Ml,1)
    tOrigin = UTCDateTime(dt.datetime.fromtimestamp(firstElem.origintime))
    tStart = tOrigin-delt_time
    tEnd = tOrigin+delt_time
    fName = dt.datetime.fromtimestamp(firstElem.origintime).strftime("%d-%H%M-%SL.S%Y%m")
    timeStr = dt.datetime.fromtimestamp(firstElem.origintime).strftime("%Y %d%m %H%M %S.%f")[:20]
    st=Stream()
    with open(f"output/{fName}","w") as f:
        line = " "+timeStr +"L  " + lat + "  " + lon +" " + depth + "  MI   " + "8 .50" + " " + ml +"LMl  " + ml+" MI  " +ml+"LMI"+" 1\n"
        #f.write(line)
        line += " GAP=GGG   MI              2.0     2.00 3.2                                    E\n"
        #f.write(line)
        line += " 2014-03-13-1605-19M.TEST__015                                                 6\n"
        #f.write(line)
        line+=" STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU VELO AIN AR TRES W  DIS CAZ7\n"
        #f.write(line)

        ArrAz = []
        ArrSta = []
        for _,ph in part_df.iterrows():
            print(ph)

            if not (math.isnan(ph.itime)) and ph.iphase!="AML":
                staPart = stations.loc[stations['sta']==ph.sta]
                sta = staPart.iloc[0]

                azimuth, back_azimuth, distance_2d = g.inv(sta.lon, sta.lat, ph.lon, ph.lon)
                if azimuth<0:
                    azimuth=360+azimuth
                ArrAz.append(azimuth)
                #print(distance_2d)
                az = right_justify(toFixed(azimuth))
                dist = right_justify(toFixed(distance_2d/1000))
                time=dt.datetime.fromtimestamp(ph.itime)
                timeStr = time.strftime("%H%M %S.%f")[:10]
                line += " "+ph.sta.ljust(4,' ')+" "+getChanString(ph.chan)+" E"+ph.iphase+"   1   "+timeStr+f"                                            {dist} {az} \n"
                #f.write(line)
                st+=client.get_waveforms(network="*",location="*",station=ph.sta,channel=ph.chan,starttime=tStart,endtime=tEnd)

        ArrAz.sort()
        #print(ArrAz)
        maxGap = 0
        for i in range(1,len(ArrAz),1):
            delt = ArrAz[i]-ArrAz[i-1]
            if delt>maxGap:
                maxGap=delt
        delt = 360 - ArrAz[len(ArrAz)-1]
        if delt>maxGap:
            maxGap = delt
        #print(maxGap)
        maxGap = toFixed(maxGap)
        strMaxGap = right_justify(maxGap)
        line=line.replace("GGG",strMaxGap)
        f.write(line)
        st.write(f"msd/{fName}.msd", format="MSEED")



pth='ural_bul.csv'
pth_stat = 'stations.csv'
stations = pd.read_csv(pth_stat)
df = pd.read_csv(pth)
first=df.head(1)
eventid = first.eventid[0]

CreateNordicFile(eventid,df)

evtList = list(set(df.eventid))
print(evtList)
#for evt in evtList:
#    CreateNordicFile(evt,df)




