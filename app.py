#!/usr/bin/env python
# coding: utf-8

# In[1]:

import streamlit as st
import plotly.express as px
from attr import define
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import serial
import keyboard as kb
import time
import math
from scipy.signal import savgol_filter
import wfdb
import scipy
ind=0

st.set_page_config(page_title="BME 300 Capstone Project",
                    page_icon=":syringe:",
                    layout="wide")
# In[2]:


# Data Collection:
dat_file=r'C:\Drive\BME 300\icelandic-16-electrode-electrohysterogram-database-1.0.0\files\ice008_p_2of4'
#dat_file=r'C:\Drive\BME 300\icelandic-16-electrode-electrohysterogram-database-1.0.0\files\ice031_l_3of3'
wfdb.rdheader(dat_file)
ind+=1


# In[3]:


record=wfdb.rdrecord(dat_file,channels=[1,3])
# record


# In[4]:


signals, fields= wfdb.rdsamp(dat_file,sampfrom=200000,sampto=720000,channels=[1,2])
signals1, fields1= wfdb.rdsamp(dat_file,sampfrom=0,channels=[1,2])
# signals.shape


# In[5]:


# fields


# In[6]:


#for i in signals[:,1]:
    #print(i)
#plt.plot(signals[:,0]) 
#plt.plot(signals[:,1])


# In[7]:


secs=len(signals1[:,1])/(200)
mins=secs/60
y=signals[:,1]
y=list(y)
x=np.arange(0,secs,1/200)
#x.shape
#x=x[250000:670000,]
x=x[200000:720000,]
#x.shape
x=list(x)
len(y)


# In[8]:


y_unfiltered=y


# In[9]:


y=savgol_filter(y_unfiltered,4500,2)


#plt.figure(figsize=(10, 4))

#plt.subplot(221)
#plt.plot(x, y_unfiltered)
#plt.title("EHG Signal with Noise")
#plt.margins(0, .05)

#plt.subplot(222)
#plt.plot(x, y)
#plt.title("Filtered EHG Signal")
#plt.margins(0, .05)

#plt.tight_layout()
#plt.show()

ehg_y=list(y)
ehg_y_unfiltered=list(y_unfiltered)
ehg_np=[x,ehg_y,ehg_y_unfiltered]
ehg_np=np.array(ehg_np)
ehg_np=ehg_np.T
ehg_df=pd.DataFrame(ehg_np,columns=["Time","Filtered_EHG","Raw_EHG"])
# ehg_df

# In[10]:


#Data PreProcessing:
freq=200 # 20Hz sampling frequency
minutes=3
slice_time= minutes*60 # 3 minutes
slice_size= slice_time*freq

sample_size=len(x)
slice_num=math.floor(sample_size/slice_size)
x_sliced= [x[i*slice_size:(i+1)*slice_size] for i in range(0,slice_num)]
y_sliced= [y[i*slice_size:(i+1)*slice_size] for i in range(0,slice_num)]

x_sliced=np.array(x_sliced)
y_sliced=np.array(y_sliced)  #now time and voltage has been successfully divided into segments of 30 mins or 36,000 samples

def findPeakIndex(y):
    y=np.array(y)
    c=y.argmax()

    return c

def percDiff(val1,val2):
    valM=abs(val1-val2)
    valP=abs(val1+val2)/2
    perDiff=(valM/valP)*100
    return perDiff


# In[11]:


# Data Analysis:
final_Decision=list()
Threshold=70
for i in range(0,len(x_sliced)):
    FD=0
    PFC=0
    fPFC=0
    PRPD=0
    peaksFD=0
    freqFD=0
    t=x_sliced[i,:]
    volt_raw=y_sliced[i,:]
    volt=savgol_filter(volt_raw,4500,2)
    #print(volt)
    # print('i='+ str(i))
    # phase-1
    sd_threshold=0.1
    sd=np.std(volt)
    if sd<sd_threshold:
        FD=0;
    else:
        #phase-2
        percent=5
        segments=math.floor(100/percent)
        segment_size= (percent/100)*len(t)
        segment_size=int(segment_size)
        #print("segment size")
        #print(segment_size)
        t05=[t[p*segment_size:(p+1)*segment_size] for p in range(0,segments)]
        volt05=[volt[p*segment_size:(p+1)*segment_size] for p in range(0,segments)]
        #print(volt05)
        for j in range(0,len(t05)):
            #print('j='+ str(j))
            tp=t05[j]
            voltp=volt05[j]
            #print(len(tp))

            chunk_size=100
            chunk_num=int(len(tp)/chunk_size)
            #print("chunk num")
            #print(chunk_num)
            t_chunked=[tp[p*chunk_size:(p+1)*chunk_size] for p in range(0,chunk_num)]
            volt_chunked=[voltp[p*chunk_size:(p+1)*chunk_size] for p in range(0,chunk_num)]
            #print(volt_chunked)
            #Peaks Analysis
            for k in range(0,len(volt_chunked)-1):
                #print('k='+ str(k))
                #print(volt_chunked[k])
                #print(volt_chunked[k+1])
                ind1=findPeakIndex(volt_chunked[k])
                ind2=findPeakIndex(volt_chunked[k+1])
                peak1=volt_chunked[k][ind1]
                peak2=volt_chunked[k+1][ind2]
                PDR=percDiff(peak1,peak2)
                if PDR>14:
                    PFC+=1
            PRPD=(PFC*100)/len(volt_chunked) 
            if PRPD>25:
                peaksFD+=5
            
            #Frequency Analysis
            for k in range(0,len(volt_chunked)-2):
                ind1=findPeakIndex(volt_chunked[k])
                ind2=findPeakIndex(volt_chunked[k+1])
                ind3=findPeakIndex(volt_chunked[k+2])
                time1=t_chunked[k][ind1]
                time2=t_chunked[k+1][ind2]
                time3=t_chunked[k+2][ind3]
                freq1=time2-time1
                freq2=time3-time2
                fPDR=percDiff(freq1,freq2)
                if fPDR>14:
                    fPFC+=1
            fPRPD=(fPFC*100)/len(volt_chunked) 
            if fPRPD>25:
                freqFD+=5
        FD=(peaksFD+freqFD)
    #Final Decision Analysis
    # print(FD)
    if FD>=Threshold:
        final_Decision.append(1)
    else:
        final_Decision.append(0)
        
preterm_perc=sum(final_Decision)/len(final_Decision)*100
# print(preterm_perc)
        


# In[12]:


import peakutils
y=np.array(y)
y_baseline=peakutils.baseline(y)

baseline_deflection=(np.trapz(y_baseline))/len(y)
# baseline_deflection


# In[13]:


import biosppy
ehg_bio=biosppy.signals.emg.emg(y)
ehg_freqs=len(ehg_bio['onsets'])/x[-1]*60
# ehg_freqs


# In[18]:


freq=200 # 20Hz sampling frequency
minutes=10
slice_time= minutes*60  # 10 minutes
slice_size= slice_time*freq

sample_size=len(x)
slice_num=math.floor(sample_size/slice_size)


# In[19]:


y_pressure=38.71-2.764*y+0.7966*y**2
y_pressure_baseline=peakutils.baseline(y_pressure)
y_dist=y_pressure-y_pressure_baseline
y_dist_sliced= [y_dist[i*slice_size:(i+1)*slice_size] for i in range(0,slice_num)]
MVUarray=[]
for i in range(0,len(y_dist_sliced)):
    MVUarray.append(sum(y_dist_sliced[i])/slice_time)
MVUarray=np.array(MVUarray)
MVU=np.mean(MVUarray)


# In[20]:


patient= {}
patient['Preterm_Contraction']=preterm_perc
patient['Baseline_Deflection']=baseline_deflection
patient['EHG_frequency']=ehg_freqs
patient['Montevideo_Units']=MVU


# In[21]:


# patient


# In[87]:


patient_df=pd.DataFrame(patient,index=[ind])
# patient_df


# In[88]:


## patient_df.to_csv(r'C:\Drive\BME 300\Feature_Extraction_using_Icelandic_Dataset.csv',mode='a',index=True,header=False)


# In[89]:

st.title(":syringe: Preterm Detector")
st.markdown("##")
st.dataframe(patient_df)
st.markdown("---")


preterm=patient_df["Preterm_Contraction"].mean()
baseline_def=patient_df["Baseline_Deflection"].mean()
freq_ehg=patient_df["EHG_frequency"].mean()
Montevideo=patient_df["Montevideo_Units"].mean()

preterm=round(preterm,2)
freq_ehg=round(freq_ehg,2)
Montevideo=round(Montevideo,2)
baseline_def=round(baseline_def,2)

col1, col2 = st.columns(2)

with col1:
    st.subheader("Percentage of Contractions that indicate preterm labor:")
    st.subheader(f"{preterm} %")
with col1:
    st.subheader("Baseline Deflection Observed from EHG:")
    st.subheader(f"\n {baseline_def}")
with col2:
    st.subheader("Frequency of onsets in EHG/Uterine EMG:")
    st.subheader(f"{freq_ehg} Hz")
with col2:
    st.subheader("Montevideo Units calculated from Uterine Pressure graph:")
    st.subheader(f"{Montevideo} mm Hg")

st.markdown("---")

fig_raw = px.line(ehg_df, x="Time", y="Raw_EHG", title='Raw EHG curve')
fig = px.line(ehg_df, x="Time", y="Filtered_EHG", title='Filtered EHG curve')

col01, col02 = st.columns(2)
with col01:
    st.plotly_chart(fig_raw)
with col02:
    st.plotly_chart(fig)

st.markdown("---")




# In[ ]:

print("Done")


