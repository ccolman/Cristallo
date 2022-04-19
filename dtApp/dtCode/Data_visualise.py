from flask import render_template, request, redirect, url_for,send_file
from dtApp import app
from dtApp import date
import plotly
import plotly.graph_objs as go
import numpy as np
import json,time
from pathlib import Path
from openpyxl import Workbook   
from dtLib.lpmupdate.circleFitting import circle_fit
from plotly.subplots import make_subplots

@app.route('/Visualise',methods=['GET','POST'])
def Vis():
    # Read in Time History
    path='dtApp/dtData/data_th_impulse.csv'
    ModeRanges=[[5,7.5],[17.5,19],[25,28]]
    f=open(path)
    TH = np.genfromtxt(f, delimiter=',',skip_header=1)
    f.close()
    time = np.array([x[0] for x in TH]) # Time Vector
    h1 = np.array([x[1] for x in TH]) # Velocity floor 1 [m/s]
    h2 = np.array([x[2] for x in TH]) # Velocity floor 2 [m/s]
    h3 = np.array([x[3] for x in TH]) # Velocity floor 3 [m/s]
    
    DATA=[h1,h2,h3]
    dt=time[1]-time[0]
    sf=1/dt
    LDATA=len(DATA)
    # Convert to FRF
    FRF=[]
    for i in range(LDATA):
        Y=DATA[i]
        resp=np.fft.fft(Y)
        freq=np.fft.fftfreq(time.shape[-1],d=1/sf)
        resp,freq=resp[:int(len(resp)/2)],freq[:int(len(resp)/2)]
        FRF.append(resp)
    # Get linear frequency
    freq = np.real(freq)    # frequency vector [Hz]
    w = 2*np.pi*freq
    dfreq=freq[1]-freq[0]
    theta = np.linspace(0, 2*np.pi, num=1000)
    delta_freq_point = 1
    
    # Circle Fit
    Mode_shape,Damping,Frequency=[],[],[]
    for ModeIndex in range(len(ModeRanges)):
        freq_mode = freq[int(ModeRanges[ModeIndex][0]/dfreq):int(ModeRanges[ModeIndex][1]/dfreq)]
        
        r_mode,xc_mode,yc_mode=np.zeros(LDATA),np.zeros(LDATA),np.zeros(LDATA)
        fn_mode,eta_mode,shape_mode=np.zeros(LDATA),np.zeros(LDATA),np.zeros(LDATA)
        for i in range(LDATA):
            y = FRF[i]
            y_mode = y[int(ModeRanges[ModeIndex][0]/dfreq):int(ModeRanges[ModeIndex][1]/dfreq)]
            y_mode_real = np.real(y_mode)
            y_mode_imag = np.imag(y_mode)
            modalfit= circle_fit(y_mode_real, y_mode_imag) # fit circle to mobility of sensor 1 around mode 1
            xc_mode[i] = modalfit['xc']
            yc_mode[i] = modalfit['yc']
            r_mode[i] = modalfit['r']
            # Natural Frequency 
            idx = np.argmax(np.absolute(y_mode))
            fn_mode[i] = freq_mode[idx]
            # Damping Ratio
            fc = freq_mode[idx-delta_freq_point]
            fd = freq_mode[idx+delta_freq_point]
            thetac = np.angle(y_mode[idx])-np.angle(y_mode[idx-delta_freq_point])
            thetad = np.angle(y_mode[idx])-np.angle(y_mode[idx+delta_freq_point])
            eta_mode[i] = np.absolute(fc**2-fd**2)/(fn_mode[i]**2*(np.absolute(np.tan(thetac/2))+np.absolute(np.tan(thetad/2))))
        # Mode Shape
        for i in range(LDATA):
            shape_mode[i]=((2*np.pi*np.mean(fn_mode))**2)*np.mean(eta_mode)*2*r_mode[i]*np.sign(xc_mode[i])
            
        # Normalize and store
        Mode_shape.append(shape_mode/np.linalg.norm(shape_mode, ord=np.inf))
        Frequency.append(fn_mode)
        Damping.append(eta_mode*100)
    # Plot 1 - FRF
    fig = make_subplots(rows=1, cols=1)
    for i in range(LDATA):
        fig.add_scatter(x=freq,y=np.abs(FRF[i]), name='Frequency Response Function', mode = 'lines', row=1, col=1)
    fig.update_xaxes(range=[3,30])
    fig.update_xaxes(title_text="Frequency [Hz]")
    plot1 = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    
    # Plot 2 - Mode Shapes
    fig = make_subplots(rows=1, cols=3)
    for i in range(3):
        fig.add_scatter(y=[0,1,2,3],x=[0]+Mode_shape[i].tolist(), mode = 'lines', row=1, col=i+1)
    plot2 = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    
    # Send
    form={'Wn':np.mean(np.array(Frequency),axis=1),'Zeta':np.mean(np.array(Damping),axis=1),
          'plot1':plot1,'plot2':plot2}
    return render_template('TempViewData.html', date=date, form=form)


@app.route('/Vis_save',methods=['GET','POST'])
def Vis_save():
    
    return render_template('home.html', date=date)