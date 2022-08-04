from __future__ import absolute_import, division, print_function
from builtins import *  # @UnusedWildImport
from mcculw import ul
from mcculw.enums import ScanOptions, FunctionType, Status, ULRange
#from mcculw.device_info import DaqDeviceInfo
from mcculw.enums import InterfaceType
from mcculw.enums import ChannelType
from mcculw.ul import get_status, daq_in_scan

from ctypes import cast, POINTER,  c_ulong

import time
from time import sleep # use Qtimer.timeout?
#import sys

import numpy as np
from numpy import array, arange
from numpy.ctypeslib import as_array

from PyQt5.QtCore import pyqtSignal, pyqtSlot
from pyqtgraph import PlotWidget
from PyQt5 import QtWidgets
from PyQt5.QtCore import QObject, QThread
import os

from tkinter import Tk
from tkinter.filedialog import askdirectory

try:
    import OpenGL
    pg.setConfigOption('useOpenGL', True)
    pg.setConfigOption('enableExperimental', True)
except Exception as e:
    print(f"Enabling OpenGL failed with {e}. Will result in slow rendering. Try installing PyOpenGL.")



class Worker(QObject):
    
    finished = pyqtSignal()
    datahandle = pyqtSignal(int) # int in py3 are unlimited size
    
    @pyqtSlot()
    def getdata(self):
        global status, curr_index, plot_limit

        while status != Status.IDLE:
            if curr_index > plot_limit:
                # use curr_index, because range(..., curr_index) ensures that it stops at last channel
                # dont have to reshape when using vispy (just match the n in vispy to plot_limit
                # currently assume ctype_array is already global, we will see after that
                self.datahandle.emit(curr_index)
                sleep(0.1) # for 200ksps, recommend plot at this speed     
            status, curr_count, curr_index = get_status(0, FunctionType.DAQIFUNCTION)
            
        self.finished.emit()

    
class PyQtGraphTest(PlotWidget):
    # Signal to indicate new data acquisition
    # Note: signals need to be defined inside a QObject class/subclass
    def __init__(self,app,ctypes_array,points_per_ch,num_ch,total_count):
    
        super().__init__()
        
        self.ctypes_array = ctypes_array
        self.qapp = app
        self.total_count = total_count
        self.num_ch = num_ch
        self.points_per_ch = points_per_ch
        self.framelength = num_ch*points_per_ch
        self.buffer = [0 for _ in range(self.framelength)]
        self.setWindowTitle('signals')
        self.resize(1080, 720)
        #self.setYRange(40000,160000, padding = 0)
        self.setMouseEnabled(x=False, y=False)
        self.x_range = arange(points_per_ch)
        
        self.colormap = ['c','g','y','r','w','b']
        self.spectrum = [self.plot(pen = self.colormap[i],downsample=4) for i in range(num_ch)]

        self.t = time.time()
        # Connect the signal   
        # 1 - create Worker and Thread inside the Form
        self.obj = Worker()  # no parent!
        self.thread = QThread()  # no parent!   
        # 2 - Connect Worker`s Signals to Form method slots to post data.
        self.obj.datahandle.connect(self.onDataReady)   
        # 3 - Move the Worker object to the Thread object
        self.obj.moveToThread(self.thread)   
        # 4 - Connect Worker Signals to the Thread slots
        self.obj.finished.connect(self.thread.quit)  
        # 5 - Connect Thread started signal to Worker operational slot method
        self.thread.started.connect(self.obj.getdata)   
        # * - Thread finished signal will close the app if you want!
        self.thread.finished.connect(self.onFinished)   
        # 6 - Start the thread
        self.thread.start()
        
    def onFinished(self):

        self.close() # avoid window close error
        self.qapp.quit() #exit the app to advance to next line in run_aquisition
         
    def onDataReady(self,index):
        self.buffer = array(as_array(self.ctypes_array,(self.total_count,)))[index-self.framelength:index].reshape(self.points_per_ch,self.num_ch)
        self.setTitle(str(int(time.time()-self.t))+' s')
        for i in range(self.num_ch):
            self.spectrum[i].setData(x=self.x_range,y=self.buffer[:,i])
        
        
# since these are methods, may be run them in "init" of canvas, before self.show?
def config_first_detected_device(board_num, dev_id_list=None):
    """Adds the first available device to the UL.  If a types_list is specified,
    the first available device in the types list will be add to the UL.

    Parameters
    ----------
    board_num : int
        The board number to assign to the board when configuring the device.

    dev_id_list : list[int], optional
        A list of product IDs used to filter the results. Default is None.
        See UL documentation for device IDs.
    """
    ul.ignore_instacal()
    devices = ul.get_daq_device_inventory(InterfaceType.ANY)
    if not devices:
        raise Exception('Error: No DAQ devices found')

    device = devices[0]
    if dev_id_list:
        device = next((device for device in devices
                       if device.product_id in dev_id_list), None)
        if not device:
            err_str = 'Error: No DAQ device found in device ID list: '
            err_str += ','.join(str(dev_id) for dev_id in dev_id_list)
            raise Exception(err_str)

    # Add the first DAQ device to the UL with the specified board number
    ul.create_daq_device(board_num, device)

def run_aquisition(filename, sample_rate = 100000, seconds = 60, minutes = 1, channel_list = [0,1,2,3,4,5], ppf = 1000):
    global  status, curr_index, plot_limit
    # By default, the example detects and displays all available devices and
    # selects the first device listed. Use the dev_id_list variable to filter
    # detected devices by device ID (see UL documentation for device IDs).
    # If use_device_detection is set to False, the board_num variable needs to
    # match the desired board number configured with Instacal.
    use_device_detection = True
    dev_id_list = []
    board_num = 0
    rate = sample_rate
    # if we input 60s, the actual recording period is 57.6s 
    # note to calculate this, was 7680, for 80k 5ch
    # https://www.mccdaq.de/pdfs/manuals/Mcculw_WebHelp/Users_Guide/Analog_Input_Boards/USB-1808.htm
    points_per_channel = seconds * rate

    try:
        if use_device_detection:
            config_first_detected_device(board_num, dev_id_list)

        # daq_dev_info = DaqDeviceInfo(board_num)
        # if not daq_dev_info.supports_analog_input:
        #     raise Exception('Error: The DAQ device does not support '
        #                     'analog input')
        num_chans = len(channel_list)
       
        # daq_in_scan inputs
        chan_list = channel_list
        chan_type_list = [ChannelType.ANALOG_DIFF for i in range(num_chans)]
        if 5 not in chan_list:
            gain_list = [ULRange.UNI10VOLTS for i in range(num_chans)]
        else:
            gain_list = [ULRange.UNI10VOLTS for i in range(num_chans-1)]
            gain_list.append(ULRange.UNI5VOLTS)
                
        
        pretrig_count = 0
        total_count = points_per_channel * num_chans 
        
        
        # here, define how many pts you want to see in one plot
        # if hangs or too fast, try undsersampling
        plot_limit = num_chans*ppf
                    
        # Use the win_buf_alloc_32 method for devices with a resolution > 16
        memhandle = ul.win_buf_alloc_32(total_count)
        # Convert the memhandle to a ctypes array.
        ctypes_array = cast(memhandle, POINTER(c_ulong))

        # Note: the ctypes array will no longer be valid after win_buf_free is
        # called.
        
        # Check if the buffer was successfully allocated
        if not memhandle:
            raise Exception('Error: Failed to allocate memory')                       
        
        # Start the scan
        daq_in_scan(
            board_num, chan_list,chan_type_list, gain_list, num_chans,
            rate, pretrig_count, total_count , memhandle, ScanOptions.BACKGROUND)
        
        # Start updating the displayed values
        status, curr_count, curr_index = get_status(
            board_num, FunctionType.DAQIFUNCTION)
               
        # we should call the class of plot over here, after creating memhandle
        app = QtWidgets.QApplication([])
        plot_window = PyQtGraphTest(app,ctypes_array,points_per_ch=ppf,num_ch=num_chans,total_count=total_count)
        plot_window.show()        
        app.exec_()        # start and block at the line, crucial
   
        # Stop the background operation (this is required even if the
        # scan completes successfully)
        ul.stop_background(board_num, FunctionType.DAQIFUNCTION)        
        
        # at last, get data from memory, this line wont change the ctypes_array elements
        # repeat it result in copy of same array
        np.save((filename),array(as_array(ctypes_array,(total_count,))).reshape(-1,num_chans))
        
        
    except Exception as e:
        print('\n', e)
    #add keyboard interruption here, when run in spyder, use the red stop recording button in kernel
    except KeyboardInterrupt:
        print("Keyboard interrupted")
    finally:
        if memhandle:
            # Free the buffer in a finally block to prevent a memory leak.
            ul.win_buf_free(memhandle)
        if use_device_detection:
            ul.release_daq_device(board_num)
        

if __name__ == '__main__':
    # start this process use run file directly, restart kernel before running
    #define your "global" 
    Tk().withdraw()
    path = askdirectory(title='Select Folder') # shows dialog box and return the path
    os.chdir(path) # use this to set the path already
    sample_rate = 1000 *200
    
    seconds = 60
    
    minutes = 4 # only accept integer
    
    channel_list = [0,1,2,3,4,5]
    
    #channel_list = [0]
    
    ppf = 1000    
    
    logname = 'Onlygzb_2' + '_'
    
    t = time.time()
    
    for i in range(minutes):
        filename = logname + str(i)
        run_aquisition(filename,sample_rate = sample_rate, 
                       seconds = seconds, minutes = minutes,
                       channel_list=channel_list,ppf=ppf)
        
        
    print('Scan completed successfully in ', time.time()-t, ' seconds')
    

    
    
