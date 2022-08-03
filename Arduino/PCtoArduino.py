# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 12:18:30 2021

@author: luri1
"""

import sys
from PyQt5.QtWidgets import QGridLayout, QGroupBox, QLabel, QApplication, QPushButton, QWidget, QHBoxLayout
import serial
import time
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from PyQt5.QtCore import QObject, QThread
from time import sleep

class Worker(QObject):    
    global arduino
    datahandle = pyqtSignal(str) # int in py3 are unlimited size   
    @pyqtSlot()
    def getdata(self):
        while True:
            if arduino.in_waiting != 0: #
                decoded_bytes = arduino.readline().decode("utf-8")
                self.datahandle.emit(decoded_bytes)
            sleep(0.5) # for 200ksps, recommend plot at this speed                

class PCtoArduino(QWidget):
    def __init__(self, qapp,arduino):
        super().__init__()
        self.setWindowTitle('Select All/Select None')
        self.window_width, self.window_height = 1200, 800
        self.setMinimumSize(self.window_width, self.window_height)
        self.status_text = ["Closed","Open"]
        self.v1_status = False
        self.v2_status = False
        self.v3_status = False
        self.v4_status = False
        self.pressure = 0
        self.app = qapp
        self.arduino = arduino
        
        #button box, for PC sender
        grid_A = QGridLayout()

        self.func1_btn = QPushButton("1")
        self.func2_btn = QPushButton("2")
        self.func3_btn = QPushButton("3")
        self.func4_btn = QPushButton("4")
        self.P_release_btn = QPushButton("Release ALL")
        self.P_close_btn = QPushButton("Close ALL")
        
        self.func1_btn.clicked.connect(self.button1pressed)
        self.func2_btn.clicked.connect(self.button2pressed)
        self.func3_btn.clicked.connect(self.button3pressed)
        self.func4_btn.clicked.connect(self.button4pressed)
        self.P_release_btn.clicked.connect(self.releasepressed)
        self.P_close_btn.clicked.connect(self.closepressed)

        gb_a = QGroupBox(self)
        gb_a.setTitle("Send singal to Arduino")
        grid_A.addWidget(self.func1_btn, 0, 0)
        grid_A.addWidget(self.func2_btn, 0, 1)
        grid_A.addWidget(self.func3_btn, 1, 0)
        grid_A.addWidget(self.func4_btn, 1, 1)
        grid_A.addWidget(self.P_release_btn,2,0,1,-1)
        grid_A.addWidget(self.P_close_btn,3,0,1,-1)
        gb_a.setLayout(grid_A)
        
        
        
        # arduino monitor        
        grid_B = QGridLayout()
        lab_1 = QLabel("Valve 1 Status: ")
        lab_2 = QLabel("Valve 2 Status: ")
        lab_3 = QLabel("Valve 3 Status: ")
        lab_4 = QLabel("Valve 4 Status: ")
        lab_P = QLabel("Pressure: ")
        
        self.v1_status_lab = QLabel(self.status_text[0])
        self.v2_status_lab = QLabel(self.status_text[0])
        self.v3_status_lab = QLabel(self.status_text[0])
        self.v4_status_lab = QLabel(self.status_text[0])
        
        self.pressure_text = QLabel(str(self.pressure))
        P_unit = QLabel(" mBar")
        gb_b = QGroupBox(self)
        gb_b.setTitle("Valve status")
        grid_B.addWidget(lab_1, 0, 0)
        grid_B.addWidget(self.v1_status_lab, 0, 1)
        grid_B.addWidget(lab_2, 1, 0)
        grid_B.addWidget(self.v2_status_lab, 1, 1)
        grid_B.addWidget(lab_3,2,0)
        grid_B.addWidget(self.v3_status_lab, 2, 1)
        grid_B.addWidget(lab_4, 3, 0)
        grid_B.addWidget(self.v4_status_lab, 3, 1)
        grid_B.addWidget(lab_P, 4,0)
        grid_B.addWidget(self.pressure_text, 4,1)
        grid_B.addWidget(P_unit, 4,2)
        gb_b.setLayout(grid_B)
        
        qhbox = QHBoxLayout(self)
        qhbox.addWidget(gb_a)
        qhbox.addWidget(gb_b)
        
        # Connect the signal   
        # 1 - create Worker and Thread inside the Form
        self.obj = Worker()  # no parent!
        self.thread = QThread()  # no parent!   
        # 2 - Connect Worker`s Signals to Form method slots to post data.
        self.obj.datahandle.connect(self.onDataReady)   
        # 3 - Move the Worker object to the Thread object
        self.obj.moveToThread(self.thread)   
        # 5 - Connect Thread started signal to Worker operational slot method
        self.thread.started.connect(self.obj.getdata)    
        # 6 - Start the thread
        self.thread.start()
        
    def button1pressed(self):
        self.sendchar(char="a")
        
    def button2pressed(self):
        self.sendchar(char="b")
        
    def button3pressed(self):
        self.sendchar(char="c")
        
    def button4pressed(self):
        self.sendchar(char="d")
        
    def releasepressed(self):
        self.sendchar(char="y")
    
    def closepressed(self):
        self.sendchar(char="z")
                
    def sendchar(self, char,start_char='<',end_char='>'):
        last_msg=""
        if char != last_msg:
            char_to_send = start_char+char+end_char
            self.arduino.write(char_to_send.encode("utf-8"))
            last_msg=char_to_send
            time.sleep(0.01)
            
    def onDataReady(self,datastring):
        self.v1_status_lab.setText(self.status_text[int(datastring[1])])
        self.v2_status_lab.setText(self.status_text[int(datastring[2])])
        self.v3_status_lab.setText(self.status_text[int(datastring[3])])
        self.v4_status_lab.setText(self.status_text[int(datastring[4])])
        self.pressure_text.setText(datastring.split(']')[0][5:])
        #print(datastring)
        
    def closeEvent(self,event):
        self.thread.quit()
        event.accept()
        self.app.quit()


if __name__ == '__main__':
    #QApplication.setAttribute(Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    global arduino 
    arduino = serial.Serial(port='COM4', baudrate=115200, timeout=None)
    app = QApplication(sys.argv)
    myapp = PCtoArduino(app, arduino)
    myapp.show()
    app.exec_()
    arduino.close()

    
    