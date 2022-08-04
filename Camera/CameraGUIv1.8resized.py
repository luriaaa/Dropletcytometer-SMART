from PyQt5 import QtCore, QtGui, QtWidgets
import cv2
import sys
from PyQt5.QtWidgets import QWidget, QLabel, QApplication, QLineEdit, QPushButton
from PyQt5.QtCore import QThread, Qt, pyqtSignal, pyqtSlot,QTimer
from PyQt5.QtGui import QImage, QPixmap, QColor
# https://github.com/elerac/EasyPySpin just to make life easier with pyspin
import EasyPySpin

import numpy as np
import time

class CamGUI(QWidget):
    def __init__(self):
        super().__init__()
        # your "global" parameters
        self.exposure = 5000
        self.gain = 30
        self.gamma = 0.25
        self.fps = 200
        self.laser_trace_on = False # tracker to determine whether show the circle for laser spot
        self.max_exp = 100000
        self.min_exp = 1
        self.max_gain = 47
        self.min_gain = 0
        self.max_gamma = 1
        self.min_gamma = 0.25
        self.max_fps = 200
        self.min_fps = 1
        self.laser_tracer_coord_x = 363
        self.laser_tracer_coord_y = 495
        self.laser_tracer_dia = 36
        self.display_width = 1152
        self.display_height = 864
        self.setWindowTitle("Blackfly U3-S16C")
        self.resize(1440, 900)
        self.video = QLabel(self)
        self.video.setGeometry(QtCore.QRect(10, 10, self.display_width, self.display_height))

        # start the video.capture, put up one frame
        self.init_camera()
              
        Qlab_w = 51
        Qlab_h = 31
        Qline_w = 61
        Qline_h = 31
        
        rangetext_w = 71
        rangetext_h = 21
        
        Qbutton_w = 101
        Qbutton_h = 31
        
        left_align = 1180
        top_align = 10
        spacing_v = 20
        spacing_h = 20
        
        # GUI widgets
        self.exposure_label = QLabel("Exposure", self)
        self.exposure_label.setGeometry(QtCore.QRect(left_align, top_align, Qlab_w, Qlab_h))
        self.exp_val = QLineEdit(str(self.exposure), self)
        self.exp_val.setGeometry(QtCore.QRect(left_align+Qlab_w+spacing_h, top_align, Qline_w, Qline_h))
        self.exp_range = QLabel("us (1~100000)", self)
        self.exp_range.setGeometry(QtCore.QRect(left_align+Qlab_w+Qline_w+spacing_h*2, top_align, rangetext_w+spacing_h, rangetext_h))

        
        self.gain_label = QLabel("Gain", self)
        self.gain_label.setGeometry(QtCore.QRect(left_align, top_align+spacing_v+Qlab_h, Qlab_w, Qlab_h))
        self.gain_val = QLineEdit(str(self.gain), self)
        self.gain_val.setGeometry(QtCore.QRect(left_align+Qlab_w+spacing_h, top_align+spacing_v+Qlab_h, Qline_w, Qline_h))
        self.gain_range = QLabel("(0~48)", self)
        self.gain_range.setGeometry(QtCore.QRect(left_align+Qlab_w+Qline_w+spacing_h*2, top_align+spacing_v+Qlab_h, rangetext_w, rangetext_h))
        
        self.gamma_label = QLabel("Gamma", self)
        self.gamma_label.setGeometry(QtCore.QRect(left_align, top_align+spacing_v*2+Qlab_h*2, Qlab_w, Qlab_h))
        self.gamma_val = QLineEdit(str(self.gamma), self)
        self.gamma_val.setGeometry(QtCore.QRect(left_align+Qlab_w+spacing_h, top_align+spacing_v*2+Qlab_h*2, Qline_w, Qline_h))
        self.gamma_range = QLabel("(0.25~1)", self)
        self.gamma_range.setGeometry(QtCore.QRect(left_align+Qlab_w+Qline_w+spacing_h*2, top_align+spacing_v*2+Qlab_h*2, rangetext_w, rangetext_h))
        
        self.fps_label = QLabel("fps", self)
        self.fps_label.setGeometry(QtCore.QRect(left_align, top_align+spacing_v*3+Qlab_h*3, Qlab_w, Qlab_h))
        self.fps_val = QLineEdit(str(self.fps), self)
        self.fps_val.setGeometry(QtCore.QRect(left_align+Qlab_w+spacing_h, top_align+spacing_v*3+Qlab_h*3, Qline_w, Qline_h))
        self.fps_range = QLabel("(1~225)", self)
        self.fps_range.setGeometry(QtCore.QRect(left_align+Qlab_w+Qline_w+spacing_h*2, top_align+spacing_v*3+Qlab_h*3, rangetext_w, rangetext_h))
        
        self.preset_bf = QPushButton("Bright field", self)
        self.preset_bf.setGeometry(QtCore.QRect(left_align, top_align+spacing_v*4+Qlab_h*4, Qbutton_w, Qbutton_h))
        self.preset_df = QPushButton("Dark field", self)
        self.preset_df.setGeometry(QtCore.QRect(left_align+spacing_h+Qbutton_w, top_align+spacing_v*4+Qlab_h*4, Qbutton_w, Qbutton_h))
        self.set_cam = QPushButton("Set CAM Param", self)
        self.set_cam.setGeometry(QtCore.QRect(left_align, top_align+spacing_v*5+Qlab_h*5, Qbutton_w*2+spacing_h, Qbutton_h))

        
        self.tracer_coord_x_label = QLabel("Laser X", self)
        self.tracer_coord_x_label.setGeometry(QtCore.QRect(left_align, top_align+spacing_v*6+Qlab_h*6, Qlab_w, Qlab_h))
        self.tracer_coord_x = QLineEdit(str(self.laser_tracer_coord_x), self)
        self.tracer_coord_x.setGeometry(QtCore.QRect(left_align+Qlab_w+spacing_h, top_align+spacing_v*6+Qlab_h*6, Qline_w, Qline_h))
        self.tracer_coord_y_label = QLabel("Laser Y", self)
        self.tracer_coord_y_label.setGeometry(QtCore.QRect(left_align, top_align+spacing_v*7+Qlab_h*7, Qlab_w, Qlab_h))
        self.tracer_coord_y = QLineEdit(str(self.laser_tracer_coord_y), self)
        self.tracer_coord_y.setGeometry(QtCore.QRect(left_align+Qlab_w+spacing_h, top_align+spacing_v*7+Qlab_h*7, Qline_w, Qline_h))
        self.tracer_dia_label = QLabel("Diameter", self)
        self.tracer_dia_label.setGeometry(QtCore.QRect(left_align, top_align+spacing_v*8+Qlab_h*8, Qlab_w, Qlab_h))
        self.tracer_dia = QLineEdit("36", self)
        self.tracer_dia.setGeometry(QtCore.QRect(left_align+Qlab_w+spacing_h, top_align+spacing_v*8+Qlab_h*8, Qline_w, Qline_h))
        self.laser_tracer = QPushButton("Laser Tracer", self)
        self.laser_tracer.setGeometry(QtCore.QRect(left_align, top_align+spacing_v*9+Qlab_h*9, Qbutton_w*2+spacing_h, Qbutton_h))



        # connect button to functions
        self.set_cam.clicked.connect(self.set_cam_params)
        self.laser_tracer.clicked.connect(self.laser_trace_tracker)
        self.preset_bf.clicked.connect(self.preset_brightfield)
        self.preset_df.clicked.connect(self.preset_darkfield)

        # timer to update each frame
        self.timer = QTimer()
        self.timer.timeout.connect(self.nextFrameSlot)
        self.timer.start(40) # 25 fps, useful range

    def init_camera(self):
        self.cap = EasyPySpin.VideoCapture(0)
        if not self.cap.isOpened():
            print("Camera can't open\nexit")
            return -1

        # Set the camera parameters (use Easypyspin wrapper)
        self.cap.set(cv2.CAP_PROP_EXPOSURE, -1)  # -1 sets exposure_time to auto
        self.cap.set(cv2.CAP_PROP_GAIN, -1)  # -1 sets gain to auto
        self.cap.set(cv2.CAP_PROP_GAMMA, -1)  # -1 sets gain to auto
        self.cap.set(cv2.CAP_PROP_FPS, 50)  # FPS explicitly set to 50
        ret, cv_img = self.cap.read() # read once
        qt_img = self.convert_cv_qt(cv_img)
        self.video.setPixmap(qt_img)
        #self.cap.release()

    def nextFrameSlot(self):
        # update the camera setup
        self.cap.set(cv2.CAP_PROP_EXPOSURE, self.exposure)  # -1 sets exposure_time to auto
        self.cap.set(cv2.CAP_PROP_GAIN, self.gain)  # -1 sets exposure_time to auto
        self.cap.set(cv2.CAP_PROP_GAMMA, self.gamma)  # -1 sets exposure_time to auto
        self.cap.set(cv2.CAP_PROP_FPS, self.fps)  # -1 sets exposure_time to auto
        rval, frame = self.cap.read()
        qt_frame = self.convert_cv_qt(frame)
        self.video.setPixmap(qt_frame)

    def convert_cv_qt(self, cv_img):
        """Convert from an opencv image to QPixmap"""

        if self.laser_trace_on :
            self.laser_tracer_coord_x = int(self.tracer_coord_x.text())
            self.laser_tracer_coord_y = int(self.tracer_coord_y.text())
            if (isinstance(self.laser_tracer_coord_x , int) and isinstance(self.laser_tracer_coord_y, int)) and (0<self.laser_tracer_coord_x < 1152 and 0<self.laser_tracer_coord_y<864): 
                    self.laser_tracer_dia = int(self.tracer_dia.text())
                    cv_img = cv2.circle(cv_img, (self.laser_tracer_coord_x,self.laser_tracer_coord_y), 1,
                                        (0, 0, 255), 2)
                    cv_img = cv2.circle(cv_img, (self.laser_tracer_coord_x,self.laser_tracer_coord_y), self.laser_tracer_dia,
                                        (0, 0, 255), 2)

        rgb_image = cv2.cvtColor(cv_img, cv2.COLOR_BAYER_GR2BGR)
        rgb_image = cv2.resize(rgb_image, (self.display_width,self.display_height))
        h, w, ch = rgb_image.shape
        bytes_per_line = ch * w
##        self.display_width = 1152
##        self.display_height = 864
        convert_to_Qt_format = QtGui.QImage(rgb_image.data, w, h, bytes_per_line, QtGui.QImage.Format_RGB888)
        p = convert_to_Qt_format.scaled(self.display_width, self.display_height, Qt.KeepAspectRatio)
        return QPixmap.fromImage(p)

    def set_cam_params(self):
        if self.exposure != int(self.exp_val.text()):
            if self.exposure > self.max_exp:
                self.exposure = self.max_exp
            elif self.exposure < self.min_exp:
                self.exposure = self.min_exp
            else:
                self.exposure = int(self.exp_val.text())
            self.cap.set(cv2.CAP_PROP_EXPOSURE, self.exposure)

        if self.gain != int(self.gain_val.text()):
            if self.gain > self.max_gain:
                self.gain = self.max_gain
            elif self.gain < self.min_gain:
                self.gain = self.min_gain
            else:
                self.gain = int(self.gain_val.text())
            self.cap.set(cv2.CAP_PROP_GAIN, self.gain)

        if self.gamma != float(self.gamma_val.text()):
            if self.gamma > self.max_gamma:
                self.gamma = self.max_gamma
            elif self.gamma < self.min_gamma:
                self.gamma = self.min_gamma
            else:
                self.gamma = float(self.gamma_val.text())
            self.cap.set(cv2.CAP_PROP_GAMMA, self.gamma)

        if self.fps != int(self.fps_val.text()):
            if self.fps > self.max_fps:
                self.fps = self.max_fps
            elif self.fps < self.min_fps:
                self.fps = self.min_fps
            elif self.exposure * self.fps > 1000000:
                self.fps = -1
            else:
                self.fps = int(self.fps_val.text())
            self.cap.set(cv2.CAP_PROP_FPS, self.fps)

        print("exposure: ", self.exposure)
        print("gain: ", self.gain)
        print("gamma: ", self.gamma)
        print("fps: ", self.fps)
        print("--------------------")

    def preset_brightfield(self):

        self.exposure = 5000
        self.cap.set(cv2.CAP_PROP_EXPOSURE, self.exposure)

        self.gain = 30
        self.cap.set(cv2.CAP_PROP_GAIN, self.gain)

        self.gamma = 0.25
        self.cap.set(cv2.CAP_PROP_GAMMA, self.gamma)

        self.fps = 50
        self.cap.set(cv2.CAP_PROP_FPS, self.fps)

        print("exposure: ", self.exposure)
        print("gain: ", self.gain)
        print("gamma: ", self.gamma)
        print("fps: ", self.fps)
        print("--------------------")

    def preset_darkfield(self):

        self.exposure = 10
        self.cap.set(cv2.CAP_PROP_EXPOSURE, self.exposure)

        self.gain = 4
        self.cap.set(cv2.CAP_PROP_GAIN, self.gain)

        self.gamma = 0.25
        self.cap.set(cv2.CAP_PROP_GAMMA, self.gamma)

        self.fps = 50
        self.cap.set(cv2.CAP_PROP_FPS, self.fps)

        print("exposure: ", self.exposure)
        print("gain: ", self.gain)
        print("gamma: ", self.gamma)
        print("fps: ", self.fps)
        print("--------------------")

    def laser_trace_tracker(self):
        self.laser_trace_on = not self.laser_trace_on


    # Override the close window event to shut down camera properly
    def closeEvent(self,event):
        self.cap.release()
        event.accept()



if __name__ == "__main__":

    app = QApplication(sys.argv)
    MainWindow = CamGUI()
    MainWindow.show()
    sys.exit(app.exec_())
