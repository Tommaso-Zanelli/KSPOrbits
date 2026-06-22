#!/usr/bin/env python3
import time

# Timer object
class timer:
    def __init__(self):
        self.t1 = 0.0
        self.t2 = 0.0
        self.dt = 0.0
    def tic(self):
        self.t1 = time.time()
    def toc(self, outStr = None):
        self.t2 = time.time()
        self.dt = self.t2 - self.t1
        if outStr != None:
            print(outStr, self.dt)
        return self.dt

# Constants
tiny = 5.56268464626801e-309
HUGE = 1.7976931348623157e+308
degToRad = 0.01745329251994329547
radToDeg = 57.29577951308232286