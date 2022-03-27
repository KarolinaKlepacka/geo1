# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:36:58 2022

@author: kklep
"""
from math import sin, cos, sqrt, degrees, atan
class Transformacje:
    def __init__(self, model: str = "wgs84"):
        #https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = (2 * self.flattening - self.flattening ** 2)
        print(model,self.b)
        
    def xyz2plh(self, X, Y, Z):
        
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev =atan(Z/(r*(1-self.ecc2)))
        lat_next = 0
        epsilon = 0.000001/206265
        while abs(lat_prev - lat_next)>epsilon:
            lat_prev = lat_next
            N    = self.a/ sqrt(1 - self.ecc2 * (sin(lat_prev))**2)
            hel  = (r/cos(lat_prev))- N
            lat_next = atan((Z/r) *(((1 - self.ecc2 * N/(N + hel))**(-1))))
        lat = lat_next
        lon = atan(Y/X)
        N = self.a/ sqrt(1 - self.ecc2 * (sin(lat))**2);

        return degrees(lat), degrees(lon), hel
        return phi, lam, hel


if __name__ == "__main__":
    #utworzenie obiektu
    geo = Transformacje(model = "grs80")
    #dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170 
    phi, lam, hel = geo.xyz2plh(X, Y, Z)
    print(phi, lam, hel)