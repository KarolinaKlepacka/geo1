# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:36:58 2022

@author: kklep
"""
from math import pi, sin, cos, sqrt, degrees, atan, tan, asin, acos, atan2
import numpy as np
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
        
    def xyz2plh(self, X, Y, Z): #output = 'dec_degree'):
        '''
        Funkcja przelicza współrzędne prostokątne na krzywoliniowe.

        Argumenty:
        ----------
        X - współrzędna X punktu | typ: float lub int
        Y - współrzędna X punktu | typ: float lub int
        Z - współrzędna X punktu | typ: float lub int

        Wyniki:
        -------
        phi - szerokość geograficzna punktu | typ: lista złożona z wartości float
        lam - długość geograficzna punktu  | typ: lista złożona z wartości float
        h - wysokość punktu               | typ: float

        Dodatkowy opis:
        ---------------
        Funkcja zamieni wartości współrzędnych prostokątnych średnich na liczby, pod warunkiem, że są to cyfry.
        Funkcja zwraca fi i lambdę w systemie dms

        '''
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev =atan(Z/(r*(1-self.ecc2)))
        lat_next = 0
        epsilon = 0.000001/206265
        while abs(lat_prev - lat_next)>epsilon:
            lat_prev = lat_next
            N    = self.a/ sqrt(1 - self.ecc2 * (sin(lat_prev))**2)
            h  = (r/cos(lat_prev))- N
            lat_next = atan((Z/r) *(((1 - self.ecc2 * N/(N + h))**(-1))))
        lat = lat_next
        lon = atan(Y/X)
        N = self.a/ sqrt(1 - self.ecc2 * (sin(lat))**2);
        return degrees(lat), degrees(lon), h
        return phi, lam, h
    def deg2rad(self, kat):
        katr = kat * pi / 180
        return(katr) 
    def plh2XYZ(self, phi, lam, h):
        """
        Parameters
        ----------
        phi : [float] : szerokość geodezyjna [rad]
        lam : [float] : długość geodezyjna [rad]
        h : [float] : wysokość elipsoidalna [m]

        Returns
        -------
        X: współrzędna prostokątna X punktu | typ: float
            
        Y: współrzędna prostokątna Y punktu | typ: float
            
        Z: współrzędna prostokątna Z punktu | typ: float
        
        Dodatkowy opis:
        ---------------
        Wartości będą zwracane w kolejności: X, Y, Z
        """

        a = self.a
        ecc2 = self.ecc2

        phi = self.deg2rad(phi)
        lam = self.deg2rad(lam)
        
        N =  a / (sqrt(1 - ecc2 * (sin(phi)) ** 2))
        
        X = (N + h) * cos(phi) * cos(lam)
        Y = (N + h) * cos(phi) * sin(lam)
        Z = (N * (1 - ecc2) + h) * sin(phi)
        
        return(X, Y, Z)

    
    def srednia(self, wartosci):
        """
        Funkcja liczy średnią wartość z elementów w liscie
        
        Parameters:
        ----------
        wartosci : [float] : lista wartosci
        
        Returns:
        --------
        srednia : [float] : średnia arytmetyczna elementów z listy 
        
        """
        suma = 0
        ile = 0
        for wartosc in wartosci:
            suma += wartosc
            ile += 1
        srednia = float(suma / ile)
        return(srednia)
    
    def R_neu(self, phi, lam):
        """
        Funkcja, która, przyjmujac współrzedne krzywoliniowe utworzy macierz obrotu 
        potrzebną do przeliczenia współrzędnych do układu współrzędnych neu
    
        INPUT:
        ----------
        fi : [float] : wspołrzędna fi punktu początkowego układu lokalnego
        l : [float] :wspołrzędna l punktu początkowego układu lokalnego
    
        OUTPUT:
        -------
        R : [array of float64] : macierz obrotu
    
        """
        N=[(-np.sin(phi) * np.cos(lam)), (-np.sin(phi) * np.sin(lam)), (np.cos(phi))]
        E=[(-np.sin(lam)), (np.cos(lam)),  (0)]
        U=[( np.cos(phi) * np.cos(lam)), ( np.cos(phi) * np.sin(lam)), (np.sin(phi))]
        R=np.transpose(np.array([N,E,U]))
        return  R, N, E, U
    
    def neu(self, R, v):
        """
        Funckja obliczająca wektor w układzie neu
    
        Parameters:
        -----------
        R : R : [array of float64] : macierz obrotu
        v : [array of float64] : wektor w układzie XYZ
        
        Returns:
        -------
        delta_neu : [array of float64] : współrzedne topocentryczne (North (N), East (E), Up (U))
    
        """
        neu=np.zeros(v.shape)
        for a in range(v.shape[0]):
            for b in range(3):
                for c in range(3):
                    neu[a,c]+=v[a,b]*R[c,b]
                    
        return neu
    def fl2xy(self, phi, lam):
        """
        Funkcja przelicza współrzędne geodezyjne na współrzędne układu 2000 i współrzędne układu 1992.

        Parameters
        ----------
        phi : szerokość geodezyjna [rad]
        lam : długość geodezyjna [rad]
        a   : [float] : dłuższa półoś elipsoidy [m]
        e2  : [float] : mimośrod elipsoidy [niemianowana]
        Returns
        -------
        x2000 : współrzędna w układzie 2000 [m]
        y2000 : współrzędna w układzie 2000 [m].
        x92 : współrzędna w układzie 1992 [m]
        y92 : współrzędna w układzie 1992 [m] 

        """
        # lam0=21*pi/180
        if (lam > 13.5 and lam) < 16.5:
            s = 5
            lam0 = 15
        elif (lam > 16.5 and lam < 19.5):
            s = 6
            lam0 = 18
        elif (lam > 19.5 and lam < 22.5):
            s = 7
            lam0 = 21
        elif (lam > 22.5 and lam < 25.5):
            s = 8
            lam0 = 24
        phi=self.deg2rad(phi)
        lam=self.deg2rad(lam)
        lam0=self.deg2rad(lam0)
        b2=(self.a**2)*(1-self.ecc2)
        ep2=((self.a**2-b2))/(b2)
        t=tan(phi)
        n2=ep2*((cos(phi))**2)
        N =  self.a / (sqrt(1 - self.ecc2 * (sin(phi)) ** 2))
        A0=1-(self.ecc2/4)-((3*(self.ecc2**2))/64)-((5*(self.ecc2**3))/256);
        A2=(3/8)*(self.ecc2+(self.ecc2**2/4)+((15*(self.ecc2**3))/128));
        A4=(15/256)*((self.ecc2**2)+((3*(self.ecc2**3))/4));
        A6=(35*(self.ecc2**3))/3072;
        si=self.a*(A0*(phi)-A2*sin(2*phi)+A4*sin(4*phi)-A6*sin(6*phi));
        dl = lam-lam0
        xgk=si+((dl**2)/2)*N*sin(phi)*cos(phi)*(1+((dl**2)/12)*((cos(phi))**2)*(5-t**2+9*n2+4*(n2**2))+((dl**4)/360)*((cos(phi))**4)*(61-58*(t**2)+t**4+270*n2-330*n2*(t**2)))
        ygk=dl*N*cos(phi)*(1+((dl**2)/6)*((cos(phi))**2)*(1-t**2+n2)+((dl**4)/120)*((cos(phi))**4)*(5-18*(t**2)+t**4+14*n2-58*n2*(t**2)))
        x2000 = xgk * 0.999923
        y2000 = ygk * 0.999923 + ((lam0*180/pi)/3) * 1000000 + 500000
        
        L092 = 19* (pi/180)
        dl1 = lam - L092
        xgk1=si+((dl1**2)/2)*N*sin(phi)*cos(phi)*(1+((dl1**2)/12)*((cos(phi))**2)*(5-t**2+9*n2+4*(n2**2))+((dl1**4)/360)*((cos(phi))**4)*(61-58*(t**2)+t**4+270*n2-330*n2*(t**2)))
        ygk1=dl1*N*cos(phi)*(1+((dl1**2)/6)*((cos(phi))**2)*(1-t**2+n2)+((dl1**4)/120)*((cos(phi))**4)*(5-18*(t**2)+t**4+14*n2-58*n2*(t**2)))
        x92 = xgk1 *0.9993 - 5300000
        y92 = ygk1 * 0.9993 + 500000
        return x2000, y2000, x92, y92
    def azymut(self, N, E):
        """   
        Funkcja wyznacza kąt azymutu na podstawie współrzędnych topocentrycznych
    
        Parameters
        -------
        N  : [float] : wpółrzedna topocentryczna N (north) [m]
        E  : [float] : wpółrzedna topocentryczna E (east) [m]
   
        Returns
        -------
        Az : [float] : azymut [rad]
    
        """  
        Az = atan2(E, N)
    
        return Az
    def kat_elewacji(self, N, E, U):
        """   
        Funkcja wyznacza kąt elewacji na podstawie współrzędnych topocentrycznych
    
        Parameters
        -------
        N  : [float] : wpółrzedna topocentryczna N (north) [m]
        E  : [float] : wpółrzedna topocentryczna E (east) [m]
        U  : [float] : wpółrzedna topocentryczna U (up) [m] 
   
        Returns
        -------
        elewacja : [float] : kąt elewacji [rad]
    
        """ 
        elewacja = acos(U/sqrt(N**2 + E**2 + U**2))
        return elewacja
    def dl2D(self,X_sr, Y_sr, X, Y):
        """   
        Funkcja wyznacza odległość na płaszczyźnie
        na podstawie współrzędnych płaskich prostokątnych
    
        Parameters
        -------
        X_sr  : [float] : średnie X wszystkich punktów [m]
        Y_sr : [float] : średnie Y wszystkich punktów[m]
        X  : [float] : współrzędna X punktu  [m]
        Y  : [float] : współrzędna Y punktu [m]
   
        Returns
        -------
        d : [float] : odległość na płaszczyźnie [m]
    
        """
        d2d = []
        for X, Y in zip(X,Y):
            d = sqrt((X_sr - X)**2 + (Y_sr - Y)**2)
            d2d.append(d)
        return d2d
    def dl3D(self,N_sr, E_sr, U_sr, N, E, U):
            
        '''
        Funkcja liczy długość wektora przyjmując elementy z listy jako współrzędne końcowe wektora, a współrzędne średnie jako współrzędne początkowe.

        Argumenty:
        ----------
        N_sr - współrzędna N średnia                     | typ: float lub int                       E_sr - współrzędna E średnia                     | typ: float lub int
        U_sr - współrzędna U średnia                     | typ: float lub int
        dN - lista współrzędnych N dla kolejnych punktów | typ: lista wartości float lub int
        dE - lista współrzędnych E dla kolejnych punktów | typ: lista wartości float lub int
        dU - lista współrzędnych U dla kolejnych punktów | typ: lista wartości float lub int

        Wyniki:
        -------
        Dlugosc - lista odległośći punktów od X, Y, Z średnich | typ: lista wartości float
    
        Dodatkowy opis:
        ---------------
    
        ---
        '''
        Dlugosc = []
    
        for N, E, U in zip(N, E, U):
            D = sqrt((N_sr - N) ** 2 + (E_sr - E) ** 2 + (U_sr - U) ** 2)
            Dlugosc.append(D)
        return(Dlugosc)
if __name__ == "__main__":
    #utworzenie obiektu
    geo = Transformacje(model = "grs80")
    plik = "wsp_inp.txt"
    tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)
    
    odczyt = np.ones((12,3))
    for i,n in enumerate(tablica):
        phi,lam,h = geo.xyz2plh(tablica[i,0],tablica[i,1],tablica[i,2])
        odczyt[i,:] = [phi,lam,h]
    print(odczyt)
    
    o = np.ones((12,3))
    for i,n in enumerate(odczyt):
        X, Y, Z = geo.plh2XYZ(odczyt[i,0],odczyt[i,1],odczyt[i,2])
        o[i,:] = [X,Y,Z]
    print(o)

    fi_sr = geo.srednia(odczyt[:,0])
    lam_sr = geo.srednia(odczyt[:,1])
    X_sr = geo.srednia(o[:,0])
    Y_sr = geo.srednia(o[:,1])
    Z_sr = geo.srednia(o[:,2])
    [fi_sr,lam_sr,h_sr]=geo.xyz2plh(X_sr, Y_sr, Z_sr)
    R, N, E, U=geo.R_neu(fi_sr,lam_sr)
    i=0
    v=np.array(np.zeros((o.shape[0],3)))
    for i in range(0,o.shape[0]):
        v[i,0]=X_sr-o[i,0]
        v[i,1]=Y_sr-o[i,1]
        v[i,2]=Z_sr-o[i,2]
        i=i+1
        
    print(N, E, U)
    neu=geo.neu(R, v)
    print(neu)
    az = np.ones((12,1))
    for i, n in enumerate(neu):
        Az = geo.azymut(neu[i,0], neu[i,1])
        az[i,:] = [Az]

    print(az)
    elewacja = np.ones((12,1))
    for i, n in enumerate(neu):
        el = geo.kat_elewacji(neu[i,0], neu[i,1],neu[i,2])
        elewacja[i,:] = [el]

    print(elewacja)
    
    
    d2d = geo.dl2D(X_sr, Y_sr, o[:,0], o[:,1])
    print(d2d)
    arr1 = np.array([d2d])
    odl2D = arr1.T
    print(odl2D)
    
    N_sr = geo.srednia(neu[:,0])  
    E_sr = geo.srednia(neu[:,1]) 
    U_sr = geo.srednia(neu[:,2]) 
    Dlugosc = geo.dl3D(N_sr, E_sr, U_sr, neu[:,0], neu[:,1], neu[:,2])
    arr = np.array([Dlugosc])
    odleglosc=arr.T
    print(odleglosc)
    r = np.ones((12,4))
    for i, n in enumerate(odczyt):
        x2000, y2000, x92, y92 = geo.fl2xy(odczyt[i,0], odczyt[i,1])
        r[i,:]=[x2000, y2000, x92, y92]
    print(r)
    Dane = np.hstack((odczyt,o,neu,r,az,elewacja,odl2D,odleglosc))
    print(Dane)
    np.savetxt("wsp_out.txt", Dane, delimiter='  ', fmt = ['%10.8f', '%10.8f', '%10.5f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.3f', '%10.3f', '%10.3f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], header = 'Konversja współrzednych geodezyjnych \ Karolina Klepacka', comments= 'fi[st]        lambda[st]      h[m]          X[m]          Y[m]                 Z[m]              N          E          U            X 2000 [m]       Y 2000[m]    X 1992[m]    Y 1992[m]  Azymut[rad] Kąt elewacji[rad] 2D        3D \n')
