# -*- coding: utf-8 -*-
"""
Created on Sun May 12 14:59:55 2019

@author: Emil
"""
import datetime
import pickle


# class for containing datasets
class WongData():
    def __init__(self, data, parametres, initial_conditions):
        # data is the data set indexed firstly by g
        # with entries z,u,q,E,V,A
        
        self.data = data
        self.params = parametres
        self.inits = initial_conditions
        now = datetime.datetime.now()
        self.name = "WongData "+now.strftime("%d-%m-%Y %H-%M")

    # saves class instance to file
    def save(self):
        """save class as self.name.txt"""
        file = open(self.name+'.txt','wb')
        pickle.dump(self,file)
        file.close()

        
    def z(self,g_ind):
        return self.data[g_ind][0]
    
    def u(self,g_ind):
        return self.data[g_ind][1]
    
    def q(self,g_ind):
        return self.data[g_ind][2]
    
    def E(self,g_ind):
        return self.data[g_ind][3]
    
    def V(self,g_ind):
        return self.data[g_ind][4]
    
    def A(self,g_ind):
        return self.data[g_ind][5]
        
    def printparams(self):
        #prints the parametres used for the data set
        g_vec = self.params[0]
        dx = self.params[1]
        Lx = self.params[2]
        dt = self.params[3]
        tmax = self.params[4]
        times = self.params[5]
        xlist = self.params[6]
        
        print("PARAMETRES FOR "+self.name)
        print("g_vec: "+str(min(g_vec))+" to "+str(max(g_vec))+" with "+str(len(g_vec))+" points.")
        print("dx = "+str(dx))
        print("Lx = "+str(Lx))
        print("dt = "+str(dt))
        print("tmax = "+str(tmax))
        print("dx/dt = "+str(dx/dt))
        
    def printinits(self):
        # prints the initial conditions used for the data set
        zx = self.inits[0]
        ux = self.inits[1]
        qx = self.inits[2]
        mass = self.inits[3]
        E0, V0, A0 = self.inits[4], self.inits[5], self.inits[6]

        print("INITIAL CONDITIONS FOR "+self.name)
        print("zx ="+str(zx))        
        print("ux ="+str(ux))
        print("qx ="+str(qx))
        print("mass ="+str(mass))
        print("E0 ="+str(E0))
        print("V0 ="+str(V0))
        print("A0 ="+str(A0))