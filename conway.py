import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from random import randint
from random import random as rand
import kivy
from kivy.app import App
from kivy.uix.gridlayout import GridLayout
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.button import Button
from kivy.uix.label import Label
from kivy.uix.checkbox import CheckBox
from kivy.graphics import Color, Rectangle
from kivy.clock import Clock
from kivy.uix.widget import Widget
from kivy.properties import NumericProperty
from functools import partial
from kivy.config import Config
Config.set('graphics', 'width', '900')
Config.set('graphics', 'height', '900')
N= 120
Nc = 2


class World(GridLayout):
    def __init__(self,rows,cols):
        super(World, self).__init__(rows=rows,cols=cols)
        self.rows = rows
        self.cols = cols
        self.index=None
        self.cells=[Cell()]*rows*cols
 
    def update(self,*largs):
        for c in self.cells:
            c.countAlive(self.cells)
            c.nextState=0
            #if c.state==1:c.aliveNeig+=1
            if c.aliveNeig==3:
                c.nextState=1
            elif c.aliveNeig==2:
                c.nextState=c.state
            elif c.aliveNeig<2 or c.aliveNeig>3:
                c.nextState=0

        for c in self.cells:
            c.setAlive()
            c.canvas.clear()
            with c.canvas:
                  
                if c.state==1:
                    self.c0=1
                    self.c1=1
                    self.c2=1
                else:
                    self.c0=0
                    self.c1=0
                    self.c2=0
                Color(self.c0,self.c1,self.c2)
                Rectangle(pos=c.center, size=(50, 50))
            #print c.color
     

class ConwayApp(App):
    def build(self):
        cellZone = World(N,N)
        cellZone.size=[100,100]
        #cellZone.pos=[1,1]
        index=-1
        for i in range(cellZone.rows):
            for j in range(cellZone.cols):
                cell = Cell(i,j)
                index+=1
                cell.index = index
                cell.neighbours()
                cell.state=randint(0,1)
                cell.nextState=cell.state
                cellZone.add_widget(cell)
                cellZone.cells[index]=cell
        cellZone.update()
        for c in cellZone.cells:
            c.countAlive(cellZone.cells)

        dt=0.1
        Clock.schedule_interval(cellZone.update,dt)
        return cellZone




class Cell(Widget):
    c0 = NumericProperty()
    c1 = NumericProperty()
    c2 = NumericProperty()
    def __init__(self,i=0,j=0):
        self.color=1.0,1.0,1.0
        self.irow=i
        self.icol=j
        self.state=0   
        self.nextState=0
        self.iNorth='*'
        self.iEast='*'
        self.iSouth='*'
        self.iWest='*'
        self.iNorthEast='*'
        self.iNorthWest='*'
        self.iSouthEast='*'
        self.iSouthWest='*'
        self.index=-1
        self.aliveNeig=0
        self.color = Color(rand(), 1, 1)
        super(Cell, self).__init__()
    
    def neighbours(self):
        if self.irow-1 >= 0 :
            self.iNorth = self.index-N
        else :
            self.iNorth = '*'
        if self.irow+1 <= N-1 :
            self.iSouth = self.index+N
        else:
            self.iSouth = '*'
        if self.icol-1 >= 0 :
            self.iWest =self.index-1
        else :
            self.iWest='*'
        if self.icol+1 <= N-1 :
            self.iEast = self.index+1
        else:
            self.iEast='*'
        if self.irow-1 >= 0 and self.icol+1 <= N-1:
            self.iNorthEast=self.index-N+1
        else:
            self.iNorthEast='*' 
        if self.irow-1 >= 0 and self.icol-1 >= 0 :
            self.iNorthWest=self.index-N-1
        else:
            self.iNorthWest='*'
        if self.irow+1 <= N-1 and self.icol+1 <= N-1:
            self.iSouthEast=self.index+N+1
        else:
            self.iSouthEast='*'
        if self.irow+1 <= N-1 and self.icol-1 >= 0:
            self.iSouthWest=self.index+N-1
        else:
            self.iSouthWest='*'

    def countAlive(self,tab):
        self.aliveNeig=0
        neig = [self.iNorth,self.iSouth,self.iWest,self.iEast,\
                self.iNorthEast,self.iNorthWest,self.iSouthEast,self.iSouthWest]
        for n in neig:
            if n!='*' and tab[n].state==1:
                self.aliveNeig+=1
        
        
    def __str__(self):
        s = str(self.iNorthWest)+' '+str(self.iNorth)+' '+str(self.iNorthEast)+'\n'
        s+= str(self.iWest)+' '+str(self.index)+' '+str(self.iEast)+'\n'
        s+= str(self.iSouthWest)+' '+str(self.iSouth)+' '+str(self.iSouthEast)
        return s

        
    def setAlive(self):
        self.state=0
        if self.nextState==1:
            self.state=1
        elif self.nextState==0:
            self.state=0
        self.nextState=0

if __name__=='__main__':
    ConwayApp().run()

