import json
#kivy
from kivy.app import App
from kivy.uix.floatlayout import FloatLayout
from kivy.clock import Clock
from kivy.graphics import Color, Rectangle,Ellipse,Line
from kivy.uix.widget import Widget
from kivy.uix.button import Button
from kivy.uix.label import Label
from kivy.core.text import Text
from kivy.uix.behaviors import ButtonBehavior


class kPoint(ButtonBehavior, Widget):
    def __init__(self,**kwargs):
        super(kPoint,self).__init__(**kwargs)
        self.name = kwargs['name']
        self.selected = False
        self.color = (0.5,0,0)
        self.lmax = 30
        self.bind(on_press = self.set_selected)
        


    def return_position(self,*args):
        return self.name

    def set_selected(self,*args):
        self.selected = True
        self.color = (0,0,0.5)


class kLabel(Widget):
    def __init__(self,**kwargs):
        super(kLabel,self).__init__(**kwargs)


class kLine(Widget):
    def __init__(self,**kwargs):
        super(kLine,self).__init__()
        p0 = kwargs['p0']
        p1 = kwargs['p1']
        self.points = p0.pos[0] , p0.pos[1] ,\
                      p1.pos[0] , p1.pos[1] 

class AddWallButton(Button):
    def __init__(self, **kwargs):
        super(AddWallButton,self).__init__(**kwargs)
        self.name = 'add_wall'

        


class ReadHouseApp(App):
    def build(self):
        layout = FloatLayout(size=(10,10))
        awb = AddWallButton(size_hint=(0.1,0.1),text='Add Wall')
        awb.bind(on_press = self.add_wall)
        layout.add_widget(awb)
        self.pts = {}
        scale = 100
        self.newWall = False
        for p in d['points']:
            pp = d['points'][p]
            self.pts[p] = kPoint(name= str(p), posi = (pp[0],pp[1]), lmax = 30)
        for p in self.pts:
            layout.add_widget(self.pts[p])
        self.wls = d['internal_walls']
        self.lay=layout
        Clock.schedule_interval(self.update,0.01)
        return layout

    def update(self,*args):
        wl = []
        if self.newWall:
            for p in self.pts:
                if (self.pts[p].selected) and (p not in wl):
                    wl.append(p)
            if len(wl)>=2:
                print 'new wall '+wl[0]+wl[1]
                self.wls[wl[0]+wl[1]]=[wl[0] , wl[1]]
                for p in self.pts:
                    self.pts[p].selected = False
                    self.pts[p].color =(0.5,0,0) 
                    self.newWall = False

        for x in self.lay.children:
            if type(x).__name__ == 'kLine':self.lay.remove_widget(x)

        for l in self.wls.keys():
            p0 =  self.pts[self.wls[l][0]]
            p1 =  self.pts[self.wls[l][1]]
            lk = kLine(p0 = p0 , p1 = p1)
            self.lay.add_widget(lk)            

    def add_wall(self,*args):
        print " adding a wall between "
        self.newWall = True

    def save_file(self,*args):
        print 'saving json'


class House(dict):
    def __init__(self, **kwargs):
        if kwargs.has_key('file'):
            sre = {}
            sre = json.load(open(kwargs['file'],'r'))
            super(House,self).__init__(sre)

if __name__ == '__main__':    
    d= House(file='romarine.json')
    ReadHouseApp().run()

