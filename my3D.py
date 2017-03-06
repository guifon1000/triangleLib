#kivy
from kivy.app import App
from kivy.uix.floatlayout import FloatLayout
from kivy.clock import Clock

#kivy3
from kivy3 import Renderer, Scene
from kivy3 import PerspectiveCamera

#geometry
from kivy3.extras.geometries import BoxGeometry
from kivy3 import Material, Mesh

import gmsh2kv as my

class My3D(App):
    def _adjust_aspect(self, *args):
        rsize = self.renderer.size
        aspect = rsize[0] / float(rsize[1])
        self.renderer.camera.aspect = aspect

    def rotate_cube(self, *dt):
        self.cube.rotation.x += 2
        self.cube.rotation.y += 2
        self.cube.rotation.z += 2
        self.cube.pos.z -= 0.01
        print self.cube.pos.z
    def build(self):
        layout = FloatLayout()

        # create renderer
        self.renderer = Renderer()

        # create scene
        scene = Scene()

        # create default cube for scene
        geo = my.Msh(file = 'wall.msh',scale = 0.0002)
        mat = Material(
		color = (0.,0.5,0.),
                diffuse = (1,0,0),
                specular = (1,1,0))
        self.cube = Mesh(geometry = geo, material = mat) # default pos == (0,0,0)
        self.cube.pos.z = -2.

        # create camera for scene
        self.camera = PerspectiveCamera(
                fov = 75,    #distance from the screen
                aspect=0,    # "screen" ratio
                near=1,      # nearest rendered point
                far=100    # furthest rendered point
        )
        # start rendering the scene and camera
        scene.add(self.cube)
        self.renderer.render(scene,self.camera)

        # set renderer ratio if its size changes
        # e. g. when added to parent
        self.renderer.bind(size=self._adjust_aspect)

        layout.add_widget(self.renderer)
        Clock.schedule_interval(self.rotate_cube, 0.05)
        return layout

My3D().run()
