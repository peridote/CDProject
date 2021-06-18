from SceneGenerator import *
import math 
import random

def scale(vec, s):
    vec[0] *= s
    vec[1] *= s
    vec[2] *= s
    return vec

s = 1

scene = generateScene('ClothScene', camPosition=[0,15,30], camLookat=[0,0,0])
addParameters(scene, h=0.005, maxIterVel=5, contactTolerance=0.05)
friction = 0.1
restitution = 0.6

addRigidBody(scene, '../models/cube.obj', 2, translation = [0, -10, 0], coScale=[100, 1, 100], scale=[100, 1, 100], dynamic=0)
"""
addRigidBody(scene, '../models/torus.obj', 4, coScale=[2, 1, 2], scale=[2, 2, 2], translation=[0,5,0], 
             friction=friction, rest = restitution, dynamic=0)
"""
sp = []
for i in range(51):
    sp.append(51*i)
    sp.append(51*i+50)
    sp.append(i)
    sp.append(2550+i)
addTriangleModel(scene, '../models/plane_50x50.obj', translation=[0, 14, 0], 
                 scale=[10,10,10], friction=friction, rest = restitution, staticParticles=sp)
#addRigidBody(scene, '../models/sphere.obj', 1, translation=[1, 16, 1], coScale=[1, 1, 1], density = 400, scale=[1, 1, 1], rest = 0.1, dynamic=1, isdraw = False)

#jointScale = [0.1,0.1,10]
#joint = addRigidBody(scene, '../models/cube.obj', 0, coScale=jointScale, scale=jointScale, translation=[5, 10, 0], dynamic=1)
#addRigidBodyParticleBallJoint(scene, joint, 50)
#joint = addRigidBody(scene, '../models/cube.obj', 2, coScale=jointScale, scale=jointScale, translation=[5, 10, -5], dynamic=1)
#addRigidBodyParticleBallJoint(scene, joint, 2600)
#addRigidBodyParticleBallJoint(scene, joint, 1325)

#slider = addRigidBody(scene, '../models/cube.obj', 0, coScale=[0.1,0.1,0.1], scale=[0.1,0.1,0.1], translation=[0, 10, -5], dynamic=0)
#addTargetVelocityMotorSliderJoint(scene, slider, joint, axis=[1,0,0], target= 10.0)

armadilloScale = [0.5,0.5,0.5]
bunnyScale = [1,1,1]
armadilloScale = scale(armadilloScale, s)
bunnyScale = scale(bunnyScale, s)
density = 400
friction = 0.1
restitution = 0.6
t = [0, 13, 0]

#addRigidBody(scene, '../models/armadillo.obj', 5, coFile='', coScale=armadilloScale, translation=[-2, 13, 0], scale=armadilloScale, dynamic=1)
#addRigidBody(scene, '../models/bunny_10k.obj', 5, coFile='../sdf/bunny_10k.csdf', coScale=bunnyScale, translation=t, scale=bunnyScale, dynamic=1, density=density)
for i in range(5):
    for j in range(5):
        addRigidBody(scene, '../models/sphere.obj', 1, translation=[-3.0+j, 15.5+i, 3.0-i], coScale=[0.3, 0.3, 0.3], density = 400, scale=[0.3, 0.3, 0.3], rest = 0.1, dynamic=1)
#addRigidBody(scene, '../models/sphere.obj', 1, translation=[2, 11, 0], coScale=[0.5, 0.5, 0.5], scale=[0.5, 0.5, 0.5], dynamic=1)

addRigidBody(scene, '../models/hollowhalfsphere.obj', 8, translation=[0, 8, 0], coScale=[4.5, 4.5, 4.5], density = 400, scale=[4.5, 4.5, 4.5], rest = 0.6, dynamic=0)

wallScale = [12, 10, 0.5]
addRigidBody(scene, '../models/cube.obj', 2, translation = [0, -2, -1], coScale=wallScale, scale=wallScale, dynamic=0)
addRigidBody(scene, '../models/cube.obj', 2, translation = [0, -2, 1.25], coScale=wallScale, scale=wallScale, dynamic=0, isdraw= True)
wallScale = [0.5, 10, 2]
addRigidBody(scene, '../models/cube.obj', 2, translation = [-6, -2, 0], coScale=wallScale, scale=wallScale, dynamic=0, isdraw= True)
addRigidBody(scene, '../models/cube.obj', 2, translation = [5.85, -2, 0], coScale=wallScale, scale=wallScale, dynamic=0, isdraw= True)

# piles
restitution = 0.6
num_piles_x = 7
num_piles_y = 7
dx_piles = 5.0
dy_piles = 5.0
startx_piles = -0.5 * (num_piles_x - 1.0)*dx_piles
starty_piles = -0.5 * (num_piles_y - 1.0)*dy_piles
pileAxis = [1, 0 ,0]
pileAngle = math.radians(90);
pileScale = [0.125,2.5,0.125]
pileScale = scale(pileScale, s)

current_y = starty_piles
for i in range(0,num_piles_y):
    current_x = startx_piles
    if i%2 == 1:
        t = [-3.9-0.65, 2-1.01*i, 0]
        addRigidBody(scene, '../models/cylinder.obj', 3, coScale=pileScale, scale=pileScale, axis = pileAxis, angle = pileAngle,
                         translation=t, dynamic=0, rest=restitution)
    for j in range(0,num_piles_x):
        if i%2 == 0:
            t = [-3.9+1.3*j, 2-1.01*i, 0]
            addRigidBody(scene, '../models/cylinder.obj', 3, coScale=pileScale, scale=pileScale, axis = pileAxis, angle = pileAngle,
                         translation=t, dynamic=0, rest=restitution)
        else:
            t = [-3.9+0.65+1.3*j, 2-1.01*i, 0]
            addRigidBody(scene, '../models/cylinder.obj', 3, coScale=pileScale, scale=pileScale, axis = pileAxis, angle = pileAngle,
                         translation=t, dynamic=0, rest=restitution)

writeScene(scene, 'ClothScene.json')
