from SceneGenerator import *
import math 
import random

def scale(vec, s):
    vec[0] *= s
    vec[1] *= s
    vec[2] *= s
    return vec

s = 1

scene = generateScene('MyScene', camPosition=[0,20,60], camLookat=[10,0,0])
addParameters(scene, h=0.005, maxIterVel=5, contactTolerance=0.02)

# floor
floorScale=[500, 1, 500]
floorScale = scale(floorScale, s)
floorT = [0,-0.5,0]
floorT = scale(floorT, s)
addRigidBody(scene, '../models/cube.obj', 2, coScale=floorScale, 
             scale=floorScale, translation=floorT, 
             dynamic=0, rest= 0.5)

# zenga
blknum = 2
height = 10
blkScale = [1, 1 , 3]
blkScale = scale(blkScale, s)
coblkScale = add_vector(blkScale, [0.01, 0.01, 0.01])
blkScaleT = [3, 1, 1]
blkScaleT = scale(blkScaleT, s)
coblkScaleT = add_vector(blkScaleT, [0.01, 0.01, 0.01])
x = 0
for i in range(height):
    if i%2 == 0:
        x = 0
    else:
        x = 1
    for j in range(blknum):
        if x == 0:
            blkT = [10+2*j, 0.5+i, 0]
            blkT = scale(blkT, s)
            addRigidBody(scene, '../models/cube_5.obj', 2, coScale=coblkScale,
                     scale=blkScale, translation=blkT,
                     dynamic=1, density = 400, rest = 0.5, friction = 0.1)
        else:
            blkT = [11, 0.5+i, -1+2*j]
            blkT = scale(blkT, s)
            addRigidBody(scene, '../models/cube_5.obj', 2, coScale=coblkScaleT,
                     scale=blkScaleT, translation=blkT,
                     dynamic=1, density = 400, rest = 0.5, friction = 0.1)
"""
# seesaw
seesawScale = [10, 0.5, 3]
seesawScale = scale(seesawScale, s)
seesawT = [-4, 4.0, 0]
seesawT = scale(seesawT,s)
ssBody = addRigidBody(scene, '../models/cube.obj', 2, coScale=seesawScale, 
             scale=seesawScale, translation=seesawT, 
             dynamic=1, rest= 0.5)
ssJointScale = [0.1, 1, 3]
ssJointScale = scale(ssJointScale, s)
ssJointT = add_vector(seesawT, [0,-1,0])
ssJointT = scale(ssJointT, s)
ssJointBody = addRigidBody(scene, '../models/cube_5.obj', 2, coScale=ssJointScale, 
             scale=ssJointScale, translation=ssJointT, 
             dynamic=0, rest= 0.5)
seesawJointT = add_vector(ssJointT, [0,0.75,0])
addHingeJoint(scene, ssBody, ssJointBody, seesawJointT, [0, 0, 1])
kinblkScale = [1, 0.5, 3]
kinblkScale = scale(kinblkScale, s)
kinblkT = add_vector(seesawT,[-5,-0.5,0])
kinblkT = scale(kinblkT, s)
addRigidBody(scene, '../models/cube_5.obj', 2, coScale=kinblkScale, 
             scale=kinblkScale, translation=kinblkT, 
             dynamic=0, rest= 0.5)
# block
blockScale = [2, 2, 2]
blockScale = scale(blockScale, s)
blockT = [0, 14, 0]
blockT = scale(blockT, s)
addRigidBody(scene, '../models/cube_5.obj', 2, coScale=blockScale, 
             scale=blockScale, translation=blockT, 
             dynamic=1, density = 5000, rest= 0.5)

#bomb
armadilloScale = [0.8,0.8,0.8]
bunnyScale = [4,4,4]
armadilloScale = scale(armadilloScale, s)
bunnyScale = scale(bunnyScale, s)
density = 500
friction = 0.1
restitution = 0.6
rotaxis = [0,1,0]
rotangle = 180.0
t = [-7.5, 5.5, 0]
addRigidBody(scene, '../models/armadillo.obj', 5, coFile='', coScale=armadilloScale, 
                         translation=t, scale=armadilloScale, axis = rotaxis, angle = rotangle,
                         dynamic=1, density=density, friction=friction, rest=restitution)
"""
writeScene(scene, 'MyScene.json')
