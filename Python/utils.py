import pyPBAD,math
import numpy as np

def D2R(angle):
    return angle*math.pi/180

def getBodyJointTrans(body,id):
    return np.array(body.getJoint(id).trans)

def setBodyJointTrans(body,id,trans):
    joint=body.getJoint(id)
    joint.trans=trans
    body.setJoint(id, joint)

def getBodyJointControl(body,id):
    return np.array(body.getJoint(id).control)

def setBodyJointControl(body,id,control):
    joint=body.getJoint(id)
    joint.control=control
    body.setJoint(id, joint)