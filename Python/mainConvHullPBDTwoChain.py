import pyPBAD

from utils import *
from visualizer import *

if __name__=='__main__':
    floor=pyPBAD.BBoxExact([-150,-150,-8],[30,30,-6])

    bodies=[]
    #first chain
    bodies.append(
        pyPBAD.ArticulatedLoader.createChain(int(pyPBAD.JOINT_TYPE.TRANS_3D) | int(pyPBAD.JOINT_TYPE.ROT_3D_EXP), 20, 0.5, 0.24,
                                          D2R(30), 0, 0, 0, 3, 0, 0, 0))
    #second chain
    bodies.append(
        pyPBAD.ArticulatedLoader.createChain(int(pyPBAD.JOINT_TYPE.TRANS_3D) | int(pyPBAD.JOINT_TYPE.ROT_3D_EXP), 20, 0.5, 0.24,
                                          D2R(30), 0, 0, 0, 3, 0, 0, 0))
    #third chain
    bodies.append(
        pyPBAD.ArticulatedLoader.createChain(int(pyPBAD.JOINT_TYPE.TRANS_3D) | int(pyPBAD.JOINT_TYPE.ROT_3D_EXP), 20, 0.5, 0.24,
                                          D2R(30), 0, 0, 0, 3, 0, 0, 0))
    #forth chain
    bodies.append(
        pyPBAD.ArticulatedLoader.createChain(int(pyPBAD.JOINT_TYPE.ROT_3D_EXP), 20, 0.5, 0.24,
                                          D2R(30), 0, 0, 0, 3, 0, 0, 0))
    #combine
    body=pyPBAD.ArticulatedBody()
    pyPBAD.ArticulatedUtils(body).combine(bodies)
    #offset
    trans=getBodyJointTrans(body,22)
    trans[2,3]=1
    setBodyJointTrans(body,22,trans)
    #offset
    trans=getBodyJointTrans(body,43)
    trans[2,3]=2
    setBodyJointTrans(body,43,trans)
    #offset
    trans=getBodyJointTrans(body,64)
    trans[2,3]=3
    setBodyJointTrans(body,64,trans)
    
    #simulator
    sim=pyPBAD.ConvHullPBDSimulator(0.033)
    sim.setArticulatedBody(body)
    sim.addShape(floor)
    sim.gravity=[0,0,-9.81]
    sim.setOutput(True)

    #setup kinematic 1
    def PFunc1(time,d):
        return np.array([math.cos(time*2),math.sin(time*2),0])*4.5
    def DFunc1(time,d):
        return np.array([-math.sin(time*2),math.cos(time*2),0])*9
    for k in range(1,3):
        #control
        ctrl=getBodyJointControl(body,k)
        for d in range(ctrl.shape[0]):
            ctrl[d]=1
        setBodyJointControl(body,k,ctrl)
        print(f"Joints[{k}] has {body.getJoint(k).nrDOF()} DOF")
        #param
        param=sim.getJointPhysicsParameter(k)
        param.kp=1000
        param.kd=10
        param.setTarPFunc(PFunc1)
        param.setTarDFunc(DFunc1)
        sim.setJointPhysicsParameter(k,param)

    #setup kinematic 2
    def PFunc2(time,d):
        return np.array([math.cos(time*2),math.sin(time*2),0])*3
    def DFunc2(time,d):
        return np.array([-math.sin(time*2),math.cos(time*2),0])*6
    for k in range(22,24):
        #control
        ctrl=getBodyJointControl(body,k)
        for d in range(ctrl.shape[0]):
            ctrl[d]=1
        setBodyJointControl(body,k,ctrl)
        print(f"Joints[{k}] has {body.getJoint(k).nrDOF()} DOF")
        #param
        param=sim.getJointPhysicsParameter(k)
        param.kp=1000
        param.kd=10
        param.setTarPFunc(PFunc2)
        param.setTarDFunc(DFunc2)
        sim.setJointPhysicsParameter(k,param)

    #setup kinematic 3
    def PFunc3(time,d):
        return np.array([math.cos(time*2),math.sin(time*2),0])*1.5
    def DFunc3(time,d):
        return np.array([-math.sin(time*2),math.cos(time*2),0])*3
    for k in range(43,45):
        #control
        ctrl=getBodyJointControl(body,k)
        for d in range(ctrl.shape[0]):
            ctrl[d]=1
        setBodyJointControl(body,k,ctrl)
        print(f"Joints[{k}] has {body.getJoint(k).nrDOF()} DOF")
        #param
        param=sim.getJointPhysicsParameter(k)
        param.kp=1000
        param.kd=10
        param.setTarPFunc(PFunc3)
        param.setTarDFunc(DFunc3)
        sim.setJointPhysicsParameter(k,param)

    #run simulator
    runFrames=False
    if runFrames:
        frames=[]
        for i in range(10):
            print(f"Computing frame {i}")
            sim.step()
            frames.append(sim.pos)
    else: frames=None

    #visualize
    vis = Visualizer(sim, frames)
    vis.mainLoop()