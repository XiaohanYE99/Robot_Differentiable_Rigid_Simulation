import pyPBAD
import numpy as np
import torch
from RBRS import *
import time
from torch.utils.tensorboard import SummaryWriter

torch.manual_seed(42) 

writer = SummaryWriter('./log_spi')
pi=math.pi
device = torch.device('cuda')
load=False
dt=0.01
sim=pyPBAD.ConvHullPBDSimulator(dt)
#sim.setOutput(True)

body=pyPBAD.ArticulatedLoader.createInitJoints(8,1000)
floor=pyPBAD.BBoxExact(np.array([-200,-200,-5]),np.array([200,200,-3]))
nrD=body.nrDOF()#DOF of body
nrJ=body.nrJ()#Joint number of body

it_num=1
horizon=1000#horizon per iteration

sim.setArticulatedBody(body)
sim.addShape(floor)
sim.gravity=np.array([0,0,-9.81])
Mesh=torch.autograd.Variable(torch.from_numpy(sim.getConvexPoints()).to(device),requires_grad=True)

#set joint class, do not perform collision detection on joints belonging to the same class
c=np.zeros(nrJ)
for i in range(0,nrJ):
    c[i]=0
sim.setclass(c)

loss_fn=torch.nn.MSELoss()
Run=RBRSLayer.apply  
Getpoints=PosLayer.apply

Target=torch.autograd.Variable(torch.zeros(Mesh.size(),dtype=torch.double).to(device),requires_grad=True)

def init_Mesh(path):
    lastMesh=torch.autograd.Variable(torch.from_numpy(np.load(path).reshape(-1)).to(device),requires_grad=True)
    print(lastMesh.size())
    DeltaMesh=lastMesh-Mesh
    sim.updateConvexPoints(DeltaMesh.detach().cpu().numpy())
    #print('moveMesh: ',DeltaMesh.max(),DeltaMesh.min())

init_Mesh('../data/sphere.npy')
for it in range(it_num):
    traj=[]
    sim.reset()
    pos=sim.pos
    lastPos=sim.pos
    pos=torch.from_numpy(pos).to(device)
    lastPos=torch.from_numpy(lastPos).to(device)
    sim.pos=pos.detach().cpu().numpy()
    #sim.setVel((pos.detach().cpu().numpy()-lastPos.detach().cpu().numpy())/dt)
    sim.vel=np.array((15,15,0))
    loss=0
    traj.append(pos.detach().cpu().numpy())
    for i in range(horizon):
        print(i)
        P=torch.zeros(nrD).to(device)
        D=torch.zeros(nrD).to(device)
        newPos=Run(sim,None,P,D,Mesh,pos,lastPos)
        points=Getpoints(sim,newPos)
        lastPos=pos
        pos=newPos
        if i%1==0:
            traj.append(newPos.detach().cpu().numpy())
        loss+=loss_fn(Target,points)
    
    loss=loss.to(device)
    loss.backward()

    print('it= ',it, ' loss= ',loss)
    #np.save('control/control_%d.npy'%it,control.detach().cpu().numpy())
    #writer.add_scalar('loss', loss, it)
    

