import torch
import math
import numpy as np
from torch.autograd import Function
import time

device = torch.device('cuda')

class RBRSLayer(Function):

    @staticmethod
    def forward(ctx, sim, task, P, D, movemesh, pos,lastPos,requires_grad=True):
        sim.stepwithbackward(False)
        #sim.step()
        
        newPos=sim.pos

        p_pt = sim.getdPTarget()
        p_dt = sim.getdDTarget()
        p_d = sim.getdD()
        p_l = sim.getdL()
        p_ll = sim.getdLL()
        
        scale=1e1
        partial_pt = torch.from_numpy(p_pt).to(device)/scale
        partial_dt = torch.from_numpy(p_dt).to(device)/scale
        partial_d = torch.from_numpy(p_d).to(device)/scale
        partial_l = torch.from_numpy(p_l).to(device)/scale
        partial_ll = torch.from_numpy(p_ll).to(device)/scale
        
        ctx.save_for_backward(partial_pt, partial_dt,partial_d,partial_l,partial_ll)
        return torch.from_numpy(newPos).to(device)

    @staticmethod
    def backward(ctx, grad_output):
        dp, dd, dD,dl, dll= ctx.saved_tensors
        grad_output=grad_output.to(device)
        return None, None, torch.matmul(grad_output.unsqueeze(0),dp).squeeze(1) , torch.matmul(grad_output.unsqueeze(0),dd).squeeze(1), torch.matmul(grad_output.unsqueeze(0),dD).squeeze(1),torch.matmul(grad_output.unsqueeze(0),dl).squeeze(1) , torch.matmul(grad_output.unsqueeze(0),dll).squeeze(1)
    
class PosLayer(Function):

    @staticmethod
    def forward(ctx, sim, pos, requires_grad=True):
        points=sim.getGlobalPoints(True)
        p = sim.getdPos().transpose()
        partial = torch.from_numpy(p).to(device)
        ctx.save_for_backward(partial)
        return torch.from_numpy(points).to(device)

    @staticmethod
    def backward(ctx, grad_output):
        d = ctx.saved_tensors[0]
        grad_output=grad_output.to(device)
        return None, torch.matmul(grad_output.unsqueeze(0),d).squeeze(1)