import pyPBAD
from utils import *
import pyTinyVisualizer as vis

class Visualizer:

    def __init__(self, sim, frames=None, headless=False):
        self.headless = headless
        self.drawer = vis.Drawer(['--headless','1' if self.headless else '0'])
        self.export = vis.CameraExportPlugin(vis.GLFW_KEY_2,vis.GLFW_KEY_3,"camera.dat")
        self.capturer = vis.CaptureGIFPlugin(vis.GLFW_KEY_1,"record.gif",self.drawer.FPS(),True)
        self.drawer.addPlugin(self.export)
        self.drawer.addPlugin(self.capturer)

        self.sim = sim
        self.cBody=[.3,.3,.3]
        self.cBB=[.3,.2,.1]
        self.bodyShapeC=vis.CompositeShape()
        self.envShapeC=vis.CompositeShape()
        self.bodyShape = self.visualizeArticulated(self.sim.getBody())
        self.envShape = self.visualizeEnvironment(self.sim.getShapes())
        self.bodyWireShape = self.visualizeArticulated(self.sim.getBody(), True)
        self.envWireShape = self.visualizeEnvironment(self.sim.getShapes(), True)
        self.bodyWireShape.setLineWidth(2)
        self.envWireShape.setLineWidth(2)
        self.drag = vis.MeshShape()

        self.drawer.addShape(self.bodyShapeC)
        self.drawer.addShape(self.envShapeC)
        self.drawer.addCamera3D(90,[0.,0.,1.])
        self.drawer.getCamera3D().setManipulator(vis.FirstPersonCameraManipulator(self.drawer.getCamera3D()))
        self.drawer.addLightSystem(1024,10,False)
        self.drawer.getLight().lightSz(10)
        self.sim.step()
        bb=self.sim.getContact().getBB()
        range=max(np.max(np.abs(bb.minCorner)),np.max(np.abs(bb.maxCorner)))
        for x in [-1,1]:
            for y in [-1,1]:
                for z in [-1,1]:
                    self.drawer.getLight().addLight(np.array([x,y,z])*range,[.5,.5,.5],[.5,.5,.5],[.5,.5,.5])

        #parameter
        self.nrJ=len(self.sim.getJoints())
        self.wire=False
        self.lastWire=True
        self.doSim=False
        self.nSub=1
        self.frames=frames
        self.frame_id=0
        self.maxStiffness=1e6
        self.output=False
        self.crossTerm=True
        self.CRBA=False

    def frame(self):
        if self.wire!=self.lastWire:
            self.lastWire=self.wire
            if self.wire:
                self.bodyShapeC.addShape(self.bodyWireShape)
                self.envShapeC.addShape(self.envWireShape)
                self.bodyShapeC.removeChild(self.bodyShape)
                self.envShapeC.removeChild(self.envShape)
            else:
                self.bodyShapeC.addShape(self.bodyShape)
                self.envShapeC.addShape(self.envShape)
                self.bodyShapeC.removeChild(self.bodyWireShape)
                self.envShapeC.removeChild(self.envWireShape)
        pos=self.sim.pos if self.frames is None else self.frames[self.frame_id]
        if self.wire:
            self.updateArticulatedBody(self.bodyWireShape,pos)
        else: self.updateArticulatedBody(self.bodyShape,pos)
        if not self.doSim:
            return
        if self.frames is None:
            for i in range(self.nSub):
                self.sim.step()
        else:
            self.frame_id+=1
            if self.frame_id>=len(self.frames):
                self.frame_id=0
        #drag
        if len(self.sim.getJoints())>self.nrJ:
            joint = self.sim.getJoints()[-1]
            info = pyPBAD.PBDArticulatedGradientInfo(self.sim.getBody(),pos)
            self.drag.clear()
            self.drag.addVertex(joint.posA(info))
            self.drag.addVertex(joint.posB(info))
            self.drag.addIndexSingle(0)
            self.drag.addIndexSingle(1)
            self.drag.setMode(vis.GL_LINES)
            self.drag.setColorDiffuse(vis.GL_LINES,.5,1.,.5)
            self.drag.setLineWidth(2)

    def key(self, wnd, key, scan, action, mods, captured):
        if captured:
            return
        elif key==vis.GLFW_KEY_R and action==vis.GLFW_PRESS:
            self.doSim = not self.doSim
        elif key==vis.GLFW_KEY_E and action==vis.GLFW_PRESS:
            self.wire = not self.wire
        elif key==vis.GLFW_KEY_F and action==vis.GLFW_PRESS:
            pov = vis.Povray("scene")
            self.drawer.drawPovray(pov)

    def mouse(self, wnd, button, action, mods, captured):
        if captured or self.wire:
            return
        elif button==vis.GLFW_MOUSE_BUTTON_2 and action==vis.GLFW_PRESS:
            joints = self.sim.getJoints()
            while len(joints)>self.nrJ:
                joints = joints[:,-1]
            self.sim.setJoints(joints)
            self.IAlpha=100
            ray=self.drawer.getCameraRay()
            k=0
            for i in range(self.sim.getBody().nrJ()):
                if self.sim.getBody().getJoint(i).mesh is not None:
                    jointShape = self.bodyShape.getChild(k)
                    hasI,self.IAlpha = jointShape.rayIntersect(ray,self.IAlpha)
                    if hasI:
                        info = pyPBAD.PBDArticulatedGradientInfo(self.sim.getBody(),self.sim.pos)
                        joint = pyPBAD.SoftJoint()
                        joint.jidA=i
                        #linearLimit
                        linearLimit=np.zeros((3,3))
                        linearLimit[0,:]=np.inf
                        linearLimit[1,:]=1
                        linearLimit[2,:]=self.sim.getJointPhysicsParameter(i).kc
                        joint.linearLimit=linearLimit
                        #transB
                        transB=np.zeros((3,4))
                        transB[:,3]=ray[:3]+ray[3:]*self.IAlpha
                        joint.transB=transB
                        #transA
                        transA=np.zeros((3,4))
                        transA[:,3]=info.T[:,4*i:4*i+3].T@(transB[:,3]-info.T[:,4*i+3])
                        joint.transA=transA
                        #drag
                        self.drag.clear()
                        self.drag.addVertex(joint.posA(info))
                        self.drag.addVertex(joint.posB(info))
                        self.drag.addIndexSingle(0)
                        self.drag.addIndexSingle(1)
                        self.drag.setMode(vis.GL_LINES)
                        self.drag.setColorDiffuse(vis.GL_LINES,.5,1.,.5)
                        self.drag.setLineWidth(2)
                        #add
                        self.drawer.addShape(self.drag)
                        joints.append(joint)
                        self.sim.setJoints(joints)
                    k+=1
        elif button==vis.GLFW_MOUSE_BUTTON_2 and action==vis.GLFW_RELEASE:
            self.drawer.removeShape(self.drag)
            joints = self.sim.getJoints()
            while len(joints)>self.nrJ:
                joints = joints[:-1]
            self.sim.setJoints(joints)

    def motion(self, wnd, x, y, captured):
        if captured:
            return
        elif len(self.sim.getJoints())>self.nrJ:
            ray=self.drawer.getCameraRay()
            info = pyPBAD.PBDArticulatedGradientInfo(self.sim.getBody(),self.sim.pos)
            joints = self.sim.getJoints()
            transB=np.zeros((3,4))
            transB[:,3]=ray[:3]+ray[3:]*self.IAlpha
            joints[-1].transB=transB
            # drag
            self.drag.clear()
            self.drag.addVertex(joints[-1].posA(info))
            self.drag.addVertex(joints[-1].posB(info))
            self.drag.addIndexSingle(0)
            self.drag.addIndexSingle(1)
            self.drag.setMode(vis.GL_LINES)
            self.drag.setColorDiffuse(vis.GL_LINES, .5, 1., .5)
            self.drag.setLineWidth(2)
            #add
            self.sim.setJoints(joints)

    def setup(self):
        self.drawer.getCamera3D().getManipulator().imGuiCallback()
        vis.imgui.Begin("Simulator",0)
        _,self.wire=vis.imgui.Checkbox("wireframe",self.wire)
        _,self.doSim=vis.imgui.Checkbox("simulating",self.doSim)
        vis.imgui.SameLine(0,-1)

        if vis.imgui.Button("reset position",vis.imgui.ImVec2(0,0)):
            self.sim.pos=[0]*self.sim.getBody().nrDOF()
            self.sim.vel=[0]*self.sim.getBody().nrDOF()
            self.sim.time=0
        vis.imgui.SameLine(0,-1)

        recording=self.capturer.recording()
        modified,recording=vis.imgui.Checkbox("record",recording)
        if modified:
            if recording:
                self.capturer.startRecording()
            else: self.capturer.stopRecording()
        if vis.imgui.Button("save camera",vis.imgui.ImVec2(0,0)):
            self.export.saveCamera()
        if vis.imgui.Button("load camera",vis.imgui.ImVec2(0,0)):
            self.export.loadCamera()
        vis.imgui.Separator()

        #friction coefficient
        fri=self.sim.getJointPhysicsParameter(0).friction
        modified,fri=vis.imgui.SliderFloat("friction coefficient",fri,0.1,1.0,"%.3f",0)
        if modified:
            for i in range(self.sim.getBody().nrJ()):
                p=self.sim.getJointPhysicsParameter(i)
                p.friction=fri
                self.sim.setJointPhysicsParameter(i,p)

        #contact stiffness
        stiff=self.sim.getHeuristcGuessStiffness()
        modified,stiff=vis.imgui.SliderFloat("contact stiffness",stiff,0,self.maxStiffness,"%.3f",0)
        if modified:
            self.sim.setHeuristcGuessStiffness(stiff)

        #gravity
        g=self.sim.gravity
        modified,g=vis.imgui.DragFloat3("gravity",g,0.01,-10,10,"%.3f",0)
        if modified:
            self.sim.gravity=g
        if vis.imgui.Button("reset gravity",vis.imgui.ImVec2(0,0)):
            self.sim.gravity=[0,0,-9.81]

        #output
        modified,self.output=vis.imgui.Checkbox("output",self.output)
        if modified:
            self.sim.setOutput(self.output)

        #CRBA
        if isinstance(self.sim,pyPBAD.PBDSimulator):
            vis.imgui.SameLine(0,-1)
            modified,self.crossTerm=vis.imgui.Checkbox("crossTerm",self.crossTerm)
            if modified:
                self.sim.setCrossTerm(self.crossTerm,self.CRBA)
            vis.imgui.SameLine(0,-1)
            modified,self.CRBA=vis.imgui.Checkbox("CRBA",self.CRBA)
            if modified:
                self.sim.setCrossTerm(self.crossTerm,self.CRBA)
        vis.imgui.End()

    def mainLoop(self):
        def frame(node):
            self.frame()
        def key(wnd,key,scan,action,mods,captured):
            self.key(wnd,key,scan,action,mods,captured)
        def mouse(wnd,button,action,mods,captured):
            self.mouse(wnd,button,action,mods,captured)
        def motion(wnd,x,y,captured):
            self.motion(wnd,x,y,captured)
        def setup():
            self.setup()
        self.drawer.addPlugin(vis.ImGuiPlugin(setup))
        self.drawer.setFrameFunc(frame)
        self.drawer.setKeyFunc(key)
        self.drawer.setMouseFunc(mouse)
        self.drawer.setMotionFunc(motion)
        self.drawer.mainLoop()

    def updateArticulatedBody(self,node:vis.CompositeShape,pos):
        j=0
        pos = pyPBAD.PBDArticulatedGradientInfo(self.sim.getBody(),pos).T
        for i in range(self.sim.getBody().nrJ()):
           if self.sim.getBody().getJoint(i).mesh is not None:
                s = node.getChild(j)
                M = pos[:,4*i:4*i+4]
                s.setLocalRotate(M[:,:3])
                s.setLocalTranslate(M[:,3])
                j+=1

    def visualizeArticulated(self, b:pyPBAD.ArticulatedBody, wire:bool=False, convex:bool=False):
        checker = vis.drawChecker(5,[1.,1.,1.],[.7,.7,.7], 9, 11, vis.GL_RGB)
        shape = vis.CompositeShape()
        for i in range(b.nrJ()):
            joint = b.getJoint(i)
            if joint.mesh is not None:
                stJ = vis.Bullet3DShape()
                if convex and not isinstance(joint.mesh, pyPBAD.ConvexHullExact):
                    vss,iss = joint.mesh.getMesh()
                    mesh = pyPBAD.ConvexHullExact(vss)
                    s = self.visualizeShapeExact(mesh, wire)
                else: s = self.visualizeShapeExact(joint.mesh, wire)
                stJ.addShape(s)
                if not wire:
                    s.setTextureDiffuse(checker)
                s.setColorDiffuse(vis.GL_TRIANGLES,self.cBody[0],self.cBody[1],self.cBody[2])
                s.setColorDiffuse(vis.GL_LINES,self.cBody[0],self.cBody[1],self.cBody[2])
                shape.addShape(stJ)
        return shape

    def visualizeEnvironment(self, ss, wire:bool=False):
        shape=vis.CompositeShape()
        for o in ss:
            shape.addShape(self.visualizeShapeExact(o,wire))
        return shape

    def visualizeShapeExact(self, s:pyPBAD.ShapeExact, wire:bool):
        if isinstance(s,pyPBAD.MeshExact):
            return self.visualizeMeshExact(s,wire)
        elif isinstance(s,pyPBAD.SphericalBBoxExact):
            return self.visualizeSphericalBBoxExact(s,wire)
        elif isinstance(s,pyPBAD.BBoxExact):
            return self.visualizeBBoxExact(s,wire)
        elif isinstance(s,pyPBAD.CompositeShapeExact):
            return self.visualizeCompositeExact(s,wire)
        else:
            print(f"visualization of {s} is not supported!")
            assert False

    def visualizeBBoxExact(self, m:pyPBAD.BBoxExact, wire:bool):
        checker=vis.drawGrid(25,0.001,0,[.8,.9,1.],[.1,.1,.1],9,11,vis.GL_RGB)
        shape=vis.Bullet3DShape()
        box=vis.makeBox(1,not wire,(m.maxCorner-m.minCorner)/2)
        box.setCastShadow(False)
        box.setTextureDiffuse(checker)
        box.setColorDiffuse(vis.GL_POINTS,.8,.9,1.)
        box.setColorDiffuse(vis.GL_LINES,.8,.9,1.)
        box.setColorDiffuse(vis.GL_TRIANGLES,.8,.9,1.)
        shape.setLocalTranslate((m.maxCorner+m.minCorner)/2)
        shape.addShape(box)
        return shape

    def visualizeSphericalBBoxExact(self, m:pyPBAD.SphericalBBoxExact, wire:bool):
        shape=vis.Bullet3DShape()
        sbox=vis.makeSphericalBox(8,not wire,m.radius(),(m.maxCorner-m.minCorner)/2)
        shape.setLocalTranslate((m.maxCorner+m.minCorner)/2)
        shape.addShape(sbox)
        return shape

    def visualizeMeshExact(self, m:pyPBAD.MeshExact, wire:bool):
        tri = vis.MeshShape()
        tcss = m.tcss()
        vss = m.vss()
        iss = m.iss()
        for i in range(len(vss)):
            if len(tcss)>0:
                tri.addVertex(vss[i], tcss[i])
            else: tri.addVertex(vss[i])
        if wire:
            for i in range(len(iss)):
                tri.addIndex([iss[i][0],iss[i][1]])
                tri.addIndex([iss[i][1],iss[i][2]])
                tri.addIndex([iss[i][2],iss[i][0]])
            tri.setMode(vis.GL_LINES)
        else:
            for i in range(len(iss)):
                tri.addIndex(iss[i])
            tri.setMode(vis.GL_TRIANGLES)
            tri.computeNormals()
        return tri

    def visualizeCompositeExact(self, m:pyPBAD.CompositeShapeExact, wire:bool):
        localT=np.eye(4)
        shape=vis.Bullet3DShape()
        geoms=m.getGeoms()
        trans=m.getTrans()
        materials=m.getMaterials()
        for i in range(len(geoms)):
            shapeI=vis.Bullet3DShape()
            if len(materials)>i:
                wire=materials[i].useWireframe
            shapeI.addShape(self.visualizeShapeExact(geoms[i],wire))
            if len(materials)>i:
                #ambient
                shapeI.setColorAmbient(vis.GL_POINTS,materials[i].ambient[0],materials[i].ambient[1],materials[i].ambient[2])
                shapeI.setColorAmbient(vis.GL_LINES,materials[i].ambient[0],materials[i].ambient[1],materials[i].ambient[2])
                shapeI.setColorAmbient(vis.GL_TRIANGLES,materials[i].ambient[0],materials[i].ambient[1],materials[i].ambient[2])
                #diffuse
                shapeI.setColorAmbient(vis.GL_POINTS,materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2])
                shapeI.setColorAmbient(vis.GL_LINES,materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2])
                shapeI.setColorAmbient(vis.GL_TRIANGLES,materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2])
                #specular
                shapeI.setColorAmbient(vis.GL_POINTS,materials[i].specular[0],materials[i].specular[1],materials[i].specular[2])
                shapeI.setColorAmbient(vis.GL_LINES,materials[i].specular[0],materials[i].specular[1],materials[i].specular[2])
                shapeI.setColorAmbient(vis.GL_TRIANGLES,materials[i].specular[0],materials[i].specular[1],materials[i].specular[2])
                #shininess
                shapeI.setShininess(vis.GL_POINTS,materials[i].shininess)
                shapeI.setShininess(vis.GL_LINES,materials[i].shininess)
                shapeI.setShininess(vis.GL_TRIANGLES,materials[i].shininess)
                #texture
                if materials[i].texFile != "":
                    tex=vis.Texture.load(materials[i].texFile)
                    print(f"reading texture from {materials[i].texFile} w={tex.width()} h={tex.height()}")
                    shapeI.setTextureDiffuse(tex)
            #assign
            localT[:3,:4]=trans[i]
            shapeI.setLocalTransform(localT)
            shape.addShape(shapeI)
        return shape