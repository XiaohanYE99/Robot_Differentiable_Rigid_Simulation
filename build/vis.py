import vtk

# 读取 VTK 文件
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName("GJKConfig.vtk")
reader.Update()

# 创建 Mapper 和 Actor
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(reader.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)

# 创建渲染器和渲染窗口
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)

renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)

# 创建交互器并启动渲染窗口
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

# 自定义交互器样式，允许使用WSAD键调整图像大小
class CustomInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self):
        self.AddObserver("KeyPressEvent", self.keyPressEvent)

    def keyPressEvent(self, obj, event):
        key = obj.GetKeySym()
        if key == "w":
            renderer.GetActiveCamera().Zoom(1.1)
        elif key == "s":
            renderer.GetActiveCamera().Zoom(0.9)
        elif key == "a":
            renderer.GetActiveCamera().Dolly(1.1)
        elif key == "d":
            renderer.GetActiveCamera().Dolly(0.9)
        
        renderWindow.Render()
        return

# 设置自定义交互器样式
style = CustomInteractorStyle()
renderWindowInteractor.SetInteractorStyle(style)

# 启动渲染
renderWindow.Render()
renderWindowInteractor.Start()