set cmakecmd="C:\Program Files\CMake\bin\cmake"
set msbuildpath="C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\amd64\"
if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\amd64\MSBuild.exe" (
	set msbuildcmd="C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\amd64\MSBuild.exe"
) else if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\MSBuild\Current\Bin\amd64\MSBuild.exe" (
	set msbuildcmd="C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\MSBuild\Current\Bin\amd64\MSBuild.exe"
)
set currpath=%cd%

::Step 1: VCPKG
if not exist C:\vcpkg (
	cd /d C:\
	git clone https://github.com/Microsoft/vcpkg.git
	cd /d C:\vcpkg
	bootstrap-vcpkg.bat
)
cd /d C:\vcpkg
SET "PATH=C:\vcpkg\installed\x64-windows\bin;C:\vcpkg\installed\x64-windows\debug\bin;%PATH%"
vcpkg install tinyxml2:x64-windows
vcpkg install assimp:x64-windows
vcpkg install eigen3:x64-windows
vcpkg install boost:x64-windows
vcpkg install cgal:x64-windows

::Step 2: TinyVisualizer
if not exist C:\TinyVisualizer (
	cd /d C:\
	git clone https://github.com/gaoxifeng/TinyVisualizer.git
)
cd /d C:\TinyVisualizer
git pull
git submodule update --init --recursive
if not exist C:\TinyVisualizer-build mkdir C:\TinyVisualizer-build
cd /d C:\TinyVisualizer-build
%cmakecmd% ..\TinyVisualizer
%msbuildcmd% ALL_BUILD.vcxproj /property:Configuration=Release
%msbuildcmd% INSTALL.vcxproj /property:Configuration=Release

::Step 3: libdifferentiable
if not exist C:\libdifferentiable (
	cd /d C:\
	git clone https://runningblade@bitbucket.org/runningblade/libdifferentiable.git
)
cd /d C:\libdifferentiable
git checkout PBAD
git pull
git submodule update --init --recursive
if not exist C:\libdifferentiable-build mkdir C:\libdifferentiable-build
cd /d C:\libdifferentiable-build
%cmakecmd% ..\libdifferentiable
%msbuildcmd% ALL_BUILD.vcxproj /property:Configuration=Release

::go back
cd /d %currpath%