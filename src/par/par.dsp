# Microsoft Developer Studio Project File - Name="par" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=par - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "par.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "par.mak" CFG="par - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "par - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "par - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "par - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "../bio" /I ".." /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "par - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /I "../bio" /I ".." /D "_LIB" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "NO_HASHMAP" /D "_MOMENTS" /D "_LITTLE_ENDIAN" /D "_COMPARE" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "par - Win32 Release"
# Name "par - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\atom.cpp
# End Source File
# Begin Source File

SOURCE=.\barycentric.cpp
# End Source File
# Begin Source File

SOURCE=.\gridfun.cpp
# End Source File
# Begin Source File

SOURCE=.\mcconextractor.cpp
# End Source File
# Begin Source File

SOURCE=.\mtxlib.cpp
# End Source File
# Begin Source File

SOURCE=.\nnbucket.cpp
# End Source File
# Begin Source File

SOURCE=.\nnstruct.cpp
# End Source File
# Begin Source File

SOURCE=.\particlefun.cpp
# End Source File
# Begin Source File

SOURCE=.\point3d.cpp
# End Source File
# Begin Source File

SOURCE=.\protein.cpp
# End Source File
# Begin Source File

SOURCE=.\rand.cpp
# End Source File
# Begin Source File

SOURCE=.\reg3data.cpp
# End Source File
# Begin Source File

SOURCE=.\spherefun.cpp
# End Source File
# Begin Source File

SOURCE=.\surface3d.cpp
# End Source File
# Begin Source File

SOURCE=.\vecmath.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\atom.h
# End Source File
# Begin Source File

SOURCE=.\barycentric.h
# End Source File
# Begin Source File

SOURCE=.\conextractor.h
# End Source File
# Begin Source File

SOURCE=.\cubes.h
# End Source File
# Begin Source File

SOURCE=.\dynarray.h
# End Source File
# Begin Source File

SOURCE=.\elements.h
# End Source File
# Begin Source File

SOURCE=.\function.h
# End Source File
# Begin Source File

SOURCE=.\geom.h
# End Source File
# Begin Source File

SOURCE=.\gridfun.h
# End Source File
# Begin Source File

SOURCE=.\mcconextractor.h
# End Source File
# Begin Source File

SOURCE=.\mtxlib.h
# End Source File
# Begin Source File

SOURCE=.\nnbucket.h
# End Source File
# Begin Source File

SOURCE=.\nnkdtree.h
# End Source File
# Begin Source File

SOURCE=.\nnstruct.h
# End Source File
# Begin Source File

SOURCE=.\parsurface.h
# End Source File
# Begin Source File

SOURCE=.\particlefun.h
# End Source File
# Begin Source File

SOURCE=.\point3d.h
# End Source File
# Begin Source File

SOURCE=.\protein.h
# End Source File
# Begin Source File

SOURCE=.\rand.h
# End Source File
# Begin Source File

SOURCE=.\reg3data.h
# End Source File
# Begin Source File

SOURCE=.\spherefun.h
# End Source File
# Begin Source File

SOURCE=.\surface3d.h
# End Source File
# Begin Source File

SOURCE=.\triangulate.h
# End Source File
# Begin Source File

SOURCE=.\vecmath.h
# End Source File
# Begin Source File

SOURCE=.\vtkMarchingCubesCases.h
# End Source File
# End Group
# End Target
# End Project
