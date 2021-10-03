# Microsoft Developer Studio Project File - Name="betti" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=betti - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "betti.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "betti.mak" CFG="betti - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "betti - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "betti - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "betti - Win32 Release"

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
# ADD CPP /nologo /MD /W3 /GX /O2 /I "../par" /I "../nr" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "NO_HASHMAP" /D "_MOMENTS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "betti - Win32 Debug"

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
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /I "../par" /I "../nr" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "NO_HASHMAP" /D "_MOMENTS" /YX /FD /GZ /c
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

# Name "betti - Win32 Release"
# Name "betti - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\actree.cpp
# End Source File
# Begin Source File

SOURCE=.\arc.cpp
# End Source File
# Begin Source File

SOURCE=.\contourtree.cpp
# End Source File
# Begin Source File

SOURCE=.\cut.cpp
# End Source File
# Begin Source File

SOURCE=.\disjointset.cpp
# End Source File
# Begin Source File

SOURCE=.\dual.cpp
# End Source File
# Begin Source File

SOURCE=.\interval_tree.cpp
# End Source File
# Begin Source File

SOURCE=.\moment.cpp
# End Source File
# Begin Source File

SOURCE=.\multicon.cpp
# End Source File
# Begin Source File

SOURCE=.\node.cpp
# End Source File
# Begin Source File

SOURCE=.\volume.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\actree.h
# End Source File
# Begin Source File

SOURCE=.\arc.h
# End Source File
# Begin Source File

SOURCE=.\basic.h
# End Source File
# Begin Source File

SOURCE=.\betti.h
# End Source File
# Begin Source File

SOURCE=.\bin.h
# End Source File
# Begin Source File

SOURCE=.\contourtree.h
# End Source File
# Begin Source File

SOURCE=.\CriticalPoint.h
# End Source File
# Begin Source File

SOURCE=.\cut.h
# End Source File
# Begin Source File

SOURCE=.\disjointset.h
# End Source File
# Begin Source File

SOURCE=.\dual.h
# End Source File
# Begin Source File

SOURCE=.\interval_tree.h
# End Source File
# Begin Source File

SOURCE=.\linkedlist.h
# End Source File
# Begin Source File

SOURCE=.\moment.h
# End Source File
# Begin Source File

SOURCE=.\momentum.h
# End Source File
# Begin Source File

SOURCE=.\multicon.h
# End Source File
# Begin Source File

SOURCE=.\node.h
# End Source File
# Begin Source File

SOURCE=.\pqueue.h
# End Source File
# Begin Source File

SOURCE=.\volume.h
# End Source File
# End Group
# End Target
# End Project
