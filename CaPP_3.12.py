##    CaPP A.B                       ##
##                                   ##
##    Copyright 2020,                ##
##    Andreas Larsen                 ##
##    andreas.larsen@bioch.ox.ac.uk  ##

## This file is part of CaPP.                                           ##
##                                                                      ##
## CaPP is free a software: you can redistribute it and/or modify       ##
## it under the terms of the GNU General Public License as published by ##
## the Free Software Foundation, either version 3 of the License, or    ##
## (at your option) any later version.                                  ##
##                                                                      ##
## CaPP is distributed in the hope that it will be useful,              ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of       ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        ##
## GNU General Public License for more details.                         ##
##       <http://www.gnu.org/licenses/>                                 ##
##                                                                      ##

## If you use this software in                                          ##
## your work, please write:                                             ##
## CaPP, source code freely available at github.com/Niels-Bohr-Institute-XNS-StructBiophys/CaPP ##

### import Libraries
import re
import os
import sys
import math
import linecache
import re
import platform
import ntpath
import numpy as np
try:
    import wx
except:
    print("")
    print("******************************************************************")
    print(" CaPP failed to import wxPython - is it correctly installed...?   ")
    print("******************************************************************")
    print("")
    sys.exit(1)
import wx.lib.scrolledpanel


### ignore "Future Warnings"
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
print('x' in np.arange(5))   #returns False, without Warning

### check OS and use correct binary. Check location of binary 
if platform.system() == "Darwin":
    capp_version = "capp_mac"
elif platform.system() == "Windows":
    capp_version = "capp/capp.exe"
    if programpath.count(" ") > 0:
        print("")
        print("****************************************************************************************")
        print("Found path is: ", programpath)
        print("Path must not contain spaces!")
        print("****************************************************************************************")
        sys.exit(1)
else:
    capp_version = "capp"
programpath = os.path.dirname(os.path.realpath(__file__))
capp_executable_location = programpath + "/" + capp_version
print(capp_executable_location)
if os.path.exists(capp_executable_location):
    print("CaPP initiated correctly")
elif capp_version == "capp":
    print("")
    print("****************************************************************************************")
    print(" Cannot find the executable: %s. It should be in the same folder as CaPP.py " % capp_version)
    print(" If you are not running MacOS, then:")
    print("  1) compile the source code Mainfunction.c with output called capp")
    print("     e.g. gcc Mainfunction.c -o capp -lm")
    print("  2) place the executable, capp, in the same folder as CaPP.py")
    print("  3) re-run CaPP with: python CaPP.py")
    print("****************************************************************************************")
    print("")
    sys.exit(1)
elif capp_version == "capp/capp.exe":
    print("")
    print("****************************************************************************************")
    print("Could not find Windows executables. Compile (gcc Mainfunction.c -o capp -lm) and place capp.exe here: ", programpath + "/" + capp_version)
    print("****************************************************************************************")
    sys.exit(1)
else:
    print("")
    print("****************************************************************************************")
    print(" Cannot find the executable: %s. It should be in the same folder as CaPP.py " % capp_version)
    print("****************************************************************************************")
    print("")
    sys.exit(1)

## Import plotting libraries
import matplotlib
matplotlib.interactive(True)
matplotlib.use('WXAgg')
import pylab

## import fitting libraries
from scipy.optimize import curve_fit

## Define main class and text
class MainCls(wx.Frame):
    
    ##################################### Function that creates the GUI ###########################################
    def __init__(self, parent, id):

        ### Overall frame for widgets
        width = wx.SystemSettings.GetMetric(wx.SYS_SCREEN_X)
        height = wx.SystemSettings.GetMetric(wx.SYS_SCREEN_Y)
        position = (width/3.0, height/30.0)
        wx.Frame.__init__(self, parent, id, title = 'CaPP - Calculate p(r) for Proteins (PDB format)', pos = position)
        BoxSizer = wx.BoxSizer(wx.VERTICAL)
        self.Panel = wx.lib.scrolledpanel.ScrolledPanel(self, -1)
        self.Panel.SetupScrolling(False, True)
        self.Bind(wx.EVT_CLOSE, self.CloseWindowFnc) # close pylab explicitely for to avoid crash
        BoxSizer.AddSpacer(3) # vertical spacing
        
        ### Widgets for PDB-file import
        PDBTextSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDBTextSizer.AddSpacer(10)
        PDBText = wx.StaticText(self.Panel, -1, 'Location of PDB file:', size = (330, -1))
        PDBTextSizer.Add(PDBText)
        BoxSizer.Add(PDBTextSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        PDBPathSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDBPathSizer.AddSpacer(10)
        self.PDBPathTxt = wx.StaticText(self.Panel, -1, '')
        self.PDBPathStr = 'N/A'
        PDBPathSizer.Add(self.PDBPathTxt)
        BoxSizer.Add(PDBPathSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        PDBBtnSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDBBtnSizer.AddSpacer(10)
        BrowsePDBBtn = wx.Button(self.Panel, label = 'Browse')
        self.Bind(wx.EVT_BUTTON, self.BrowsePDBFnc, BrowsePDBBtn)
        PDBBtnSizer.Add(BrowsePDBBtn, 1, wx.EXPAND)
        PDBBtnSizer.AddSpacer(10)
        BoxSizer.Add(PDBBtnSizer, 0, wx.EXPAND|wx.HORIZONTAL)

        BoxSizer.AddSpacer(5)
        LinePDB = wx.StaticLine(self.Panel, -1)
        BoxSizer.Add(LinePDB, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(5)

        ### Widgets for Data-file import
        DataTextSizer = wx.BoxSizer(wx.HORIZONTAL)
        DataTextSizer.AddSpacer(10)
        DataText = wx.StaticText(self.Panel, -1, 'Location of Data file (optional, for q-range and fit):', size = (330, -1))
        DataTextSizer.Add(DataText)
        BoxSizer.Add(DataTextSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(5)
        DataPathSizer = wx.BoxSizer(wx.HORIZONTAL)
        DataPathSizer.AddSpacer(10)
        self.DataPathTxt = wx.StaticText(self.Panel, -1, '')
        self.DataPathStr = 'Non'
        DataPathSizer.Add(self.DataPathTxt)
        BoxSizer.Add(DataPathSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        
        DataBtnSizer = wx.BoxSizer(wx.HORIZONTAL)
        DataBtnSizer.AddSpacer(10)
        self.BrowseDataBtn = wx.Button(self.Panel, label = 'Browse')
        self.Bind(wx.EVT_BUTTON, self.BrowseDataFnc, self.BrowseDataBtn)
        DataBtnSizer.Add(self.BrowseDataBtn, 1, wx.EXPAND)
        DataBtnSizer.AddSpacer(10)
        BoxSizer.Add(DataBtnSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        self.BrowseDataBtn.Disable()
        
        nm_button = wx.BoxSizer(wx.HORIZONTAL)
        self.nm_button = wx.CheckBox(self.Panel, -1, 'q in data is in [1/nm] (default: [1/aa])', size = (350, -1))
        nm_button.Add(self.nm_button)
        BoxSizer.Add(nm_button, 0, wx.ALIGN_CENTER, wx.EXPAND|wx.HORIZONTAL)
        self.Bind(wx.EVT_CHECKBOX, self.nmFnc, self.nm_button)
        self.nm_button.Disable()
        
        BoxSizer.AddSpacer(5)
        LineData = wx.StaticLine(self.Panel, -1)
        BoxSizer.Add(LineData, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(5)

        ### Widgets used to choose between SAXS and SANS
        SAXS_button = wx.BoxSizer(wx.HORIZONTAL) # prepare for SAXS button
        SAXS_button.AddSpacer(10)
        self.SAXS_button = wx.RadioButton(self.Panel, -1, 'SAXS',style=wx.RB_GROUP) # make SAXS button
        SAXS_button.Add(self.SAXS_button)
        BoxSizer.Add(SAXS_button, 0, wx.EXPAND|wx.HORIZONTAL) # position SAXS button
        self.Bind(wx.EVT_RADIOBUTTON, self.EnableSAXSFnc, self.SAXS_button)
        
        SAXS_solvent_box = wx.BoxSizer(wx.HORIZONTAL) # prepare for SAXS solvent box
        SAXS_solvent_box.AddSpacer(10)
        
        self.left_SAXS = wx.StaticText(self.Panel, -1, "Sucrose content =")
        SAXS_solvent_box.Add(self.left_SAXS)
        self.left_SAXS.Disable()
        
        self.SAXS_solvent_box = wx.TextCtrl(self.Panel, -1, '', size = (40, -1)) # make SAXS solvent box
        self.SAXS_solvent_box.SetValue('0')
        SAXS_solvent_box.Add(self.SAXS_solvent_box) # add SAXS botton to its own line
        BoxSizer.Add(SAXS_solvent_box) # position SAXS solvent box
        self.SAXS_solvent_box.Disable()
        
        self.right_SAXS = wx.StaticText(self.Panel, -1, "g/100 ml (%)")
        SAXS_solvent_box.Add(self.right_SAXS)
        self.right_SAXS.Disable()
        
        SANS_button = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS button
        SANS_button.AddSpacer(10) # vertical spacing
        self.SANS_button = wx.RadioButton(self.Panel, -1, 'SANS') # make SANS button
        SANS_button.Add(self.SANS_button) # add SANS_button to its own line
        BoxSizer.Add(SANS_button, 0, wx.EXPAND|wx.HORIZONTAL) # position SANS button
        self.Bind(wx.EVT_RADIOBUTTON, self.EnableSANSFnc, self.SANS_button)

        SANS_solvent_box = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS solvent box
        SANS_solvent_box.AddSpacer(10)

        self.left_SANS = wx.StaticText(self.Panel, -1, "D2O content, solvent     =")
        SANS_solvent_box.Add(self.left_SANS)
        self.left_SANS.Disable()

        self.SANS_solvent_box = wx.TextCtrl(self.Panel, -1, '', size = (40, -1)) # make SANS solvent box
        self.SANS_solvent_box.SetValue('100')
        SANS_solvent_box.Add(self.SANS_solvent_box) # add SANS solven box to its own line
        BoxSizer.Add(SANS_solvent_box) # position SANS solvent box
        self.SANS_solvent_box.Disable()

        self.right_SANS = wx.StaticText(self.Panel, -1, "%")
        SANS_solvent_box.Add(self.right_SANS)
        self.right_SANS.Disable()

        SANS_perdeut_box_A = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS perdeuteration box A
        SANS_perdeut_box_A.AddSpacer(10)

        self.left_SANS_d = wx.StaticText(self.Panel, -1, "Perdeut., chain A,B,...,G =")
        SANS_perdeut_box_A.Add(self.left_SANS_d)
        self.left_SANS_d.Disable()

        self.SANS_perdeut_box_A = wx.TextCtrl(self.Panel, -1, '', size = (40, -1)) # make SANS perdeuteration box
        self.SANS_perdeut_box_A.SetValue('0')
        SANS_perdeut_box_A.Add(self.SANS_perdeut_box_A) # add SANS_perdeut_box_A to its own line
        BoxSizer.Add(SANS_perdeut_box_A) # position SANS perdeuteration box
        self.SANS_perdeut_box_A.Disable()

        SANS_perdeut_box_B = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS perdeuteration box B
        self.SANS_perdeut_box_B = wx.TextCtrl(self.Panel, -1, '', size = (40, -1)) # make SANS perdeuteration box
        self.SANS_perdeut_box_B.SetValue('0')
        SANS_perdeut_box_A.Add(self.SANS_perdeut_box_B) # add SANS_perdeut_box_B to same line as SANS_perdeut_box_A
        BoxSizer.Add(SANS_perdeut_box_B) # position SANS perdeuteration box B
        self.SANS_perdeut_box_B.Disable()

        SANS_perdeut_box_C = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS perdeuteration box C
        self.SANS_perdeut_box_C = wx.TextCtrl(self.Panel, -1, '', size = (40, -1)) # make SANS perdeuteration box
        self.SANS_perdeut_box_C.SetValue('0')
        SANS_perdeut_box_A.Add(self.SANS_perdeut_box_C) # add SANS_perdeut_box_C to same line as SANS_perdeut_box_A
        BoxSizer.Add(SANS_perdeut_box_C) # position SANS perdeuteration box C
        self.SANS_perdeut_box_C.Disable()

        SANS_perdeut_box_D = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS perdeuteration box D
        self.SANS_perdeut_box_D = wx.TextCtrl(self.Panel, -1, '', size = (40, -1)) # make SANS perdeuteration box
        self.SANS_perdeut_box_D.SetValue('0')
        SANS_perdeut_box_A.Add(self.SANS_perdeut_box_D) # add SANS_perdeut_box_D to same line as SANS_perdeut_box_A
        BoxSizer.Add(SANS_perdeut_box_D) # position SANS perdeuteration box D
        self.SANS_perdeut_box_D.Disable()

        SANS_perdeut_box_E = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS perdeuteration box E
        self.SANS_perdeut_box_E = wx.TextCtrl(self.Panel, -1, '', size = (40, -1)) # make SANS perdeuteration box
        self.SANS_perdeut_box_E.SetValue('0')
        SANS_perdeut_box_A.Add(self.SANS_perdeut_box_E) # add SANS_perdeut_box_E to same line as SANS_perdeut_box_A
        BoxSizer.Add(SANS_perdeut_box_E) # position SANS perdeuteration box E
        self.SANS_perdeut_box_E.Disable()

        SANS_perdeut_box_F = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS perdeuteration box E
        self.SANS_perdeut_box_F = wx.TextCtrl(self.Panel, -1, '', size = (40, -1)) # make SANS perdeuteration box
        self.SANS_perdeut_box_F.SetValue('0')
        SANS_perdeut_box_A.Add(self.SANS_perdeut_box_F) # add SANS_perdeut_box_F to same line as SANS_perdeut_box_A
        BoxSizer.Add(SANS_perdeut_box_F) # position SANS perdeuteration box F
        self.SANS_perdeut_box_F.Disable()

        SANS_perdeut_box_G = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS perdeuteration box E
        self.SANS_perdeut_box_G = wx.TextCtrl(self.Panel, -1, '', size = (40, -1)) # make SANS perdeuteration box
        self.SANS_perdeut_box_G.SetValue('0')
        SANS_perdeut_box_A.Add(self.SANS_perdeut_box_G) # add SANS_perdeut_box_G to same line as SANS_perdeut_box_A
        BoxSizer.Add(SANS_perdeut_box_G) # position SANS perdeuteration box F
        self.SANS_perdeut_box_G.Disable()

        self.right_SANS_d = wx.StaticText(self.Panel, -1, "%")
        SANS_perdeut_box_A.Add(self.right_SANS_d)
        self.right_SANS_d.Disable()

        Perdeut_text = wx.BoxSizer(wx.HORIZONTAL)
        Perdeut_text.AddSpacer(10)
        self.Perdeut_text = wx.StaticText(self.Panel, -1, 'Enter a number between 0 and 100 %, OBS:', size = (450, -1))
        Perdeut_text.Add(self.Perdeut_text)
        BoxSizer.Add(Perdeut_text, 0, wx.EXPAND|wx.HORIZONTAL)
        self.Perdeut_text.Disable()

        Perdeut_text2 = wx.BoxSizer(wx.HORIZONTAL)
        Perdeut_text2.AddSpacer(10)
        self.Perdeut_text2 = wx.StaticText(self.Panel, -1, 'For single-chain proteins (no chain labels), value in "chain A" box is global.', size = (450, -1))
        Perdeut_text2.Add(self.Perdeut_text2)
        BoxSizer.Add(Perdeut_text2, 0, wx.EXPAND|wx.HORIZONTAL)
        self.Perdeut_text2.Disable()
        BoxSizer.AddSpacer(10)

        RES_button = wx.BoxSizer(wx.HORIZONTAL)
        self.RES_button = wx.CheckBox(self.Panel, -1, 'Include resolution effects (from 4th column in data)', size = (350, -1))
        RES_button.Add(self.RES_button)
        BoxSizer.Add(RES_button, 0, wx.ALIGN_CENTER, wx.EXPAND|wx.HORIZONTAL)
        self.Bind(wx.EVT_CHECKBOX, self.resFnc, self.RES_button)
        self.RES_button.SetValue(0)
        self.RES_button.Disable()

        Line_after_contrast_widget = wx.StaticLine(self.Panel, -1) # make line after widget
        BoxSizer.Add(Line_after_contrast_widget, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(5)

        ### Widgets used to add Water Layer
        Water_layer_button = wx.BoxSizer(wx.HORIZONTAL)
        self.Water_layer_button = wx.CheckBox(self.Panel, -1, 'Add Water Layer (WL):', size = (350, -1))
        Water_layer_button.Add(self.Water_layer_button)
        BoxSizer.Add(Water_layer_button, 0, wx.ALIGN_CENTER, wx.EXPAND|wx.HORIZONTAL)
        self.Bind(wx.EVT_CHECKBOX, self.EnableWLButtonsFnc, self.Water_layer_button)

        Water_layer_contrast_box = wx.BoxSizer(wx.HORIZONTAL)
        Water_layer_contrast_box.AddSpacer(29)

        self.left_WL = wx.StaticText(self.Panel, -1, "WL contrast =")
        Water_layer_contrast_box.Add(self.left_WL)
        self.left_WL.Disable()
        self.Water_layer_contrast_box = wx.TextCtrl(self.Panel, -1, '', size = (40, -1))
        self.Water_layer_contrast_box.SetValue('6')
        Water_layer_contrast_box.Add(self.Water_layer_contrast_box)
        BoxSizer.Add(Water_layer_contrast_box)
        self.Water_layer_contrast_box.Disable()
        self.right_WL = wx.StaticText(self.Panel, -1, "% of solvent scattering")
        Water_layer_contrast_box.Add(self.right_WL)
        self.right_WL.Disable()

        BoxSizer.AddSpacer(5) # vertical spacing
        
        ### Widgets to exclude water layer from TMD
        WLText = wx.BoxSizer(wx.HORIZONTAL)
        WLText.AddSpacer(10)
        self.WLText = wx.StaticText(self.Panel, -1, 'If the protein has a Transmembrane Domain (TMD):', size = (330, -1))
        WLText.Add(self.WLText)
        BoxSizer.Add(WLText, 0, wx.EXPAND|wx.HORIZONTAL)
        self.WLText.Disable()

        No_Exclude_water_layer_button = wx.BoxSizer(wx.HORIZONTAL)
        No_Exclude_water_layer_button.AddSpacer(10)
        self.No_Exclude_water_layer_button = wx.RadioButton(self.Panel, -1, 'Do not exclude the WL from the TMD', style=wx.RB_GROUP)
        No_Exclude_water_layer_button.Add(self.No_Exclude_water_layer_button)
        BoxSizer.Add(No_Exclude_water_layer_button, 0, wx.EXPAND|wx.HORIZONTAL)
        self.No_Exclude_water_layer_button.Disable()
        self.Bind(wx.EVT_RADIOBUTTON, self.DisableBilayerThicknessFnc, self.No_Exclude_water_layer_button)

        Exclude_water_layer_OPM_button = wx.BoxSizer(wx.HORIZONTAL)
        Exclude_water_layer_OPM_button.AddSpacer(10)
        self.Exclude_water_layer_OPM_button = wx.RadioButton(self.Panel, -1, 'Exclude the WL from the TMD, using OPM(*)')
        Exclude_water_layer_OPM_button.Add(self.Exclude_water_layer_OPM_button)
        BoxSizer.Add(Exclude_water_layer_OPM_button, 0, wx.EXPAND|wx.HORIZONTAL)
        self.Exclude_water_layer_OPM_button.Disable()
        self.Bind(wx.EVT_RADIOBUTTON, self.DisableBilayerThicknessFnc, self.Exclude_water_layer_OPM_button)

        Exclude_water_layer_M_button = wx.BoxSizer(wx.HORIZONTAL)
        Exclude_water_layer_M_button.AddSpacer(10)
        self.Exclude_water_layer_M_button = wx.RadioButton(self.Panel, -1, 'Exclude WL from TMD, manually(**):')
        Exclude_water_layer_M_button.Add(self.Exclude_water_layer_M_button)
        BoxSizer.Add(Exclude_water_layer_M_button, 0, wx.EXPAND|wx.HORIZONTAL)
        self.Exclude_water_layer_M_button.Disable()
        self.Bind(wx.EVT_RADIOBUTTON, self.EnableBilayerThicknessFnc, self.Exclude_water_layer_M_button)

        Bilayer_thickness_box = wx.BoxSizer(wx.HORIZONTAL)
        Bilayer_thickness_box.AddSpacer(29)

        self.left_thick = wx.StaticText(self.Panel, -1, "Bilayer Thickness =")
        Bilayer_thickness_box.Add(self.left_thick)
        self.left_thick.Disable()

        self.Bilayer_thickness_box = wx.TextCtrl(self.Panel, -1, '', size = (40, -1))
        self.Bilayer_thickness_box.SetValue('30.0')
        Bilayer_thickness_box.Add(self.Bilayer_thickness_box)
        BoxSizer.Add(Bilayer_thickness_box, 0)
        self.Bilayer_thickness_box.Disable()

        self.right_thick = wx.StaticText(self.Panel, -1, "Angstrom")
        Bilayer_thickness_box.Add(self.right_thick)
        self.right_thick.Disable()

        OPMText = wx.BoxSizer(wx.HORIZONTAL)
        OPMText.AddSpacer(10)
        self.OPMText = wx.StaticText(self.Panel, -1, '(*) Orientations of Proteins in Membranes Database', size = (330, -1))
        OPMText.Add(self.OPMText)
        BoxSizer.Add(OPMText, 0, wx.EXPAND|wx.HORIZONTAL)
        self.OPMText.Disable()
        
        MText1 = wx.BoxSizer(wx.HORIZONTAL)
        MText1.AddSpacer(10)
        self.MText1 = wx.StaticText(self.Panel, -1, '(**) Align protein with TMD perp. to xy-plane, and TMD-center in z=0.', size = (430, -1))
        MText1.Add(self.MText1)
        BoxSizer.Add(MText1, 0, wx.EXPAND|wx.HORIZONTAL)
        self.MText1.Disable()

        BoxSizer.AddSpacer(5)
        Line_after_exclude_WL_widget = wx.StaticLine(self.Panel, -1) # make line after widget
        BoxSizer.Add(Line_after_exclude_WL_widget, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(5)

        ### Widgets for calculation buttons
        CalculateButtonSpace = wx.BoxSizer(wx.HORIZONTAL)
        CalculateButtonSpace.AddSpacer(10)
        self.CalculateButton = wx.Button(self.Panel, label = 'Calculate p(r), Rg, and P(q)')
        self.Bind(wx.EVT_BUTTON, self.cappFnc, self.CalculateButton)
        self.CalculateButton.Disable()
        CalculateButtonSpace.Add(self.CalculateButton, 1, wx.EXPAND)
        CalculateButtonSpace.AddSpacer(10)
        BoxSizer.Add(CalculateButtonSpace, 0, wx.EXPAND|wx.HORIZONTAL)
        
        BoxSizer.AddSpacer(5)
        Line_after_widget = wx.StaticLine(self.Panel, -1) # make line after widget
        BoxSizer.Add(Line_after_widget, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(5)

        ### Widgets for fitting
        FitButtonSpace = wx.BoxSizer(wx.HORIZONTAL)
        FitButtonSpace.AddSpacer(10)
        self.FitButton = wx.Button(self.Panel, label = 'Fit P(q) to data')
        self.Bind(wx.EVT_BUTTON, self.FitFnc, self.FitButton)
        FitButtonSpace.Add(self.FitButton, 1, wx.EXPAND)
        FitButtonSpace.AddSpacer(10)
        BoxSizer.Add(FitButtonSpace, 0, wx.EXPAND|wx.HORIZONTAL)
        self.FitButton.Disable()
        BoxSizer.AddSpacer(1)

        BoxSizer.AddSpacer(4)
        Skip_box = wx.BoxSizer(wx.HORIZONTAL)
        Skip_box.AddSpacer(5)
        self.left_skip = wx.StaticText(self.Panel, -1, "Skip the first ")
        Skip_box.Add(self.left_skip)
        self.left_skip.Disable()
        self.Skip_box = wx.TextCtrl(self.Panel, -1, '', size = (40, -1))
        self.Skip_box.SetValue('0')
        Skip_box.Add(self.Skip_box)
        BoxSizer.Add(Skip_box)
        self.Skip_box.Disable()
        self.right_skip = wx.StaticText(self.Panel, -1, " points in q (default: 0)")
        Skip_box.Add(self.right_skip)
        self.right_skip.Disable()
        
        fitWL_button = wx.BoxSizer(wx.HORIZONTAL)
        self.fitWL_button = wx.CheckBox(self.Panel, -1, 'Fit WL contrast? (default: fixed to chosen value)', size = (350, -1))
        fitWL_button.Add(self.fitWL_button)
        BoxSizer.Add(fitWL_button, 0, wx.ALIGN_CENTER, wx.EXPAND|wx.HORIZONTAL)
        self.fitWL_button.Disable()
        
        fitPDB2_button = wx.BoxSizer(wx.HORIZONTAL)
        self.fitPDB2_button = wx.CheckBox(self.Panel, -1, 'Fit with 2 PDBs? Location of 2nd PDB file:', size = (350, -1))
        fitPDB2_button.Add(self.fitPDB2_button)
        self.Bind(wx.EVT_CHECKBOX, self.EnableBrowsePDB2Fnc, self.fitPDB2_button)
        BoxSizer.Add(fitPDB2_button, 0, wx.ALIGN_CENTER, wx.EXPAND|wx.HORIZONTAL)
        self.fitPDB2_button.Disable()
        
        PDB2PathSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDB2PathSizer.AddSpacer(10)
        self.PDB2PathTxt = wx.StaticText(self.Panel, -1, '')
        self.PDB2PathStr = 'N/A'
        PDB2PathSizer.Add(self.PDB2PathTxt)
        BoxSizer.Add(PDB2PathSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        
        PDB2BtnSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDB2BtnSizer.AddSpacer(10)
        self.BrowsePDB2Btn = wx.Button(self.Panel, label = 'Browse')
        self.Bind(wx.EVT_BUTTON, self.BrowsePDB2Fnc, self.BrowsePDB2Btn)
        PDB2BtnSizer.Add(self.BrowsePDB2Btn, 1, wx.EXPAND)
        PDB2BtnSizer.AddSpacer(10)
        BoxSizer.Add(PDB2BtnSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        self.BrowsePDB2Btn.Disable()
        
        BoxSizer.AddSpacer(5)
        Line_after_widget = wx.StaticLine(self.Panel, -1) # make line after widget
        BoxSizer.Add(Line_after_widget, 0, wx.EXPAND|wx.HORIZONTAL)
        #BoxSizer.AddSpacer(5)
        
        ### Widgets to change resolution (binsize)
        Change_resolution_button = wx.BoxSizer(wx.HORIZONTAL)
        self.Change_resolution_button = wx.CheckBox(self.Panel, -1, 'Change resolution (binsize) of the p(r):', size = (270, -1))
        Change_resolution_button.Add(self.Change_resolution_button)
        BoxSizer.Add(Change_resolution_button, 0, wx.ALIGN_CENTER, wx.EXPAND|wx.HORIZONTAL)
        self.Bind(wx.EVT_CHECKBOX, self.EnableResolutionBoxFnc, self.Change_resolution_button)
        Resolution_box = wx.BoxSizer(wx.HORIZONTAL)        
        self.Resolution_box = wx.TextCtrl(self.Panel, -1, '', size = (40, -1))
        self.Resolution_box.SetValue('1.0')
        Change_resolution_button.Add(self.Resolution_box)
        self.Resolution_box.Disable()
        self.right_res = wx.StaticText(self.Panel, -1, "A")
        Change_resolution_button.Add(self.right_res)
        self.right_res.Disable()


        ### Widgets to choose explicit H and D
        Explicit_button = wx.BoxSizer(wx.HORIZONTAL)
        self.Explicit_button = wx.CheckBox(self.Panel, -1, 'Account for H and D explicitely? (default: no)', size = (320, -1))
        self.Bind(wx.EVT_CHECKBOX, self.ExplicitFnc, self.Explicit_button)
        Explicit_button.Add(self.Explicit_button)
        BoxSizer.Add(Explicit_button, 0, wx.ALIGN_CENTER, wx.EXPAND|wx.HORIZONTAL)        
        BoxSizer.AddSpacer(10)

        ### Conclusion of init-function
        self.Panel.SetSizerAndFit(BoxSizer)
        self.SetSizerAndFit(BoxSizer)
        self.SetAutoLayout(True)
        self.Panel.Layout()
        self.Layout()

    ##################################### FUNCTIONS related to GUI ###########################################

    ### define function for closing GUI
    def onClose(self,event):
        self.Destroy()

    ### define function for printing q unit 
    def nmFnc(self,event):
        if self.nm_button.GetValue():
            print("q in data is in [1/nm]")
        else:
            print("q in data is in [1/aa]")

    ### define function for printing if resution effects are included
    def resFnc(self,event):
        if self.RES_button.GetValue():
            print("SANS resolution effects included (from 4th column in data)")
        else:
            print("SANS resolutioin effects not included (default)")

    ### define funciton for enable/disable "Browse PDB 2" botton
    def EnableBrowsePDB2Fnc(self,event):
        if self.fitPDB2_button.GetValue():
            self.BrowsePDB2Btn.Enable()
        else:
            self.BrowsePDB2Btn.Disable()

    ### define function for enable/disable GUI related to WL
    def EnableWLButtonsFnc(self, event):
        if self.Water_layer_button.GetValue():
            self.left_WL.Enable()
            self.Water_layer_contrast_box.Enable()
            self.right_WL.Enable()
            self.Exclude_water_layer_OPM_button.Enable()
            self.Exclude_water_layer_M_button.Enable()
            self.WLText.Enable()
            self.No_Exclude_water_layer_button.Enable()
            self.OPMText.Enable()
            self.MText1.Enable()
            self.fitWL_button.Enable()
        else:
            self.left_WL.Disable()
            self.Water_layer_contrast_box.Disable()
            self.right_WL.Disable()
            self.Exclude_water_layer_OPM_button.Disable()
            self.Exclude_water_layer_M_button.Disable()
            self.Bilayer_thickness_box.Disable()
            self.WLText.Disable()
            self.No_Exclude_water_layer_button.Disable()
            self.right_thick.Disable()
            self.left_thick.Disable()
            self.OPMText.Disable()
            self.MText1.Disable()
            self.Exclude_water_layer_M_button.SetValue(False)
            self.Exclude_water_layer_OPM_button.SetValue(False)
            self.No_Exclude_water_layer_button.SetValue(False)
            self.fitWL_button.Disable()
            self.fitWL_button.SetValue(0)

    ### define function for enable/disable GUI related to resolutions of p(r), i.e binsize
    def EnableResolutionBoxFnc(self, event):
        if self.Change_resolution_button.GetValue():
            self.Resolution_box.Enable()
            self.right_res.Enable()
            Message = 'If you have calculated with other resolution before, then you need to manually delete/rename/move those output files to get the right results. Keep plot windows open for easy comparison.'
            print(Message)
            wx.MessageBox(Message, "Changing resolution/binsize of p(r)", wx.OK | wx.ICON_INFORMATION)
        else:
            self.Resolution_box.Disable()
            self.right_res.Disable()

    ### define function for explicit H/D button
    def ExplicitFnc(self, event):
        if self.Explicit_button.GetValue():
            Message = 'If you have calculated with implicit H/D before, then you need to manually delete/rename/move those output files to get the right results. Keep plot windows open for easy comparison.'
            print(Message)
            wx.MessageBox(Message, "Changing to explicit H/D", wx.OK | wx.ICON_INFORMATION)
        else:
            Message = 'If you have calculated with explicit H/D before, then you need to manually delete/rename/move those output files to get the right results. Keep plot windows open for easy comparison.'
            print(Message)
            wx.MessageBox(Message, "Changing to implicit H/D", wx.OK | wx.ICON_INFORMATION)

    ### define function for enable GUI related to bilayer thickness
    def EnableBilayerThicknessFnc(self, event):
        self.Bilayer_thickness_box.Enable()
        self.right_thick.Enable()
        self.left_thick.Enable()

    ### define function for disable GUI related to bilayer thickness
    def DisableBilayerThicknessFnc(self, event):
        self.Bilayer_thickness_box.Disable()
        self.right_thick.Disable()
        self.left_thick.Disable()
    
    ### define function for enable/disable GUI related to SANS
    def EnableSANSFnc(self, event):
        self.left_SANS.Enable()
        self.left_SANS_d.Enable()
        self.SANS_solvent_box.Enable()
        self.SANS_perdeut_box_A.Enable()
        self.SANS_perdeut_box_B.Enable()
        self.SANS_perdeut_box_C.Enable()
        self.SANS_perdeut_box_D.Enable()
        self.SANS_perdeut_box_E.Enable()
        self.SANS_perdeut_box_F.Enable()
        self.SANS_perdeut_box_G.Enable()
        self.Perdeut_text.Enable()
        self.Perdeut_text2.Enable()
        self.right_SANS.Enable()
        self.right_SANS_d.Enable()
        if self.DataPathStr == 'Non':
            self.RES_button.Disable()
        else:
            self.RES_button.Enable()
        self.left_SAXS.Disable()
        self.SAXS_solvent_box.Disable()
        self.right_SAXS.Disable()

    ### define function for enable/disable GUI related to SAXS    
    def EnableSAXSFnc(self, event):
        self.left_SAXS.Enable()
        self.SAXS_solvent_box.Enable()
        self.right_SAXS.Enable()
        self.left_SANS.Disable()
        self.left_SANS_d.Disable()
        self.SANS_solvent_box.Disable()
        self.SANS_perdeut_box_A.Disable()
        self.SANS_perdeut_box_B.Disable()
        self.SANS_perdeut_box_C.Disable()
        self.SANS_perdeut_box_D.Disable()
        self.SANS_perdeut_box_E.Disable()
        self.SANS_perdeut_box_F.Disable()
        self.SANS_perdeut_box_G.Disable()
        self.Perdeut_text.Disable()
        self.Perdeut_text2.Disable()
        self.right_SANS.Disable()
        self.right_SANS_d.Disable()
        self.RES_button.SetValue(0)
        self.RES_button.Disable()
    
    ### Define function for filebrowsing for PDB
    def BrowsePDBFnc(self,event):
        
        FileDialogWindow = wx.FileDialog(None, 'Please select PDB file...', os.getcwd(), defaultFile = '')

        if FileDialogWindow.ShowModal() == wx.ID_OK:
            self.PDBPathStr = FileDialogWindow.GetPath()
            PDBPathDisplayStr = str(self.PDBPathStr)

            while len(PDBPathDisplayStr) > 49:
                PDBPathDisplayStr = PDBPathDisplayStr[1:]

            if len(self.PDBPathStr) > 49:
                PDBPathDisplayStr = '...' + PDBPathDisplayStr

            self.PDBPathTxt.SetLabel(PDBPathDisplayStr)

            DirectoryStr = os.path.dirname(self.PDBPathStr) # get current directory
            os.chdir(DirectoryStr) # cd to current directory

            self.CalculateButton.Enable()
            self.BrowseDataBtn.Enable()
            self.nm_button.Enable()
    
        FileDialogWindow.Destroy()

    ### Define function for filebrowsing for PDB2
    def BrowsePDB2Fnc(self,event):
    
        FileDialogWindow = wx.FileDialog(None, 'Please select second PDB file...', os.getcwd(), defaultFile = '')
        
        if FileDialogWindow.ShowModal() == wx.ID_OK:
            self.PDB2PathStr = FileDialogWindow.GetPath()
            PDB2PathDisplayStr = str(self.PDB2PathStr)
            
            while len(PDB2PathDisplayStr) > 49:
                PDB2PathDisplayStr = PDB2PathDisplayStr[1:]
            
            if len(self.PDB2PathStr) > 49:
                PDB2PathDisplayStr = '...' + PDB2PathDisplayStr
        
            self.PDB2PathTxt.SetLabel(PDB2PathDisplayStr)
            
            DirectoryStr = os.path.dirname(self.PDB2PathStr) # get current directory
            os.chdir(DirectoryStr) # cd to current directory
                    
        FileDialogWindow.Destroy()
    
    ### Define function for filebrowsing for data
    def BrowseDataFnc(self, event):
        FileDialogWindow = wx.FileDialog(None, 'Please select Data-file...', os.getcwd(), defaultFile = '')
        
        if FileDialogWindow.ShowModal() == wx.ID_OK:
            self.DataPathStr = FileDialogWindow.GetPath()
            DataPathDisplayStr = str(self.DataPathStr)
            
            while len(DataPathDisplayStr) > 49:
                DataPathDisplayStr = DataPathDisplayStr[1:]
            
            if len(self.DataPathStr) > 49:
                DataPathDisplayStr = '...' + DataPathDisplayStr
        
            self.DataPathTxt.SetLabel(DataPathDisplayStr)
    
            self.FitButton.Enable()
            self.left_skip.Enable()
            self.Skip_box.Enable()
            self.right_skip.Enable()
            self.fitPDB2_button.Enable()
            if self.SANS_button.GetValue():
                self.RES_button.Enable()
        
        FileDialogWindow.Destroy()
    
    ### Define exit function
    def CloseWindowFnc(self, event):
        try:
            pylab.close()
        except:
            pass
        sys.exit(0)

    ### Define funciton to extract file name from a file path
    def path_leaf(self,path):
        head, tail = ntpath.split(path)
        return tail or ntpath.basename(head) 

    ##################################### FUNCTIONS for calculations ###########################################

    ### define function for running CaPP (c-part) to calculate p(r)
    def cappFnc(self, event):

        # open pop-up window
        width = wx.SystemSettings.GetMetric(wx.SYS_SCREEN_X)
        height = wx.SystemSettings.GetMetric(wx.SYS_SCREEN_Y)
        position = (width/3.0, height/3.0)
        self.second_window = wx.Frame(None, title="Calculating p(r) - see Progression in the Terminal window...", size=(450,50), pos = position)
        self.second_window.Show()

        # import options from GUI
        XN_choice           = ""
        solvent             = ""
        Perdeut_choice_A    = ""
        Perdeut_choice_B    = ""
        Perdeut_choice_C    = ""
        Perdeut_choice_D    = ""
        Perdeut_choice_E    = ""
        Perdeut_choice_F    = ""
        Perdeut_choice_G    = ""
        perdeut_A           = ""
        perdeut_B           = ""
        perdeut_C           = ""
        perdeut_D           = ""
        perdeut_E           = ""
        perdeut_F           = ""
        perdeut_G           = ""
        X_choice            = ""
        PrcSucrose          = ""
        WL_choice           = ""
        WL_contrast         = ""
        Exclude_WL_choice   = ""
        Bilayer_thickness   = ""
        Resolution_choice   = ""
        Resolution          = ""
        Explicit_choice     = ""

        if self.Water_layer_button.GetValue():
            WL_choice               = "-c"
            WL_contrast_tmp         = float(self.Water_layer_contrast_box.GetValue())/100.0
            WL_contrast             = "%6.4f" % WL_contrast_tmp
        if self.Exclude_water_layer_OPM_button.GetValue():
            Exclude_WL_choice       = "-d"
        if self.Exclude_water_layer_M_button.GetValue():
            Exclude_WL_choice       = "-m"
            Bilayer_thickness       = self.Bilayer_thickness_box.GetValue()
        if self.Change_resolution_button.GetValue():
            Resolution_choice       = "-r"
            Resolution              = self.Resolution_box.GetValue()
        if self.Explicit_button.GetValue():
            Explicit_choice         = "-H"
        if self.SAXS_button.GetValue():
            X_choice                = "-x"
            PrcSucrose_tmp = float(self.SAXS_solvent_box.GetValue())
            PrcSucrose              = "%5.2f" % PrcSucrose_tmp
            if PrcSucrose_tmp > 200:
                Message = 'Maximum solubility of sucrose in water is 200%, i.e 200g/100ml'
                print(Message)
                wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)
        elif self.SANS_button.GetValue():
            XN_choice               = "-s"
            solvent_tmp = float(self.SANS_solvent_box.GetValue())/100.0
            solvent                 = "%4.3f" % solvent_tmp
            perdeut_tmp = float(self.SANS_perdeut_box_A.GetValue())/100.0
            Perdeut_choice_A          = "-A"
            perdeut_A                 = "%4.3f" % perdeut_tmp 
            perdeut_tmp = float(self.SANS_perdeut_box_B.GetValue())/100.0
            if perdeut_tmp > 0.0:
                Perdeut_choice_B      = "-B"
                perdeut_B             = "%4.3f" % perdeut_tmp
            perdeut_tmp = float(self.SANS_perdeut_box_C.GetValue())/100.0
            if perdeut_tmp > 0.0:
                Perdeut_choice_C       = "-C"
                perdeut_C              = "%4.3f" % perdeut_tmp
            perdeut_tmp = float(self.SANS_perdeut_box_D.GetValue())/100.0
            if perdeut_tmp > 0.0:
                Perdeut_choice_D      = "-D"
                perdeut_D             = "%4.3f" % perdeut_tmp
            perdeut_tmp = float(self.SANS_perdeut_box_E.GetValue())/100.0
            if perdeut_tmp > 0.0:
                Perdeut_choice_E      = "-E"
                perdeut_E             = "%4.3f" % perdeut_tmp
            perdeut_tmp = float(self.SANS_perdeut_box_F.GetValue())/100.0
            if perdeut_tmp > 0.0:
                Perdeut_choice_F      = "-F"
                perdeut_F             = "%4.3f" % perdeut_tmp
            perdeut_tmp = float(self.SANS_perdeut_box_G.GetValue())/100.0
            if perdeut_tmp > 0.0:
                Perdeut_choice_G      = "-G"
                perdeut_G             = "%4.3f" % perdeut_tmp
        else:
            X_choice                = "-x"
            PrcSucrose_tmp = float(self.SAXS_solvent_box.GetValue())
            PrcSucrose              = "%5.2f" % PrcSucrose_tmp
            Message = 'SAXS chosen as default'
            print(Message)
            wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)
    
        # assemble options into a command line and run the command line
        Output = "%s/%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % (programpath, capp_version, XN_choice, solvent, Perdeut_choice_A, perdeut_A, Perdeut_choice_B, perdeut_B, Perdeut_choice_C, perdeut_C, Perdeut_choice_D, perdeut_D, Perdeut_choice_E, perdeut_E, Perdeut_choice_F, perdeut_F, Perdeut_choice_G, perdeut_G, X_choice, PrcSucrose, WL_choice, WL_contrast, Exclude_WL_choice, Bilayer_thickness, Resolution_choice, Resolution, Explicit_choice, self.PDBPathStr)
        if self.fitWL_button.GetValue() == 0:
            print("Command generated by the GUI and sent to terminal:")
            print(Output)
            print("")
        os.system(Output)
        if self.fitPDB2_button.GetValue():
            Output = "%s/%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s" % (programpath, capp_version, XN_choice, solvent, Perdeut_choice_A, perdeut_A, Perdeut_choice_B, perdeut_B, Perdeut_choice_C, perdeut_C, Perdeut_choice_D, perdeut_D, Perdeut_choice_E, perdeut_E, Perdeut_choice_F, perdeut_F, Perdeut_choice_G, perdeut_G, X_choice, PrcSucrose, WL_choice, WL_contrast, Exclude_WL_choice, Bilayer_thickness, Resolution_choice, Resolution, Explicit_choice, self.PDB2PathStr)
            if self.fitWL_button.GetValue() == 0:
                print("Command generated by the GUI and sent to terminal (2nd PDB):")
                print(Output)
                print("")
            os.system(Output)

        ## close pop-up window
        self.second_window.Destroy()

        # generate filenames, plot p(r), calculate and plot P(q)
        filename,filename_nowater,short_filename,PDBID = self.GenerateFilenamesFnc(1)
        self.PlotprFnc(filename,PDBID) # plot p(r)   
        self.CalcPqFnc(filename,filename_nowater,short_filename,PDBID) # calculate and plot P(q)
        if self.fitPDB2_button.GetValue():
            filename,filename_nowater,short_filename,PDBID = self.GenerateFilenamesFnc(2)
            self.PlotprFnc(filename,PDBID) # plot p(r)   
            self.CalcPqFnc(filename,filename_nowater,short_filename,PDBID) # calculate and plot P(q)

    ### define function for generating filenames
    def GenerateFilenamesFnc(self,PDBnumber):
        if PDBnumber == 2:
            short_filename = os.path.splitext(self.PDB2PathStr)[0]
            PDBID = os.path.splitext(os.path.basename(self.PDB2PathStr))[0]
        else:
            short_filename = os.path.splitext(self.PDBPathStr)[0]
            PDBID = os.path.splitext(os.path.basename(self.PDBPathStr))[0]

        # check if user has taken pdb with water
        if short_filename[-2:] == '_w':
            if self.Water_layer_button.GetValue():
                Message = 'Your input PDB name ends with "_w", used by CaPP to denote PDB with water layer added. Please use the PDB without water if it is PDB generated by CaPP, or rename the PDB.'
                print(Message)
                wx.MessageBox(Message, "Error: PDB name with \"_w\" ending", wx.OK | wx.ICON_INFORMATION)
            else:  
                Message = 'Warning: your input PDB name ends with "_w", used by CaPP to denote PDB with water layer added. If it is a PDB generated by CaPP, then use the PDB without water, else consider renaming the PDB to avoid potential mistakes regarding added water layer. '
                print(Message)
                wx.MessageBox(Message, "Warning: PDB name with \"_w\" ending", wx.OK | wx.ICON_INFORMATION) 
        if short_filename[-7:] == '_w_only':
                Message = 'Error: You are trying to calculate p(r) and P(q) for water layer alone. Not possible this way. p(r) for water alone has been calculated already (output file \"PDBname_w_only_pr.dat\").'
                print(Message)
                wx.MessageBox(Message, "Error: PDB name with \"_w_only\" ending", wx.OK | wx.ICON_INFORMATION) 
        
        if self.Water_layer_button.GetValue():
            filename = short_filename + "_w" 
        else:
            filename = short_filename

        if self.SANS_button.GetValue():
            string_value_B = ""
            string_value_C = ""
            string_value_D = ""
            string_value_E = ""
            string_value_F = ""
            string_value_G = ""
            value = float(self.SANS_solvent_box.GetValue())
            string_value1 = "%1.0f" % value
            value = float(self.SANS_perdeut_box_A.GetValue())
            string_value_A = "%1.0f" % value 
            value = float(self.SANS_perdeut_box_B.GetValue())
            if value > 0.0:
                string_value_B = "_B%1.0f" % value
            value = float(self.SANS_perdeut_box_C.GetValue())
            if value > 0.0:
                string_value_C = "_C%1.0f" % value
            value = float(self.SANS_perdeut_box_D.GetValue())
            if value > 0.0:
                string_value_D = "_D%1.0f" % value
            value = float(self.SANS_perdeut_box_E.GetValue())
            if value > 0.0:
                string_value_E = "_E%1.0f" % value
            value = float(self.SANS_perdeut_box_F.GetValue())
            if value > 0.0:
                string_value_F = "_F%1.0f" % value
            value = float(self.SANS_perdeut_box_G.GetValue())
            if value > 0.0:
                string_value_G = "_G%1.0f" % value
            filename = filename + "_N" + string_value1 + "_P" + string_value_A + string_value_B  + string_value_C + string_value_D + string_value_E + string_value_F + string_value_G
            filename_nowater = short_filename + "_N" + string_value1 + "_P" + string_value_A + string_value_B + string_value_C + string_value_D + string_value_E + string_value_F + string_value_G             
        else:
            value = float(self.SAXS_solvent_box.GetValue())
            string_value = "%1.0f" % value
            filename = filename + "_X" + string_value
            filename_nowater = short_filename + "_X" + string_value

        return filename,filename_nowater,short_filename,PDBID

    ### define function for plotting p(r)
    def PlotprFnc(self,filename,PDBID):
      
        #import p(r)
        pr_filename = filename + "_pr.dat"
        r,pr = np.genfromtxt(pr_filename, skip_header=10,usecols=[0,1],unpack=True)
        
        # normalize p(r)
        pr = pr / max(pr)
        
        #plot p(r)
        if self.fitWL_button.GetValue() == 0:
            self.Figure1 = pylab.figure(1)
            Subplot = pylab.subplot(111)
            if self.Water_layer_button.GetValue():
                if self.Exclude_water_layer_M_button.GetValue() or self.Exclude_water_layer_OPM_button.GetValue():
                    structure_name = PDBID + " + WL " + self.Water_layer_contrast_box.GetValue()  + "% (Excl. TMD)"
                else:
                    structure_name = PDBID + " + WL " + self.Water_layer_contrast_box.GetValue()  + "%"
            else:
                structure_name = PDBID

            if self.SANS_button.GetValue():
                structure_name = structure_name + ", SANS " + self.SANS_solvent_box.GetValue() + "% D2O, [" + self.SANS_perdeut_box_A.GetValue() + "," + self.SANS_perdeut_box_B.GetValue() + "," + self.SANS_perdeut_box_C.GetValue() + "," + self.SANS_perdeut_box_D.GetValue() + "," + self.SANS_perdeut_box_E.GetValue() + "," + self.SANS_perdeut_box_F.GetValue() + "," + self.SANS_perdeut_box_G.GetValue() + "] % Deut."
            else:
                structure_name = structure_name + ", SAXS " + self.SAXS_solvent_box.GetValue() + "% Sucr."
            Subplot.plot(r, pr, label=structure_name)
            Subplot.legend(fontsize=9)
            Subplot.set_xlabel('r [$\AA$]',fontsize=14)
            Subplot.set_ylabel('p(r)',fontsize=14)
            pylab.suptitle('Pair Distance Distribution for the structure(s)',fontsize=14)
            pylab.show()
    
    ### define function to calculate, plot and export P(q)
    def CalcPqFnc(self,filename,filename_nowater,short_filename,PDBID):

        # import p(r) protein
        pr_filename = filename + "_pr.dat"
        r = np.genfromtxt(pr_filename, skip_header=10,usecols=[0],unpack=True)
        grp = r
        hrp = r
        jrp = r
        krp = r
        grw = r
        hrw = r
        jrw = r
        krw = r
        grc = r
        hrc = r
        jrc = r
        krc = r
        gr  = r
        hr  = r
        jr  = r
        kr  = r
        prw_filename = ""
        prc_filename = ""
        prt_filename = ""
        
        if self.Water_layer_button.GetValue():
            pr_filename = filename_nowater + "_pr.dat" #protein (no water)
            grp,hrp,jrp,krp = np.genfromtxt(pr_filename, skip_header=10,usecols=[2,3,4,5],unpack=True)
            prw_filename = short_filename + "_w_only_pr.dat" #water layer
            grw,hrw,jrw,krw = np.genfromtxt(prw_filename, skip_header=10,usecols=[2,3,4,5],unpack=True)
            prc_filename = short_filename + "_cross_pr.dat" #cross terms
            grc,hrc,jrc,krc = np.genfromtxt(prc_filename, skip_header=10,usecols=[2,3,4,5],unpack=True)
            prt_filename = filename + "_pr.dat" #total (proten + water)
            gr,hr,jr,kr = np.genfromtxt(prt_filename, skip_header=10,usecols=[2,3,4,5],unpack=True)
        else:
            gr,hr,jr,kr = np.genfromtxt(pr_filename, skip_header=10,usecols=[2,3,4,5],unpack=True)

        # import mean volume
        v_mean_line = linecache.getline(pr_filename,8)
        v_mean = float(v_mean_line.split(':')[1])

        # create or import q vector
        if self.DataPathStr == 'Non':
            PointsInQ = 200
            qmin = 0.001
            qmax = 1.0
            dq = (qmax - qmin) / (PointsInQ - 1)
            q = np.arange(qmin, qmax+dq, dq) # create q vector
        else:
            #count number of lines in data file
            NumberOfLines = sum(1 for line in open(self.DataPathStr)) + 1

            #count number of header lines in data file
            file = open(self.DataPathStr,'r')
            headerlines = 0
            datalines = 0
            footerlines = 0
            totallines = 0
            columns_in_datafile = 0
            HEADER_END = 0
            while totallines < NumberOfLines:
                line = file.readline()
                numbers_str = line.split()
                totallines = totallines + 1
                try:
                    numbers_float = [float(x) for x in numbers_str]
                    columns = len(numbers_float)
                    if columns < 3 or columns > 4:
                        if HEADER_END == 0:
                            headerlines = headerlines + 1
                        else:
                            if totallines < NumberOfLines:
                                footerlines = footerlines + 1
                    else:
                        datalines = datalines + 1 
                        columns_in_datafile = columns
                    HEADER_END = 1
                except:
                    if HEADER_END == 0:
                        headerlines = headerlines + 1
                    else:
                        footerlines = footerlines + 1
            file.close()
            if datalines == 0:
                Message = 'Cannot read datafile. Datafile should have 3 or 4 columns: q, I, and sigma_q (sigma_q is optional: SANS resolution effects).'
                print(Message)
                wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)    
            else:
                if self.RES_button.GetValue():
                    if columns_in_datafile == 3:
                        Message = 'Datafile should have 4 columns: q, I, sigma_I, sigma_q - 4th for SANS resolution effects.' 
                        print(Message)    
                        wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)    
            file.close()
            q = np.genfromtxt(self.DataPathStr,skip_header=headerlines,skip_footer=footerlines,usecols=[0], unpack=True) #import q vector
            if self.nm_button.GetValue():
                q = q*0.1 # convert from nm to angstrom
            dq = q
            if self.RES_button.GetValue():
                dq = np.genfromtxt(self.DataPathStr, skip_header=headerlines,skip_footer=footerlines,usecols=[3], unpack=True) #import dq vector
                sigma = np.array([-3.,-2.,-1.,0.,1.,2.,3.])
                q_MAT = np.matrix([q, q, q, q, q, q, q])
                for i in range(0,len(sigma)):
                    q_MAT[i] = sigma[i] * dq + q
        PointsInQ = len(q)

        ### define inner function to calculate Pq 
        def CalcPqSubFnc(q,v_mean,r, grp,hrp,jrp, grw,hrw,jrw,krw, grc,jrc,hrc,krc,gr,hr,jr,kr):

            # atomic form factor carbon (used for all atoms)
            b_c = 6.0 # scattering length of carbon [in units of electron scattering lengths]
            a = [2.31,1.02,1.5886,0.865] # parameters in carbon form factor amplitude
            b = [20.8439,10.2075,0.5687,51.6512] # parameters in carbon form factor amplitude
            c = 0.2156 # parameters in carbon form factor amplitude

            # calculate for each q:
            #      carbon atomic form factor amplitude = psi_c
            #      mean Gassian sphere for excluded volume = Gauss_sphere_mean
            #      form factor P(q) = Pq
            PointsInQ = len(q)
            PointsInr = len(r)
            Pq = np.zeros(PointsInQ)
            inv4pi = 1.0/(4.0*math.pi)
            inv4pi2 = inv4pi*inv4pi
            qr = 0.0
            for i in range(0,PointsInQ):
                q2 = q[i]*q[i]
                if self.SANS_button.GetValue():
                    psi_c = 1.0
                else:
                    psi_sum = c
                    for j in range(0,len(a)):
                        psi_sum = psi_sum + a[j] * math.exp(-b[j] * q2 * inv4pi2)
                    psi_c = psi_sum/b_c # normalize to psi_c(q = 0) = 1. f_c = psi_c * b_c, where b_c = 6.0
                psi_c2 = psi_c * psi_c
                Gauss_sphere_mean = math.exp(- q2 * v_mean**(2.0/3.0) * inv4pi)
                Gauss_sphere_mean2 = Gauss_sphere_mean * Gauss_sphere_mean
                Psum = 0.0
                for j in range(1,PointsInr):
                    qr = q[i] * r[j]
                    sinc = math.sin(qr)/qr
                    if self.Water_layer_button.GetValue():
                        sump = grp[j] * psi_c2 - (hrp[j] + jrp[j]) * psi_c * Gauss_sphere_mean + krp[j]  * Gauss_sphere_mean2
                        sumw = (grw[j] - hrw[j] - jrw[j] + krw[j]) * Gauss_sphere_mean2
                        sumc = (grc[j] - jrc[j]) * psi_c * Gauss_sphere_mean - (hrc[j] - krc[j]) * Gauss_sphere_mean2
                        Psum = Psum + (sump + sumw + sumc) * sinc
                    else:
                        Psum = Psum + (gr[j] * psi_c2 - (hr[j] + jr[j]) * psi_c * Gauss_sphere_mean + kr[j] * Gauss_sphere_mean2 ) * sinc
                Pq[i] = Psum
            Pq = Pq / Pq[0] # Normalize P(q) such that P(0) = 1
            return Pq

        # Calculate P(q)
        Pq = CalcPqSubFnc(q,v_mean,r, grp,hrp,jrp, grw,hrw,jrw,krw, grc,jrc,hrc,krc, gr,hr,jr,kr)

        # Calculate beta (for decoupling aproximation)
        beta_filename = PDBID + "_beta.list"
        r_beta,dB = np.genfromtxt(beta_filename, skip_header=3,usecols=[0,1],unpack=True) # import list of distance and excess scattering length values
        A00 = np.zeros(PointsInQ)
        for i in range(0,PointsInQ):
            qr = q[i]*r_beta
            A00[i] = sum(dB*np.sin(qr)/qr)
        A00 = A00 / A00[0] # normalise such that A00[0] = 1
        A002 = A00**2
        Beta = A002/Pq

        # Export P(q) and beta
        Pq_filename = filename + "_Pq.dat"
        Pq_fid = open(Pq_filename, "w")
        Pq_fid.write("q [1/AA]     P(q)         A00(q)^2     Beta(q) = A00(q)^2/P(q)\n")
        for i in range(0,PointsInQ):
            Pq_fid.write("%e %e %e %e \n" % (q[i],Pq[i],A002[i],Beta[i]))
        Pq_fid.close()
        if self.RES_button.GetValue():
            Pq_RES_filename = filename + "_Pq_RES.dat"
            Pq_RES_fid = open(Pq_RES_filename, "w")
            Pq_RES_fid.write("q [1/AA]  P(q-3dq)  P(q-2dq)  P(q-dq)  P(q)  P(q+dq)  P(q+2dq)  P(q+3dq)\n")
            Pq_MAT  = np.matrix([Pq, Pq, Pq, Pq, Pq, Pq, Pq]) # 7xPq
            for i in range(0,len(sigma)):
                q_tmp = np.array(q_MAT)[i]
                Pq_tmp = CalcPqSubFnc(q_tmp,v_mean,r, grp,hrp,jrp, grw,hrw,jrw,krw, grc,jrc,hrc,krc,gr,hr,jr,kr)
                Pq_MAT[i] = Pq_tmp
            for i in range(0,PointsInQ):
                Pq_RES_fid.write("%f %f %f %f %f %f %f %f\n" % (q[i],Pq_MAT.item(0,i),Pq_MAT.item(1,i),Pq_MAT.item(2,i),Pq_MAT.item(3,i),Pq_MAT.item(4,i),Pq_MAT.item(5,i),Pq_MAT.item(6,i)))
            Pq_RES_fid.close()

        fs = 14 # label and title fontsize
        self.PlotPqFnc(PDBID,q,Pq,fs) # plot P(q)

    ### Define funciton for plotting Pq
    def PlotPqFnc(self,PDBID,q,Pq,fs):
        if self.fitWL_button.GetValue() == 0:
            self.Figure2 = pylab.figure(2)
            Subplot2 = pylab.subplot(111)

            if self.Water_layer_button.GetValue():
                if self.Exclude_water_layer_M_button.GetValue() or self.Exclude_water_layer_OPM_button.GetValue():
                    structure_name = PDBID + " + WL " + self.Water_layer_contrast_box.GetValue()  + "% (Excl. TMD)"
                else:
                    structure_name = PDBID + " + WL " + self.Water_layer_contrast_box.GetValue()  + "%"
            else:
                structure_name = PDBID
            if self.SANS_button.GetValue():
                structure_name = structure_name + ", SANS " + self.SANS_solvent_box.GetValue() + "% D2O, [" + self.SANS_perdeut_box_A.GetValue() + "," + self.SANS_perdeut_box_B.GetValue() + "," + self.SANS_perdeut_box_C.GetValue() + "," + self.SANS_perdeut_box_D.GetValue() + "," + self.SANS_perdeut_box_E.GetValue() + "," + self.SANS_perdeut_box_F.GetValue() + "," + self.SANS_perdeut_box_G.GetValue() + "] % Deut."
            else:
                structure_name = structure_name + ", SAXS " + self.SAXS_solvent_box.GetValue() + "% Sucr"
            Subplot2.plot(q,Pq,label=structure_name)
            Subplot2.set_xscale('log',nonposx = 'clip')
            Subplot2.set_yscale('log',nonposy = 'clip')
            Subplot2.legend(loc=3,fontsize=9)
            Subplot2.set_xlabel('q [$1/\AA$]',fontsize=fs)
            Subplot2.set_ylabel('P(q)',fontsize=fs)
            pylab.suptitle('Form Factor for the structure(s)',fontsize=fs)
            pylab.show()  

    ### Define function for Fitting Bg, scale, water layer (optional parameter), and distribution of 2nd PDB (optional parameter)
    def FitFnc(self, event):
        
        # check for data file and second PDB file
        if self.DataPathStr == 'Non':
            Message = 'Please provide a data file (q,I,dI), recalculate P(q) (to obtain the right q-values) and try again.'
            print(Message)
            wx.MessageBox(Message, "CaPP - Error", wx.OK | wx.ICON_INFORMATION)
        elif self.PDB2PathStr == 'N/A' and self.fitPDB2_button.GetValue():
            Message = '"Fit with 2 PDBs" chosen. Please provide a second PDB file, using the lower "Browse" button.'
            print(Message)
            wx.MessageBox(Message, "CaPP - Error", wx.OK | wx.ICON_INFORMATION)
        else:
            # define Pq filename
            filename = self.GenerateFilenamesFnc(1)[0] # get first return value of function, filename
            Pq_filename = filename + "_Pq.dat"
            Pq_RES_filename = filename + "_Pq_RES.dat"
            if self.fitPDB2_button.GetValue():
                filename2 = self.GenerateFilenamesFnc(2)[0] # get first return value of function, filename
                Pq_filename2 = filename2 + "_Pq.dat"
                Pq_RES_filename2 = filename2 + "_Pq_RES.dat"
        
            #count number of lines in data file
            NumberOfLines = sum(1 for line in open(self.DataPathStr)) + 1

            #count number of header and footer lines in data file
            file = open(self.DataPathStr,'r')
            headerlines = 0
            datalines = 0
            footerlines = 0
            totallines = 0
            columns_in_datafile = 0
            HEADER_END = 0
            while totallines < NumberOfLines:
                line = file.readline()
                numbers_str = line.split()
                totallines = totallines + 1
                try:
                    numbers_float = [float(x) for x in numbers_str]
                    columns = len(numbers_float)
                    if columns < 3 or columns > 4:
                        if HEADER_END == 0:
                            headerlines = headerlines + 1
                        else:
                            if totallines < NumberOfLines:
                                footerlines = footerlines + 1
                    else:
                        datalines = datalines + 1 
                        columns_in_datafile = columns
                    HEADER_END = 1
                except:
                    if HEADER_END == 0:
                        headerlines = headerlines + 1
                    else:
                        footerlines = footerlines + 1
            file.close()
            if datalines == 0:
                Message = 'Cannot read datafile. Datafile should have 3 or 4 columns: q, I, and sigma_q (sigma_q is optional: SANS resolution effects).'
                print(Message)
                wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)    
            else:
                if self.RES_button.GetValue():
                    if columns_in_datafile == 3:
                        Message = 'Datafile should have 4 columns: q, I, sigma_I, sigma_q - 4th for SANS resolution effects.' 
                        print(Message)    
                        wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)    
            #import and truncate data
            q,I,dI = np.genfromtxt(self.DataPathStr,skip_header=headerlines,skip_footer=footerlines,usecols=[0,1,2], unpack=True)
            if self.nm_button.GetValue():
                q  = q *0.1 # convert from nm to angstrom
            skip_first = int(float(self.Skip_box.GetValue())) #skip first points in q
            I_trunc = I[skip_first:]
            q_trunc = q[skip_first:]
            dI_trunc = dI[skip_first:]
            
            #variables for calculation of chi2r
            Np = 2 #number of fitted parameters (bg and scale)
            if self.fitWL_button.GetValue():
                Np = Np + 1 # also fit the water layer
            if self.fitPDB2_button.GetValue():
                Np = Np + 1 # also fit the prefactor A
            df = len(q_trunc) - Np # df = degrees of freedom
            chi2r = 0.0
            
            #variables and functions for fitting
            Pq      = [1.0]*len(q_trunc) #declare P(q)
            Pq2     = Pq #declare P(q) 2
            Pq_RES  = np.matrix([Pq, Pq, Pq, Pq, Pq, Pq, Pq]) # 7xPq
            Pq2_RES = Pq_RES

            A0 = 0.5 # part of 1st PDB, initial guess
            S0 = I_trunc[2] #Scale, initial guess
            if S0 < 0.0:
                S0 = 0.0
            B0 = 0.001 #Background, initial guess
            def func1(q_trunc,S,B):
                if self.RES_button.GetValue() == 0:
                    I = S * Pq + B
                else:
                    sum_Pq_RES = [0.0]*len(q_trunc)
                    sum_w = 0.0 
                    for j in range(-3,3+1):
                        w = math.exp(-j**2/2)
                        sum_Pq_RES  += w * np.array(Pq_RES)[j+3]
                        sum_w       += w
                    PofQ = sum_Pq_RES/sum_w
                    I = S * PofQ + B
                return I
            def func2(q_trunc,S,B,A):
                if self.RES_button.GetValue() == 0:
                    I = S * (A * Pq + (1.0-A)*Pq2) + B
                else:
                    sum_Pq_RES  = [0.0]*len(q_trunc)
                    sum_Pq2_RES = [0.0]*len(q_trunc)
                    sum_w       = 0.0
                    for j in range(-3,3+1):
                        w = math.exp(-j**2/2)
                        sum_Pq_RES   += w * np.array(Pq_RES)[j+3]
                        sum_Pq2_RES  += w * np.array(Pq2_RES)[j+3]
                        sum_w        += w
                    PofQ1 = sum_Pq_RES/sum_w
                    PofQ2 = sum_Pq2_RES/sum_w
                    I = S * (A * PofQ1 + (1.0-A)*PofQ2) + B
                return I 

            #variables used if water layer contrast (wlc) is fitted - golden range search
            wlc_string = ""
            maxite = 10
            R = 0.61803398875 # inverse golden ratio
            xmin = -10
            xmax = 20
            dx = xmax-xmin
            x1 = 0.0
            x2 = 0.0
            dx = 0.0
            xfinal = 0.0
            chi2r_1 = 99.9
            chi2r_2 = 0.0
            d_chi2r = chi2r_1 - chi2r_2
            tol = 0.01 #tolerance for change in chi2r before stopping
            WhichToCalc = 0
            ii = 0
            
            # fit data
            if self.fitWL_button.GetValue():
                if self.fitPDB2_button.GetValue():
                    Message = 'Fitting WL contrast (WL), part of 1st PDB (A), scale (S) and background (B)'
                else:
                    Message = 'Fitting water layer contrast, scale and background'
                print(Message)
                while d_chi2r > tol:
                    dx = xmax-xmin
                    
                    if WhichToCalc == 1 or WhichToCalc == 0:
                        x1 = xmax-dx*R
                        wlc_string = str(x1)
                        self.Water_layer_contrast_box.SetValue(wlc_string)
                        self.cappFnc(event=None)
                        Pq = np.genfromtxt(Pq_filename, skip_header=1+skip_first,usecols=[1], unpack=True)
                        if self.RES_button.GetValue():
                            Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7 = np.genfromtxt(Pq_RES_filename,skip_header=1+skip_first,usecols=[1,2,3,4,5,6,7], unpack=True)
                            Pq_RES  = np.matrix([Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7])                        
                        if self.fitPDB2_button.GetValue():
                            Pq2 = np.genfromtxt(Pq_filename2, skip_header=1+skip_first,usecols=[1], unpack=True)
                            if self.RES_button.GetValue():
                                Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7 = np.genfromtxt(Pq_RES_filename2, skip_header=1+skip_first,usecols=[1,2,3,4,5,6,7], unpack=True)
                                Pq2_RES  = np.matrix([Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7])                              
                            popt, pcov = curve_fit(func2,q_trunc,I_trunc,bounds=((0,-np.inf,0),(np.inf,np.inf,1)),sigma=dI_trunc,p0=[S0,B0,A0])    
                            I_fit = func2(q_trunc,*popt)
                        else:
                            popt, pcov = curve_fit(func1,q_trunc,I_trunc,sigma=dI_trunc,p0=[S0,B0])
                            I_fit = func1(q_trunc,*popt)
                        chi2r_1 = np.sum(((I_trunc - I_fit )/dI_trunc)**2)/df
                        ii = ii + 1
                        Status = " Ite = %d, Chi2r = %f, WL = %f" % (ii, chi2r_1, x1)
                        print(Status)
                
                    if WhichToCalc == 2 or WhichToCalc == 0:
                        x2 = dx*R+xmin
                        wlc_string = str(x2)
                        self.Water_layer_contrast_box.SetValue(wlc_string)
                        self.cappFnc(event=None)
                        Pq = np.genfromtxt(Pq_filename, skip_header=1+skip_first,usecols=[1], unpack=True)
                        if self.RES_button.GetValue():
                            Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7 = np.genfromtxt(Pq_RES_filename, skip_header=1+skip_first,usecols=[1,2,3,4,5,6,7], unpack=True)
                            Pq_RES  = np.matrix([Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7])
                        if self.fitPDB2_button.GetValue():
                            Pq2 = np.genfromtxt(Pq_filename2, skip_header=1+skip_first,usecols=[1], unpack=True)
                            if self.RES_button.GetValue():
                                Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7 = np.genfromtxt(Pq_RES_filename2, skip_header=1+skip_first,usecols=[1,2,3,4,5,6,7], unpack=True)
                                Pq2_RES  = np.matrix([Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7])                
                            popt, pcov = curve_fit(func2,q_trunc,I_trunc,bounds=((0,-np.inf,0),(np.inf,np.inf,1)),sigma=dI_trunc,p0=[S0,B0,A0])
                            I_fit = func2(q_trunc,*popt) 
                        else:
                            popt, pcov = curve_fit(func1,q_trunc,I_trunc,sigma=dI_trunc,p0=[S0,B0])
                            I_fit = func1(q_trunc,*popt)
                        chi2r_2 = np.sum(((I_trunc - I_fit )/dI_trunc)**2)/df
                        ii = ii + 1
                        Status = " Ite = %d, Chi2r = %f, WL = %f" % (ii, chi2r_2, x2)
                        print(Status)
                    if chi2r_1 > chi2r_2:
                        xfinal = (x1 + x2)/2
                        xmin = x1
                        d_chi2r = chi2r_1 - chi2r_2
                        x1 = x2
                        chi2r_1 = chi2r_2
                        WhichToCalc = 2
                    elif chi2r_1 < chi2r_2:
                        xfinal = (x1 + x2)/2
                        xmax = x2
                        d_chi2r = chi2r_2 - chi2r_1
                        x2 = x1
                        chi2r_2 = chi2r_1
                        WhichToCalc = 1
                    else:
                        d_chi2r = 0.0
                    if d_chi2r <= tol:
                        Message=' Fit has converged'
                        #wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)
                        #print(Message)
                    elif ii == maxite :
                        d_chi2r = 0.0
                        Message = ' Max iterations reached. Algorithm stopped'
                        #wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)
                        #print(Message)
                wlc_string = str(xfinal)
                self.Water_layer_contrast_box.SetValue(wlc_string)
                self.cappFnc(event=None)
                Pq = np.genfromtxt(Pq_filename, skip_header=1+skip_first,usecols=[1], unpack=True)
                if self.RES_button.GetValue():
                    Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7 = np.genfromtxt(Pq_RES_filename, skip_header=1+skip_first,usecols=[1,2,3,4,5,6,7], unpack=True)
                    Pq_RES  = np.matrix([Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7])
                if self.fitPDB2_button.GetValue():
                    Pq2 = np.genfromtxt(Pq_filename2, skip_header=1+skip_first,usecols=[1], unpack=True)
                    if self.RES_button.GetValue():
                        Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7 = np.genfromtxt(Pq_RES_filename2, skip_header=1+skip_first,usecols=[1,2,3,4,5,6,7], unpack=True)
                        Pq2_RES  = np.matrix([Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7])                
                    popt, pcov = curve_fit(func2,q_trunc,I_trunc,bounds=((0,-np.inf,0),(np.inf,np.inf,1)),sigma=dI_trunc,p0=[S0,B0,A0])
                    I_fit = func2(q_trunc,*popt) 
                else:
                    popt, pcov = curve_fit(func1,q_trunc,I_trunc,sigma=dI_trunc,p0=[S0,B0])
                    I_fit = func1(q_trunc,*popt)
                chi2r = np.sum(( (I_trunc - I_fit) / dI_trunc)**2)/df
                ii = ii + 1
                Status = " Ite = %d, Chi2r = %f, WL = %f" % (ii, chi2r, xfinal)
                print(Status)
                wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)
                print(Message)
            else:
                self.cappFnc(event=None)
                Pq = np.genfromtxt(Pq_filename, skip_header=1+skip_first,usecols=[1], unpack=True)
                if self.RES_button.GetValue():
                    Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7 = np.genfromtxt(Pq_RES_filename, skip_header=1+skip_first,usecols=[1,2,3,4,5,6,7], unpack=True)
                    Pq_RES  = np.matrix([Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7])
                if self.fitPDB2_button.GetValue():
                    Pq2 = np.genfromtxt(Pq_filename2, skip_header=1+skip_first,usecols=[1], unpack=True)
                    if self.RES_button.GetValue():
                        Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7 = np.genfromtxt(Pq_RES_filename2, skip_header=1+skip_first,usecols=[1,2,3,4,5,6,7], unpack=True)
                        Pq2_RES  = np.matrix([Pq_R1,Pq_R2,Pq_R3,Pq_R4,Pq_R5,Pq_R6,Pq_R7])         
                    popt, pcov = curve_fit(func2,q_trunc,I_trunc,bounds=((0,-np.inf,0),(np.inf,np.inf,1)),sigma=dI_trunc,p0=[S0,B0,A0])
                    I_fit = func2(q_trunc,*popt) 
                else:
                    popt, pcov = curve_fit(func1,q_trunc,I_trunc,sigma=dI_trunc,p0=[S0,B0])
                    I_fit = func1(q_trunc,*popt)
                chi2r = np.sum(( (I_trunc - I_fit) / dI_trunc)**2)/df
            
            # plot data and fit
            wlc = float(self.Water_layer_contrast_box.GetValue()) #water layer contrast
            if matplotlib.pyplot.fignum_exists(3):
                self.Figure3 = pylab.figure(3)
            else:
                self.Figure3 = pylab.figure(3)
                Subplot3 = pylab.subplot(111)
                Subplot3.errorbar(q, I, dI, label="Data", color='r',fmt='.',zorder=0)
            Subplot3 = pylab.subplot(111)
            if self.fitPDB2_button.GetValue():
                if self.Water_layer_button.GetValue():
                    fitvalues = 'fit:    Sc=%5.2f, Bg=%5.3f, A = %5.2f, WL=%4.1f, Chi2r=%4.2f' % (popt[0],popt[1],popt[2],wlc,chi2r)
                    Subplot3.plot(q_trunc, I_fit, label=fitvalues)
                else:
                    fitvalues = 'fit: Sc=%5.2f, Bg=%5.3f, A = %5.2f, Chi2r=%4.2f' % (popt[0],popt[1],popt[2],chi2r)
                    Subplot3.plot(q_trunc, I_fit, label=fitvalues)
                fiterrors = 'errors(std. dev): dSc=%5.2f, dBg=%5.3f, dA = %5.2f' % (math.sqrt(pcov.item(0,0)),math.sqrt(pcov.item(1,1)),math.sqrt(pcov.item(2,2)))
            else:
                if self.Water_layer_button.GetValue():
                    fitvalues = 'fit: Scale=%5.2f, Bg=%5.3f, WL=%4.1f, Chi2r=%4.2f' % (popt[0],popt[1],wlc,chi2r)
                    Subplot3.plot(q_trunc, I_fit, label=fitvalues)
                else:
                    fitvalues ='fit: Scale=%5.2f, Bg=%5.3f, Chi2r=%4.2f' % (popt[0],popt[1],chi2r)
                    Subplot3.plot(q_trunc, I_fit, label=fitvalues)
                fiterrors = 'errors(std. dev): dSc=%5.2f, dBg=%5.3f' % (math.sqrt(pcov.item(0,0)),math.sqrt(pcov.item(1,1)))
            Subplot3.set_xscale('log', nonposx = 'clip')
            Subplot3.set_yscale('log', nonposy = 'clip')
            Subplot3.legend(loc=3,fontsize=9)
            Subplot3.set_xlabel('q [$1/\AA$]',fontsize=14)
            Subplot3.set_ylabel('I(q)',fontsize=14)
            pylab.suptitle('Fit structure to data',fontsize=14)
            pylab.show()
            if chi2r > 100:
                Message = 'Very poor fit! (Chi2r>100) - is the data right? units of data in aa or nm? is the structure right? oligomerization? unspecific aggregation?'
                print(Message)
                wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)

            # export fit
            filename = os.path.splitext(self.PDBPathStr)[0]

            if self.Water_layer_button.GetValue():
                filename = filename + "_w"

            if self.fitPDB2_button.GetValue():
                filename2 = os.path.splitext(self.PDB2PathStr)[0]
                filename2 = self.path_leaf(filename2)
                filename = filename + "_" + filename2
            if self.SANS_button.GetValue():
                value = float(self.SANS_solvent_box.GetValue())
                string_value = "%1.0f" % value
                filename = filename + "_N" + string_value
            else:
                value = float(self.SAXS_solvent_box.GetValue())
                string_value = "%1.0f" % value
                filename = filename + "_X" + string_value
            
            fit_filename = filename + "_Fit.dat"

            Fit_fid = open(fit_filename, "w")
            Fit_fid.write("#   %s\n" % fitvalues)
            Fit_fid.write("#   %s\n" % fiterrors)
            if self.Water_layer_button.GetValue() and self.fitPDB2_button.GetValue():
                Fit_fid.write("#   Sc = Scale, Bg = Background, WL = excees water layer scattering density, A = part of 1st PDB, Chi2r = reduced chi square\n")
            elif self.Water_layer_button.GetValue():
                Fit_fid.write("#   Sc = Scale, Bg = Background, WL = excees water layer scattering density, Chi2r = reduced chi square\n")
            elif self.fitPDB2_button.GetValue():
                Fit_fid.write("#   Sc = Scale, Bg = Background, A = part of 1st PDB, Chi2r = reduced chi square\n")
            else:
                Fit_fid.write("#   Sc = Scale, Bg = Background, Chi2r = reduced chi square\n")

            if self.SANS_button.GetValue():
                Conditions = "SANS in " + self.SANS_solvent_box.GetValue() + "% D2O"
            else:
                Conditions = "SAXS in " + self.SAXS_solvent_box.GetValue() + "% Sucrose"

            if self.Water_layer_button.GetValue():
                if self.Exclude_water_layer_M_button.GetValue() or self.Exclude_water_layer_OPM_button.GetValue():
                    Conditions = Conditions + ", WL added with " + self.Water_layer_contrast_box.GetValue()  + "% contrast (Excluded at TMD)."
                else:
                    Conditions = Conditions + ", WL added with " + self.Water_layer_contrast_box.GetValue()  + "% contrast."

            if self.fitPDB2_button.GetValue():
                Fit_fid.write("#   PDB 1: %s\n" % self.PDBPathStr)
                Fit_fid.write("#   PDB 2: %s\n" % self.PDB2PathStr)
            else:
                Fit_fid.write("#   PDB: %s\n" % self.PDBPathStr)
            Fit_fid.write("#   Data: %s\n" % self.DataPathStr)
            Fit_fid.write("#   Conditions/Settings: %s\n" % Conditions)
            Fit_fid.write("#   q [1/AA]  \tData_I \tData_dI  \tFit_I\n")

            for i in range(0,len(q_trunc)):
                Fit_fid.write("     %e \t%e \t%e \t%e\n" % (q_trunc[i],I_trunc[i],dI_trunc[i],I_fit[i]))
            Fit_fid.close()

### Boilerplate code linking program and widgets
if __name__ == '__main__':
    App = wx.App(False)
    frame = MainCls(parent = None, id = -1)
    frame.Show()
    App.MainLoop()
