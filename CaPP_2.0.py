##        CaPP 1.0              ##
##                              ##
##    Copyright 2017,           ##
##    Andreas Larsen            ##
##    andreas.larsen@nbi.dk     ##

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

## Libraries
import re
import os
import sys
import math
import linecache
import re
import platform
try:
    import wx
except:
    print("")
    print("******************************************************************")
    print(" CaPP failed to import wxPython - is it correctly installed...?   ")
    print("******************************************************************")
    print("")
    sys.exit(1)

if platform.system() == "Darwin":
    capp_version = "capp_mac"
elif platform.system() == "Windows":
    capp_version = "capp_windows.exe"
else:
    capp_version = "capp"

if os.path.exists(capp_version):
    print("CaPP initiated correctly")
elif capp_version == "capp":
    print("")
    print("****************************************************************************************")
    print(" Cannot find the executable: %s. It should be in the same folder as CaPP.py " % capp_version)
    print(" If you are not running Mac or Windows OS, then:")
    print("  1) compile the source code Mainfunction.c with output called capp")
    print("     e.g. gcc Mainfunction.c -o capp")
    print("  2) plase the executable, capp, in the same folder as CaPP.py")
    print("  3) re-run CaPP with: python CaPP.py")
    print("****************************************************************************************")
    print("")
    sys.exit(1)

else:
    print("")
    print("****************************************************************************************")
    print(" Cannot find the executable: %s. It should be in the same folder as CaPP.py " % capp_version)
    print("****************************************************************************************")
    print("")
    sys.exit(1)

import wx.lib.scrolledpanel

programpath = os.getcwd()

## Import plotting libraries
import matplotlib
matplotlib.interactive(True)
matplotlib.use('WXAgg')
import pylab
import numpy as np

## import fitting libraries
from scipy.optimize import curve_fit

## Define main class and text
class MainCls(wx.Frame):
    def __init__(self, parent, id):


        ### Overall frame for widgets
        width = wx.SystemSettings.GetMetric(wx.SYS_SCREEN_X)
        height = wx.SystemSettings.GetMetric(wx.SYS_SCREEN_Y)
        position = (width/3.0, height/15.0)
        wx.Frame.__init__(self, parent, id, title = 'CaPP - Calculate p(r) from PDB files', pos = position)
        BoxSizer = wx.BoxSizer(wx.VERTICAL)
        self.Panel = wx.lib.scrolledpanel.ScrolledPanel(self, -1)
        self.Panel.SetupScrolling(False, True)
        self.Bind(wx.EVT_CLOSE, self.CloseWindowFnc) # close pylab explicitely for to avoid crash


        ### Widgets for PDB-file import
        BoxSizer.AddSpacer(10) #vertical spacing

        PDBTextSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDBTextSizer.AddSpacer(10)

        PDBText = wx.StaticText(self.Panel, -1, 'Location of PDB file:', size = (330, -1))
        PDBTextSizer.Add(PDBText)

        BoxSizer.Add(PDBTextSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(10)

        PDBPathSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDBPathSizer.AddSpacer(10)

        self.PDBPathTxt = wx.StaticText(self.Panel, -1, '')
        self.PDBPathStr = 'N/A'
        PDBPathSizer.Add(self.PDBPathTxt)

        BoxSizer.Add(PDBPathSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(10)

        PDBBtnSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDBBtnSizer.AddSpacer(10)

        BrowsePDBBtn = wx.Button(self.Panel, label = 'Browse')
        self.Bind(wx.EVT_BUTTON, self.BrowsePDBFnc, BrowsePDBBtn)
        PDBBtnSizer.Add(BrowsePDBBtn, 1, wx.EXPAND)
        PDBBtnSizer.AddSpacer(10)

        BoxSizer.Add(PDBBtnSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(10)

        LinePDB = wx.StaticLine(self.Panel, -1)
        BoxSizer.Add(LinePDB, 0, wx.EXPAND|wx.HORIZONTAL)
        

        ### Widgets for Data-file import
        BoxSizer.AddSpacer(10) # vertical and horizontal spacing
    
        DataTextSizer = wx.BoxSizer(wx.HORIZONTAL)
        DataTextSizer.AddSpacer(10)
        
        DataText = wx.StaticText(self.Panel, -1, 'Location of Data file (optional, for q-range and fit):', size = (330, -1))
        DataTextSizer.Add(DataText)
        
        BoxSizer.Add(DataTextSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(10)
        
        DataPathSizer = wx.BoxSizer(wx.HORIZONTAL)
        DataPathSizer.AddSpacer(10)
        
        self.DataPathTxt = wx.StaticText(self.Panel, -1, '')
        self.DataPathStr = 'Non'
        DataPathSizer.Add(self.DataPathTxt)
        
        BoxSizer.Add(DataPathSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(10)
        
        DataBtnSizer = wx.BoxSizer(wx.HORIZONTAL)
        DataBtnSizer.AddSpacer(10)
        
        BrowseDataBtn = wx.Button(self.Panel, label = 'Browse')
        self.Bind(wx.EVT_BUTTON, self.BrowseDataFnc, BrowseDataBtn)
        DataBtnSizer.Add(BrowseDataBtn, 1, wx.EXPAND)
        DataBtnSizer.AddSpacer(10)
        
        BoxSizer.Add(DataBtnSizer, 0, wx.EXPAND|wx.HORIZONTAL)


        ### Widgets used to choose units of data
        BoxSizer.AddSpacer(10)
    
        nm_button = wx.BoxSizer(wx.HORIZONTAL)
        nm_button.AddSpacer(10)
        self.nm_button = wx.CheckBox(self.Panel, -1, 'q in data is in nm? (default unit is angstrom)', size = (350, -1))
        nm_button.Add(self.nm_button)
        BoxSizer.Add(nm_button, 0, wx.ALIGN_CENTER, wx.EXPAND|wx.HORIZONTAL)
        self.Bind(wx.EVT_CHECKBOX, self.nmFnc, self.nm_button)
        BoxSizer.AddSpacer(10)
    
        LineData = wx.StaticLine(self.Panel, -1)
        BoxSizer.Add(LineData, 0, wx.EXPAND|wx.HORIZONTAL)


        ### Widgets used to choose between SAXS and SANS
        BoxSizer.AddSpacer(10)

        SAXS_button = wx.BoxSizer(wx.HORIZONTAL) # prepare for SAXS button
        SAXS_button.AddSpacer(10)
        self.SAXS_button = wx.RadioButton(self.Panel, -1, 'SAXS',style=wx.RB_GROUP) # make SAXS button
        SAXS_button.Add(self.SAXS_button)
        BoxSizer.Add(SAXS_button, 0, wx.EXPAND|wx.HORIZONTAL) # position SAXS button
        self.Bind(wx.EVT_RADIOBUTTON, self.DisableSANSFnc, self.SAXS_button)
        

        BoxSizer.AddSpacer(10)

        SANS_button = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS button
        SANS_button.AddSpacer(10) # vertical spacing
        self.SANS_button = wx.RadioButton(self.Panel, -1, 'SANS') # make SANS button
        SANS_button.Add(self.SANS_button) # ?
        BoxSizer.Add(SANS_button, 0, wx.EXPAND|wx.HORIZONTAL) # position SANS button
        self.Bind(wx.EVT_RADIOBUTTON, self.EnableSANSFnc, self.SANS_button)

        SANS_solvent_box = wx.BoxSizer(wx.HORIZONTAL) # prepare for SANS solvent box
        SANS_solvent_box.AddSpacer(29)

        self.left_SANS = wx.StaticText(self.Panel, -1, "D2O content =")
        SANS_solvent_box.Add(self.left_SANS)
        self.left_SANS.Disable()

        self.SANS_solvent_box = wx.TextCtrl(self.Panel, -1, '', size = (40, -1)) # make SANS solvent box
        self.SANS_solvent_box.SetValue('100')
        SANS_solvent_box.Add(self.SANS_solvent_box) # ?
        BoxSizer.Add(SANS_solvent_box) # position SANS solvent box
        self.SANS_solvent_box.Disable()

        self.right_SANS = wx.StaticText(self.Panel, -1, "%")
        SANS_solvent_box.Add(self.right_SANS)
        self.right_SANS.Disable()

        Line_after_contrast_widget = wx.StaticLine(self.Panel, -1) # make line after widget
        BoxSizer.Add(Line_after_contrast_widget, 0, wx.EXPAND|wx.HORIZONTAL)


        ### Widgets used to add Water Layer
        BoxSizer.AddSpacer(10)

        Water_layer_button = wx.BoxSizer(wx.HORIZONTAL)
        Water_layer_button.AddSpacer(10)
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
        self.Water_layer_contrast_box.SetValue('10')
        Water_layer_contrast_box.Add(self.Water_layer_contrast_box)
        BoxSizer.Add(Water_layer_contrast_box)
        self.Water_layer_contrast_box.Disable()

        self.right_WL = wx.StaticText(self.Panel, -1, "% of solvent scattering")
        Water_layer_contrast_box.Add(self.right_WL)
        self.right_WL.Disable()


        ### Widgets to exclude water layer from TMD
        BoxSizer.AddSpacer(10)

        WLText = wx.BoxSizer(wx.HORIZONTAL)
        WLText.AddSpacer(10)
        self.WLText = wx.StaticText(self.Panel, -1, 'If the protein has a Transmembrane Domain (TMD), consider:', size = (330, -1))
        WLText.Add(self.WLText)
        BoxSizer.Add(WLText, 0, wx.EXPAND|wx.HORIZONTAL)
        self.WLText.Disable()

        BoxSizer.AddSpacer(10)

        No_Exclude_water_layer_button = wx.BoxSizer(wx.HORIZONTAL)
        No_Exclude_water_layer_button.AddSpacer(10)
        self.No_Exclude_water_layer_button = wx.RadioButton(self.Panel, -1, 'Do not exclude the WL from the TMD', style=wx.RB_GROUP)
        No_Exclude_water_layer_button.Add(self.No_Exclude_water_layer_button)
        BoxSizer.Add(No_Exclude_water_layer_button, 0, wx.EXPAND|wx.HORIZONTAL)
        self.No_Exclude_water_layer_button.Disable()
        self.Bind(wx.EVT_RADIOBUTTON, self.DisableBilayerThicknessFnc, self.No_Exclude_water_layer_button)

        BoxSizer.AddSpacer(10)

        Exclude_water_layer_OPM_button = wx.BoxSizer(wx.HORIZONTAL)
        Exclude_water_layer_OPM_button.AddSpacer(10)
        self.Exclude_water_layer_OPM_button = wx.RadioButton(self.Panel, -1, 'Exclude the WL from the TMD, using OPM(*)')
        Exclude_water_layer_OPM_button.Add(self.Exclude_water_layer_OPM_button)
        BoxSizer.Add(Exclude_water_layer_OPM_button, 0, wx.EXPAND|wx.HORIZONTAL)
        self.Exclude_water_layer_OPM_button.Disable()
        self.Bind(wx.EVT_RADIOBUTTON, self.DisableBilayerThicknessFnc, self.Exclude_water_layer_OPM_button)

        BoxSizer.AddSpacer(10)

        Exclude_water_layer_M_button = wx.BoxSizer(wx.HORIZONTAL)
        Exclude_water_layer_M_button.AddSpacer(10)
        self.Exclude_water_layer_M_button = wx.RadioButton(self.Panel, -1, 'Exclude WL from TMD, manually:')
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
        self.OPMText = wx.StaticText(self.Panel, -1, '(*) The Orientation of Proteins in Membranes (OPM) Database', size = (330, -1))
        OPMText.Add(self.OPMText)
        BoxSizer.Add(OPMText, 0, wx.EXPAND|wx.HORIZONTAL)
        self.OPMText.Disable()

        BoxSizer.AddSpacer(2)

        Line_after_exclude_WL_widget = wx.StaticLine(self.Panel, -1) # make line after widget
        BoxSizer.Add(Line_after_exclude_WL_widget, 0, wx.EXPAND|wx.HORIZONTAL)


        ### Widgets to change resolution
        BoxSizer.AddSpacer(10)

        Change_resolution_button = wx.BoxSizer(wx.HORIZONTAL)
        Change_resolution_button.AddSpacer(10)
        self.Change_resolution_button = wx.CheckBox(self.Panel, -1, 'Change resolution of the p(r)? (Default: 3 Angstrom)', size = (350, -1))
        Change_resolution_button.Add(self.Change_resolution_button)
        BoxSizer.Add(Change_resolution_button, 0, wx.ALIGN_CENTER, wx.EXPAND|wx.HORIZONTAL)
        self.Bind(wx.EVT_CHECKBOX, self.EnableResolutionBoxFnc, self.Change_resolution_button)

        Resolution_box = wx.BoxSizer(wx.HORIZONTAL)
        Resolution_box.AddSpacer(29)

        self.left_res = wx.StaticText(self.Panel, -1, "Binsize =")
        Resolution_box.Add(self.left_res)
        self.left_res.Disable()

        self.Resolution_box = wx.TextCtrl(self.Panel, -1, '', size = (40, -1))
        self.Resolution_box.SetValue('2.0')
        Resolution_box.Add(self.Resolution_box)
        BoxSizer.Add(Resolution_box)
        self.Resolution_box.Disable()

        self.right_res = wx.StaticText(self.Panel, -1, "Angstrom")
        Resolution_box.Add(self.right_res)
        self.right_res.Disable()

        Line_after_resolution_widget = wx.StaticLine(self.Panel, -1) # make line after widget
        BoxSizer.Add(Line_after_resolution_widget, 0, wx.EXPAND|wx.HORIZONTAL)


        ### Widgets for calculation buttons
        BoxSizer.AddSpacer(10)

        CalculateRgButtonSpace = wx.BoxSizer(wx.HORIZONTAL)
        CalculateRgButtonSpace.AddSpacer(10)
        self.CalculateRgButton = wx.Button(self.Panel, label = 'Calculate Rg')
        self.Bind(wx.EVT_BUTTON, self.RgFnc, self.CalculateRgButton)
        self.CalculateRgButton.Disable()
        CalculateRgButtonSpace.Add(self.CalculateRgButton, 1, wx.EXPAND)
        CalculateRgButtonSpace.AddSpacer(10)
        BoxSizer.Add(CalculateRgButtonSpace, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(1)

        CalculateButtonSpace = wx.BoxSizer(wx.HORIZONTAL)
        CalculateButtonSpace.AddSpacer(10)
        self.CalculateButton = wx.Button(self.Panel, label = 'Calculate p(r)')
        self.Bind(wx.EVT_BUTTON, self.cappFnc, self.CalculateButton)
        self.CalculateButton.Disable()
        CalculateButtonSpace.Add(self.CalculateButton, 1, wx.EXPAND)
        CalculateButtonSpace.AddSpacer(10)
        BoxSizer.Add(CalculateButtonSpace, 0, wx.EXPAND|wx.HORIZONTAL)
        BoxSizer.AddSpacer(1)

        CalcPqButtonSpace = wx.BoxSizer(wx.HORIZONTAL)
        CalcPqButtonSpace.AddSpacer(10)
        self.CalcPqButton = wx.Button(self.Panel, label = 'Calculate P(q) from p(r)')
        self.Bind(wx.EVT_BUTTON, self.CalcPqFnc, self.CalcPqButton)
        CalcPqButtonSpace.Add(self.CalcPqButton, 1, wx.EXPAND)
        CalcPqButtonSpace.AddSpacer(10)
        BoxSizer.Add(CalcPqButtonSpace, 0, wx.EXPAND|wx.HORIZONTAL)
        self.CalcPqButton.Disable()
        BoxSizer.AddSpacer(1)

        FitButtonSpace = wx.BoxSizer(wx.HORIZONTAL)
        FitButtonSpace.AddSpacer(10)
        self.FitButton = wx.Button(self.Panel, label = 'Fit P(q) to data')
        self.Bind(wx.EVT_BUTTON, self.FitFnc, self.FitButton)
        FitButtonSpace.Add(self.FitButton, 1, wx.EXPAND)
        FitButtonSpace.AddSpacer(10)
        BoxSizer.Add(FitButtonSpace, 0, wx.EXPAND|wx.HORIZONTAL)
        self.FitButton.Disable()
        BoxSizer.AddSpacer(10)

        ### Conclusion of init-function
        self.Panel.SetSizerAndFit(BoxSizer)
        self.SetSizerAndFit(BoxSizer)
        self.SetAutoLayout(True)
        self.Panel.Layout()
        self.Layout()


### Define minor functions:
#       "onClose()",
#       "nmFnc()"
#       "EnableWLButtons()",
#       "EnableResolutionBoxFnc()",
#       "EnableBilayerThicknessFnc()",
#       "DisableBilayerThicknessFnc()",
#       "EnableSANSFnc()",
#       "DisableSANSFnc()",

    def onClose(self,event):
        self.Destroy()

    def nmFnc(self,event):
        if self.Water_layer_button.GetValue():
            print("data in nm")

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
            self.Exclude_water_layer_M_button.SetValue(False)
            self.Exclude_water_layer_OPM_button.SetValue(False)
            self.No_Exclude_water_layer_button.SetValue(False)

    def EnableResolutionBoxFnc(self, event):
        if self.Change_resolution_button.GetValue():
            self.left_res.Enable()
            self.Resolution_box.Enable()
            self.right_res.Enable()
        else:
            self.left_res.Disable()
            self.Resolution_box.Disable()
            self.right_res.Disable()

    def EnableBilayerThicknessFnc(self, event):
        self.Bilayer_thickness_box.Enable()
        self.right_thick.Enable()
        self.left_thick.Enable()

    def DisableBilayerThicknessFnc(self, event):
        self.Bilayer_thickness_box.Disable()
        self.right_thick.Disable()
        self.left_thick.Disable()

    def EnableSANSFnc(self, event):
        self.left_SANS.Enable()
        self.SANS_solvent_box.Enable()
        self.right_SANS.Enable()

    def DisableSANSFnc(self, event):
        self.left_SANS.Disable()
        self.SANS_solvent_box.Disable()
        self.right_SANS.Disable()

    def RgFnc(self, event):
        Output_Rg = "%s/%s -g %s" % (programpath, capp_version, self.PDBPathStr)
        os.system(Output_Rg)
        Message = 'See the Rg of the protein in the terminal window'
        wx.MessageBox(Message, "Radius of gyration", wx.OK | wx.ICON_INFORMATION)

    ### Define RgFnc function
    def RgFnc(self, event):
        
        print(os.getcwd())
        
        # import options from GUI
        XN_choice           = ""
        solvent             = ""
        WL_choice           = ""
        WL_contrast         = ""
        Exclude_WL_choice   = ""
        Bilayer_thickness   = ""
        Resolution_choice   = ""
        Resolution          = ""
        
        if self.Water_layer_button.GetValue():
            WL_choice               = "-c"
            WL_contrast_tmp         = float(self.Water_layer_contrast_box.GetValue())/100
            WL_contrast             = "%1.2f" % WL_contrast_tmp
        
        if self.Exclude_water_layer_OPM_button.GetValue():
            Exclude_WL_choice       = "-d"
        if self.Exclude_water_layer_M_button.GetValue():
            Exclude_WL_choice       = "-m"
            Bilayer_thickness       = self.Bilayer_thickness_box.GetValue()
    
        if self.Change_resolution_button.GetValue():
            Resolution_choice   = "-r"
            Resolution          = self.Resolution_box.GetValue()
        
        if self.SAXS_button.GetValue():
            dummy = 1 # do nothing
        elif self.SANS_button.GetValue():
            XN_choice = "-s"
            solvent_tmp = float(self.SANS_solvent_box.GetValue())/100
            solvent = "%1.3f" % solvent_tmp
        else:
            Message = 'SAXS chosen as default'
            wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)

        # assemble options into a command line
        Output = "%s/%s %s %s %s %s %s %s %s %s -g %s" % (programpath, capp_version, XN_choice, solvent, WL_choice, WL_contrast, Exclude_WL_choice, Bilayer_thickness, Resolution_choice, Resolution, self.PDBPathStr)

        # run command line
        os.system(Output)
        
        # print command line to terminal window
        print("Command generated by the GUI and sent to terminal:")
        print(Output)
        print("")
    
    ### Define CaPP function
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
        WL_choice           = ""
        WL_contrast         = ""
        Exclude_WL_choice   = ""
        Bilayer_thickness   = ""
        Resolution_choice   = ""
        Resolution          = ""

        if self.Water_layer_button.GetValue():
            WL_choice               = "-c"
            WL_contrast_tmp         = float(self.Water_layer_contrast_box.GetValue())/100
            WL_contrast             = "%1.2f" % WL_contrast_tmp

        if self.Exclude_water_layer_OPM_button.GetValue():
            #check file
            file = fopen(self.PDBPathStr,'r')
            line = file.readline()
            
            
            Message = 'ERROR: see terminal window'
            wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)
            Exclude_WL_choice       = "-d"
        if self.Exclude_water_layer_M_button.GetValue():
            Exclude_WL_choice       = "-m"
            Bilayer_thickness       = self.Bilayer_thickness_box.GetValue()

        if self.Change_resolution_button.GetValue():
                Resolution_choice   = "-r"
                Resolution          = self.Resolution_box.GetValue()

        if self.SAXS_button.GetValue():
            dummy = 1 # do nothing
        elif self.SANS_button.GetValue():
            XN_choice = "-s"
            solvent_tmp = float(self.SANS_solvent_box.GetValue())/100
            solvent = "%1.3f" % solvent_tmp
        else:
            Message = 'SAXS chosen as default'
            wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)
    
        # assemble options into a command line
        Output = "%s/%s %s %s %s %s %s %s %s %s %s" % (programpath, capp_version, XN_choice, solvent, WL_choice, WL_contrast, Exclude_WL_choice, Bilayer_thickness, Resolution_choice, Resolution, self.PDBPathStr)

        # run command line
        os.system(Output)

        # close pop-up window
        self.second_window.Destroy()

        # print command line to terminal window
        print("Command generated by the GUI and sent to terminal:")
        print(Output)
        print("")

        # Plot p(r)
        if capp_version == "capp_windows.exe":
            print("Plotting not available on Windows") # then it crashes
            Message = 'p(r) succesfully calculated - see folder with pdb file'
        else:
            # import p(r)
            filename, extension = os.path.splitext(self.PDBPathStr)
            if self.Water_layer_button.GetValue():
                pr_filename = filename + "_w_pr.dat"
            else:
                pr_filename = filename + "_pr.dat"
            r,pr = np.loadtxt(pr_filename, skiprows=10,usecols=[0,1],unpack=True)
            pr = pr / max(pr)

            # plot p(r)
            self.Figure1 = pylab.figure(1)
            Subplot = pylab.subplot(111)
            length_of_name = len(pr_filename)
            if self.Water_layer_button.GetValue():
                if self.Exclude_water_layer_M_button.GetValue() or self.Exclude_water_layer_OPM_button.GetValue():
                    structure_name = "'..." + pr_filename[length_of_name-10-9:-9] + "'" + " + WL " + self.Water_layer_contrast_box.GetValue()  + "% (Excl. TMD)"
                else:
                    structure_name = "'..." + pr_filename[length_of_name-10-9:-9] + "'" + " + WL " + self.Water_layer_contrast_box.GetValue()  + "%"
            else:
                structure_name = "'..." + pr_filename[length_of_name-10-7:-7] + "'"

            if self.SANS_button.GetValue():
                structure_name = structure_name + ", SANS " + self.SANS_solvent_box.GetValue() + "% D2O"
            else:
                structure_name = structure_name + ", SAXS "
            Subplot.plot(r, pr, label=structure_name)
            Subplot.legend(fontsize=9)
            Subplot.set_xlabel('r [$\AA$]',fontsize=14)
            Subplot.set_ylabel('p(r)',fontsize=14)
            pylab.suptitle('Pair Distance Distribution for the structure(s)',fontsize=14)
            pylab.show()

        # enable P(q) calculation button
        self.CalcPqButton.Enable()


    ### Define function for calculating form factor P(q)
    def CalcPqFnc(self, event):

        # import p(r) from file
        filename, extension = os.path.splitext(self.PDBPathStr)
        if self.Water_layer_button.GetValue():
            pr_filename = filename + "_w_pr.dat"
        else:
            pr_filename = filename + "_pr.dat"
        r,pr,gr,hr,jr,kr = np.loadtxt(pr_filename, skiprows=10,usecols=[0,1,2,3,4,5],unpack=True)
        pr = pr / max(pr)

        # import mean volume
        v_mean_line = linecache.getline(pr_filename,8)
        v_mean = float(v_mean_line[16:])

        PointsInr = len(r)
        
        # create or import q-vector
        if self.DataPathStr == 'Non':
            # create q-vector
            PointsInQ = 200
            qmin = 0.001
            qmax = 1.0
            dq = (qmax - qmin) / (PointsInQ - 1)
            q = np.arange(qmin, qmax+dq, dq)
        else:
            #import q-vector
            headerlines = 2
            q = np.genfromtxt(self.DataPathStr, skiprows=headerlines,usecols=[0], unpack=True)
            if self.nm_button.GetValue():
                q = q*0.1 # convert from nm to angstrom
        
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
        Pq = np.zeros(PointsInQ)
        inv4pi = 1.0/(4.0*math.pi)
        inv4pi2 = inv4pi*inv4pi
        qr = 0.0
        for i in range(0,PointsInQ):
            q2 = q[i]*q[i]
            if self.SANS_button.GetValue():
                psi_c = 1.0;
            else:
                psi_sum = c;
                for j in range(0,4):
                    psi_sum = psi_sum + a[j] * math.exp(-b[j] * q2 * inv4pi2)
                psi_c = psi_sum/b_c # normalize to psi_c(q = 0) = 1. f_c = psi_c * b_c, where b_c = 6.0
            psi_c2 = psi_c * psi_c
            Gauss_sphere_mean = math.exp(- q2 * v_mean**(2.0/3.0) * inv4pi);
            Gauss_sphere_mean2 = Gauss_sphere_mean * Gauss_sphere_mean
            Psum = 0.0
            for j in range(1,PointsInr):
                qr = q[i] * r[j]
                sinc = math.sin(qr)/qr
                Psum = Psum + (gr[j] * psi_c2 - (hr[j] + jr[j]) * psi_c * Gauss_sphere_mean + kr[j] * Gauss_sphere_mean2 ) * sinc
            Pq[i] = Psum
        Pq = Pq / Pq[0] # Normalize P(q) such that P(0) = 1

        # plot P(q)
        if capp_version == "capp_windows.exe":
            print("Plotting not available on Windows") # then it crashes
        else:
            self.Figure2 = pylab.figure(2)
            Subplot2 = pylab.subplot(111)

            length_of_name = len(pr_filename)
            if self.Water_layer_button.GetValue():
                if self.Exclude_water_layer_M_button.GetValue() or self.Exclude_water_layer_OPM_button.GetValue():
                    structure_name = "'..." + pr_filename[length_of_name-10-9:-9] + "'" + " + WL " + self.Water_layer_contrast_box.GetValue()  + "% (Excl. TMD)"
                else:
                    structure_name = "'..." + pr_filename[length_of_name-10-9:-9] + "'" + " + WL " + self.Water_layer_contrast_box.GetValue()  + "%"
            else:
                structure_name = "'..." + pr_filename[length_of_name-10-7:-7] + "'"
            if self.SANS_button.GetValue():
                structure_name = structure_name + ", SANS " + self.SANS_solvent_box.GetValue() + "% D2O"
            else:
                structure_name = structure_name + ", SAXS "
            Subplot2.plot(q, Pq, label=structure_name)
            Subplot2.set_xscale('log', nonposx = 'clip')
            Subplot2.set_yscale('log', nonposy = 'clip')
            Subplot2.legend(loc=3,fontsize=9)
            Subplot2.set_xlabel('q [$1/\AA$]',fontsize=14)
            Subplot2.set_ylabel('P(q)',fontsize=14)
            pylab.suptitle('Form Factor for the structure(s)',fontsize=14)
            pylab.show()

        # export P(q)
        Pq_filename = pr_filename[:-6] + "Pq.dat"
        Pq_fid = open(Pq_filename, "w")
        Pq_fid.write("q [1/AA]  P(q)\n")
        for i in range(0,PointsInQ):
            Pq_fid.write("%f %f \n" % (q[i],Pq[i]))
        Pq_fid.close()

        # disable P(q) calculation button
        self.CalcPqButton.Disable()

        # enable  Fit button
        self.FitButton.Enable()

### Define function for Fitting Bg and scale
    def FitFnc(self, event):
        
        # check for data file
        if self.DataPathStr == 'Non':
            Message = 'Please provide a data file (q,I,dI), recalculate P(q) (to obtain the right q-values) and try again.'
            wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)
            self.CalcPqButton.Enable()
            self.FitButton.Disable()
        else:
            filename, extension = os.path.splitext(self.PDBPathStr)
            if self.Water_layer_button.GetValue():
                pr_filename = filename + "_w_pr.dat"
            else:
                pr_filename = filename + "_pr.dat"
            Pq_filename = pr_filename[:-6] + "Pq.dat"
            file = open(self.DataPathStr,'r')
            #count number of header lines
            count = 0
            STOP = 0
            while STOP < 1:
                line = file.readline()
                #print("count:", count)
                #print("line:", line)
                numbers_str = line.split()
                #print("numbers_str", numbers_str)
                try:
                    numbers_float = [float(x) for x in numbers_str]
                    STOP = 1
                except:
                    count = count + 1
                    #print("count:",count)
            file.close()
            q,I,dI = np.genfromtxt(self.DataPathStr, skiprows=count,usecols=[0,1,2], unpack=True) #import data
            if self.nm_button.GetValue():
                q = q*0.1 # convert from nm to angstrom
            Pq = np.genfromtxt(Pq_filename, skiprows=1,usecols=[1], unpack=True) #import P(q)
            S0 = I[2] #Scale, initial guess
            B0 = 0.001 #Bg, initial guess
            #print('S0 is ', S0)
            def func(q,S,B):
                return S * Pq + B
            popt, pcov = curve_fit(func,q,I,sigma=dI,p0=[S0,B0])
            #print('popt is ', popt)
            # plot data and fit
            if capp_version == "capp_windows.exe":
                Message = 'Plotting not (yet) available for windows (program crashes).'
                wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)
                self.FitButton.Disable()
                self.CalcPqButton.Enable()
            else:
                self.Figure3 = pylab.figure(3)
                Subplot3 = pylab.subplot(111)
                Subplot3.plot(q, I, label="data")
                Subplot3 = pylab.subplot(111)
                chi2r = np.sum(((I - func(q,*popt) )/dI)**2)/len(q)
                #print("chi2r = ", chi2r)
                Subplot3.plot(q, func(q,*popt), label='fit: Scale=%5.2f, Bg=%5.3f, Chi2r=%4.2f' % (popt[0],popt[1],chi2r))
                Subplot3.set_xscale('log', nonposx = 'clip')
                Subplot3.set_yscale('log', nonposy = 'clip')
                Subplot3.legend(loc=3,fontsize=9)
                Subplot3.set_xlabel('q [$1/\AA$]',fontsize=14)
                Subplot3.set_ylabel('I(q)',fontsize=14)
                pylab.suptitle('Fit structure to data',fontsize=14)
                pylab.show()
                if chi2r > 100:
                    Message = 'Very poor fit! (Chi2r>100) - is the data right? units of data in aa or nm? is the structure right? oligomerization? aggregation?'
                    wx.MessageBox(Message, "CaPP", wx.OK | wx.ICON_INFORMATION)
                    self.CalcPqButton.Enable()
            
            # disable Fit button
            self.FitButton.Disable()

    ### Define function for filebrowsing for PDB
    def BrowsePDBFnc(self, event):
        
        FileDialogWindow = wx.FileDialog(None, 'Please select PDB-file...', os.getcwd(), defaultFile = '',style = wx.OPEN)

        if FileDialogWindow.ShowModal() == wx.ID_OK:
            self.PDBPathStr = FileDialogWindow.GetPath()
            PDBPathDisplayStr = str(self.PDBPathStr)

            while len(PDBPathDisplayStr) > 49:
                PDBPathDisplayStr = PDBPathDisplayStr[1:]

            if len(self.PDBPathStr) > 49:
                PDBPathDisplayStr = '...' + PDBPathDisplayStr

            self.PDBPathTxt.SetLabel(PDBPathDisplayStr)

            DirectoryStr = os.path.dirname(self.PDBPathStr)
            os.chdir(DirectoryStr)

            self.CalculateButton.Enable()
            self.CalculateRgButton.Enable()
        
        FileDialogWindow.Destroy()

### Define function for filebrowsing for data
    def BrowseDataFnc(self, event):
        FileDialogWindow = wx.FileDialog(None, 'Please select Data-file...', os.getcwd(), defaultFile = '',style = wx.OPEN)
        
        if FileDialogWindow.ShowModal() == wx.ID_OK:
            self.DataPathStr = FileDialogWindow.GetPath()
            DataPathDisplayStr = str(self.DataPathStr)
            
            while len(DataPathDisplayStr) > 49:
                DataPathDisplayStr = DataPathDisplayStr[1:]
            
            if len(self.DataPathStr) > 49:
                DataPathDisplayStr = '...' + DataPathDisplayStr
        
            self.DataPathTxt.SetLabel(DataPathDisplayStr)
            
            DirectoryStr = os.path.dirname(self.DataPathStr)
            os.chdir(DirectoryStr)
            
            self.CalculateButton.Enable()
            self.CalculateRgButton.Enable()
        
        FileDialogWindow.Destroy()
    
    ### Define exit function
    def CloseWindowFnc(self, event):
        try:
            pylab.close()
        except:
            pass
        sys.exit(0)

### Boilerplate code linking program and widgets
if __name__ == '__main__':
    App = wx.App(False)
    frame = MainCls(parent = None, id = -1)
    frame.Show()
    App.MainLoop()
