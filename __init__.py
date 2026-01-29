# This Python 3.x file uses the following encoding: utf-8
# Crossection plugin for PyMol 2.x (Windows version) Copyright Notice.
# ====================================================================
#
# The crossection plugin source code is copyrighted, but you can freely
# use and copy it as long as you don't change or remove any of the
# copyright notices.
#
# ----------------------------------------------------------------------
# Crossection plugin is Copyright (C) 2020 by Goetz Parsiegla
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Goetz Parsiegla not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# GOETZ PARSIEGLA DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------
#
# It calculates the consecutive cross-sections of a protein model and stores them in "model_profile.txt"
# The Plugin interface architecture is an adaptation from Daniel Seeligers 'Vina Autodock'-Plugin
# (J. Comput.-Aided Mol. Des. 24:417-422 (2010)). The idea for the creation of cross-sections in PyMol
# is adapted from the script to create an image with a slice region (http://www.pymolwiki.org/index.php/Gallery)
# The plugin creates images from a reference with a given surface in green and compares them to images of the
# slice surfaces at rising depths created in red. The reference follows the slice depth to correct the perspective
# decrease of the surface size with rising depth. 
# If you find bugs or have any suggestions for future versions contact me : goetz.parsiegla@imm.cnrs.fr

import sys
import os
from pymol import cmd
from pymol.cgo import *
from PIL import Image
from time import sleep
# pymol.Qt is a wrapper which provides the PySide2/Qt5/Qt4 interface
# if it has been installed in python before !
from pymol.Qt import QtWidgets, QtCore

__version__ = "1.4"

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Crossection', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None

def run_plugin_gui():
    global dialog
    if dialog is None:
        dialog = QtWidgets.QDialog() # now global dialog holds an object

    # filename of our UI file
    uifile = os.path.join(os.path.dirname(__file__), 'form.ui')

    # load the UI file into our dialog
    from pymol.Qt.utils import loadUi
    form = loadUi(uifile, dialog)
    Crossection(form) # call the plugin class and pass the form as an argument
    dialog.show()

#---------------------------------------------------------------------------------------
# defaults and globals
depth = 0 	    # start at pentration depth 0              
imwidth = 600
imheigth = 500
prot = ""
wd = ""
cmd.do("set solvent_radius, 1.6")
#---------------------------------------------------------------------------------------

#==========================================================================
#    THE MAJOR PLUGIN CLASS

class Crossection:
    def __init__(self, form):
        
        # Make Buttongroups
        Buttongroup_1 = QtWidgets.QButtonGroup()
        Buttongroup_1.addButton(form.radioButton, 0)
        Buttongroup_1.addButton(form.radioButton_2, 1)

        # Assign nonlocal variables to lineEdits 
        status_line = form.lineEdit

        # Defaults
        form.doubleSpinBox.setMinimum(0)
        form.doubleSpinBox.setMaximum(100)
        form.doubleSpinBox_2.setMinimum(0)
        form.doubleSpinBox_2.setMaximum(100)     
        form.doubleSpinBox.setValue(2) #  profile resolution in A
        form.doubleSpinBox_2.setValue(10) # profile depth in A

        intro_text = """
        <html><body>Calculates the consecutive cross-section areas of a protein model in a choosen slicing distance \
        up to a pentration depth and stores them in file profile_model.csv.<br>
        1. Load the model in PyMol<br>
        2. Import and select it in the plugin<br>
        3. Build its surface and orient the model<br>
        4. Execute 'Make Profile'</body></html>
        """
        status_text ="Ready..... (version: %s)" % __version__

        # Page buildup 
        form.textBrowser.setHtml(intro_text)
        status_line.setText(status_text)

        def set_status_line(text):
            status_line.clear()
            status_line.insert(text)

        def object_selection_changed():
            prot = form.comboBox.currentText()
            txt = 'Selected object '+prot
            set_status_line(txt)

        def image_mode_changed():
            if form.checkBox_2.isChecked():
                action = 'Enabled'
            else:
                action = 'Disabled'
            txt = 'Keep raw images'+' <'+action+'>'
            set_status_line(txt)

        def surface_mode_changed():
            if form.checkBox.isChecked():
                action = 'Enabled'
            else:
                action = 'Disabled'
            txt = 'All atoms for surface'+' <'+action+'>'
            set_status_line(txt)

    #------------------------------------------------------------------
           
        def import_objects():
            form.comboBox.clear()
            lst = cmd.get_names()
            if 'sele' in lst:
                lst.remove('sele')
            if 'cgo' in lst:
                lst.remove('cgo')
            form.comboBox.addItems(lst)
     
        def create_surface(prot):
            cmd.color("white","name C*")
            cmd.color("red","name O*")
            cmd.color("blue","name N*")
            cmd.color("yellow","name S*")
            if form.checkBox.isChecked():
                cmd.do("set surface_mode, 1")
            else:
                cmd.do("set surface_mode, 0")
            cmd.clip("atoms", 90, prot)
            cmd.show(representation="surface", selection=prot)
            cmd.show(representation="lines", selection=prot)
            
        def add_hydrogens(prot):
            cmd.h_add(selection=prot)

        def surface_object_selected():
            prot = form.comboBox.currentText()
            if prot == "":
                set_status_line('No object selected ! ')
            else:
                create_surface(prot)

        def hydrogens_object_selected():
            prot = form.comboBox.currentText()
            if prot == "":
                set_status_line('No object selected ! ')
            else:
                add_hydrogens(prot)

        def profile_object_selected():
            global prot
            prot = form.comboBox.currentText()
            if prot == "":
                set_status_line('No object selected ! ')
            else:
                set_status_line('Calculating crossections for '+str(prot)+' ... ')
                do_crossection(prot)          

        def do_crossection(prot):
            setup_slice_parameters()
            cmd.hide(representation="everything", selection="all")
            cmd.show(representation="surface", selection=prot)
            make_output_file()
            do_slicing()
            close_output_file()
            set_status_line('Done ... Profile written to : profile_'+str(prot)+'.csv')
            cleanup()
            
        def setup_slice_parameters():
            # image parameters
            cmd.viewport(imwidth, imheigth)
            cmd.bg_color(color="black")
            cmd.do("set ray_interior_color, red")
            cmd.do("set opaque_background") 
            # surface parameters
            cmd.color("white","name C*")
            cmd.color("white","name O*")
            cmd.color("white","name N*")
            cmd.color("white","name S*")
            cmd.color("white","name H*")        
            if form.checkBox.isChecked():
                cmd.do("set surface_mode, 1")
            else:
                cmd.do("set surface_mode, 0")
            # enable pure colors withouot lighteffects
            cmd.do("set ambient, 1") 
            cmd.do("set reflect, 0")
            cmd.do("set specular, off")
            # must disable depth cue and shadows
            cmd.do("unset depth_cue")
            cmd.do("unset ray_shadows")
            cmd.do("set ray_trace_mode, 0")
            cmd.do("set backface_cull, 0")
            # suppress internal cavities in crossection
            cmd.do("set cavity_cull, 30")


        def make_sphere(nx, ny, nz, nr):	
    	   # create reference sphere with cross-section 10 A2	
        	obj = [
            COLOR, 0.0, 1.0, 0.0,
            SPHERE, nx,ny,nz,nr
            ]
        	cmd.load_cgo(obj, "sphere") 	

        def do_slicing():        
            # move clipping planes on protein surface envelop
            cmd.center(prot, origin = 0)
            cmd.zoom(prot, buffer = 2.5, complete = 1)
            cmd.clip("atoms", 1.2, prot)
            view = cmd.get_view()
            clip_slab = view[16]-view[15]
            step_width = form.doubleSpinBox.value()
            nr = 1.784 # = sqr(10/pi) for area=10
            nz = (-view[11]-view[15]-(2*nr)-1.2)
            make_sphere(view[12],view[13],view[14],nr)
            cmd.translate([0,0,nz], object="sphere", camera = 1)
            print (" initial sphere position from origin(nz) = "+str(nz))
            cmd.set_view ((\
            view[0], view[1], view[2],\
            view[3], view[4], view[5],\
            view[6], view[7], view[8],\
            view[9], view[10], view[11],\
            view[12], view[13], view[14],\
            view[15], view[16], view[17] ))    
            if Buttongroup_1.checkedId() == 1:
                profile_depth = form.doubleSpinBox_2.value()
                slices= int(profile_depth/step_width)
            else:
                slices = int(abs(clip_slab/step_width)+1)    
            global depth
            print ("starting depth : "+str(depth))
            print ("step width : "+str(step_width))
            print ("number of slices : "+str(slices))
            print ("envelop depth : "+str(clip_slab))

            p=0       
            pixel_count(p)
            for p in range (slices):
                p+=1
                cmd.translate([0,0,-step_width], object="sphere", camera = 1)
                cmd.clip("move", -step_width)
                depth=depth+step_width
                pixel_count(p)
                            
        def pixel_count(p): 	
    	   #determine number of pixels in reference
            cmd.hide(representation="surface", selection=prot)
            cmd.show(representation="cgo", selection="sphere") 
            wd = os.getcwd()
            filename = ('refslice'+str(p)+'.png')
            print (str(filename))
            cmd.ray(imwidth, imheigth)
            cmd.png(str(filename))
            sleep(0.2) # wait for picture to be created
            im = Image.open(filename)
            pix = im.load()
    	
            countgreen = 0
            countall = 0
            for i in range (imwidth):
                for j in range (imheigth):
                    if pix[i,j][1]==255 and pix[i,j][0]==0 :
                        countgreen+=1
                        countall+=1
            os.remove(filename)
            print ("reference area refslice"+str(p)+" done")	

    	   #determine number of pixels in slice
    	
            cmd.hide(representation="cgo", selection="sphere")
            cmd.show(representation="surface", selection=prot)

            print ("wd = "+str(wd))
            filename = (str(prot)+'_slice_'+str(p)+'.png')
            print (str(filename))
            cmd.ray (imwidth, imheigth)
            cmd.png (str(filename))
            sleep(1.0) # wait for picture to be created	
            im = Image.open(filename)
            pix = im.load()
    	
            countred = 0
            countall = 0
            for i in range (imwidth):
                for j in range (imheigth):
                    if pix[i,j][0]==255 and pix[i,j][1]==0 :
                        countred+=1
                    countall+=1
            surface = countred / (countgreen/10)		
            print ("number of red pixels in slice : " + str(countred))
            print ("number of green pixels in reference : " + str(countgreen))
            print (" total number of  pixels in slice : " + str(countall))
            print (" surface [A2] : " + str(surface))
            print (" at penetration depth [A] : " + str(depth))
            if form.checkBox_2.isChecked():
                pass
            else:
                os.remove(filename) 
            
            #write result to profile.txt	

            outfile.write(";"+str(depth)+";"+str(surface)+"\n")
            print ("slice"+str(p)+" done")	

        def make_output_file():
            # open profile output file
            global wd
            wd=os.getcwd()
            global outfile
            outfile=open(str(wd)+"\\profile_"+str(prot)+".csv", "w")
            print ("Opened File for output : "+str(wd)+"\\profile_"+str(prot)+".csv")
        
            title1 = str("Profile for model : "+prot+"\n"+"\n")
            title2 = str(";penetration depth [A]"+";"+"cross-section [A2]"+"\n")
            outfile.write(title1)
            outfile.write(title2)
            
        def close_output_file():    
            outfile.close()		
            print ("File : profile_"+str(prot)+".csv written to " + str(wd))
        
        def cleanup():
            # enable normal colors with lighteffects
            cmd.do("set ambient, 0.1") 
            cmd.do("set reflect, 0.5")
            cmd.do("set ray_shadows")
            cmd.do("set specular, on")
            cmd.color("white","name C*")
            cmd.color("red","name O*")
            cmd.color("blue","name N*")
            cmd.color("yellow","name S*")
            # reset globals
            global depth
            depth = 0
            #delete last sphere
            cmd.delete("sphere")
            cmd.refresh()

        # Assign buttons
        form.pushButton.clicked.connect(profile_object_selected)
        form.pushButton_2.clicked.connect(hydrogens_object_selected)
        form.pushButton_3.clicked.connect(import_objects)
        form.pushButton_4.clicked.connect(surface_object_selected)
                 
        # Callback bindings
        form.checkBox.stateChanged.connect(surface_mode_changed) 
        form.checkBox_2.stateChanged.connect(image_mode_changed)
        form.comboBox.activated.connect(object_selection_changed)