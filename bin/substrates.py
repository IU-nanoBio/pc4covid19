# substrates  Out:Plots

import os, math
from pathlib import Path
from ipywidgets import Layout, Label, Text, Checkbox, Button, BoundedIntText, HBox, VBox, Box, \
    FloatText, Dropdown, SelectMultiple, interactive
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle, Ellipse, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.colors as mplc
from matplotlib import gridspec
from collections import deque
from pyMCDS import pyMCDS
import numpy as np
import scipy.io
import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html
import glob
import platform
import zipfile
from debug import debug_view 
import warnings

hublib_flag = True
if platform.system() != 'Windows':
    try:
#        print("Trying to import hublib.ui")
        from hublib.ui import Download
    except:
        hublib_flag = False
else:
    hublib_flag = False

#warnings.warn(message, mplDeprecation, stacklevel=1)
warnings.filterwarnings("ignore")

class SubstrateTab(object):

    def __init__(self):
        
        self.output_dir = '.'
        # self.output_dir = 'tmpdir'

        # These are recomputed below 
        # self.figsize_width_substrate = 15.0  # allow extra for colormap
        # self.figsize_height_substrate = 12.5

        foo=12.0
        self.width_substrate = foo  # allow extra for colormap
        self.height_substrate = foo

        self.figsize_width_svg = foo
        self.figsize_height_svg = foo
        self.width_svg = foo
        self.height_svg = foo

        self.ax0 = None
        self.ax1 = None

        self.custom_data_plotted = False

        # self.fig = plt.figure(figsize=(7.2,6))  # this strange figsize results in a ~square contour plot

        self.first_time = True
        self.modulo = 1

        self.use_defaults = True

        self.svg_delta_t = 1
        self.substrate_delta_t = 1
        self.svg_frame = 1
        self.substrate_frame = 1

        self.customized_output_freq = False
        self.therapy_activation_time = 1000000
        self.max_svg_frame_pre_therapy = 1000000
        self.max_substrate_frame_pre_therapy = 1000000

        self.svg_xmin = 0

        # Probably don't want to hardwire these if we allow changing the domain size
        # self.svg_xrange = 2000
        # self.xmin = -1000.
        # self.xmax = 1000.
        # self.ymin = -1000.
        # self.ymax = 1000.
        # self.x_range = 2000.
        # self.y_range = 2000.

        self.show_nucleus = False
        self.show_edge = True

        substrates_default_disabled_flag = True  # True = disable them by default; False=enable them

        # initial value
        self.field_index = 4
        # self.field_index = self.mcds_field.value + 4

        self.skip_cb = False

        # define dummy size of mesh (set in the tool's primary module)
        self.numx = 0
        self.numy = 0

        self.title_str = ''

        tab_height = '600px'
        tab_height = '500px'
        constWidth = '180px'
        constWidth2 = '150px'
        tab_layout = Layout(width='900px',   # border='2px solid black',
                            height=tab_height, ) #overflow_y='scroll')

        max_frames = 1   
        # self.mcds_plot = interactive(self.plot_substrate, frame=(0, max_frames), continuous_update=False)  
        # self.i_plot = interactive(self.plot_plots, frame=(0, max_frames), continuous_update=False)  
        self.i_plot = interactive(self.plot_substrate, frame=(0, max_frames), continuous_update=False)  

        # "plot_size" controls the size of the tab height, not the plot (rf. figsize for that)
        # NOTE: the Substrates Plot tab has an extra row of widgets at the top of it (cf. Cell Plots tab)
        svg_plot_size = '700px'
        svg_plot_size = '900px'
        plot_area_width = '1500px'
        plot_area_height = '900px'
        self.i_plot.layout.width = plot_area_width 
        self.i_plot.layout.height = plot_area_height 

        self.fontsize = 20
        # self.fontsize = 30


        #============  new GUI =================
        self.max_frames = BoundedIntText(value=0,description='# cell frames',layout=Layout(width='160px'))  # border='1px solid black',
        self.cells_toggle = Checkbox(description='Cells',value=True, style = {'description_width': 'initial'}, layout=Layout(width='110px', )) #border='1px solid black'))
        self.cell_edges_toggle = Checkbox(description='edge',value=self.show_edge, style = {'description_width': 'initial'}, layout=Layout(width='110px',))  #   align_items='stretch',
            

        layout1 = Layout(display='flex',
                            flex_flow='row',
                            align_items='center',
                            width='25%', )  #border='1px solid black')
        hb1=HBox([self.cells_toggle,self.cell_edges_toggle ])  # layout=layout1)
        cells_vbox=VBox([self.max_frames, hb1], layout=Layout(width='350px',border='1px solid black',))
        #--------------------------
        self.substrates_toggle=Checkbox(description='Substrates', style = {'description_width': 'initial'})

        self.field_min_max = {'assembled_virion':[0.,1.,False]  }
        # hacky I know, but make a dict that's got (key,value) reversed from the dict in the Dropdown below

        # ipywidgets 8 docs: Selection widgets no longer accept a dictionary of options. Pass a list of key-value pairs instead.
        self.field_dict = {0:'director signal', 1:'cargo signal'}

        self.substrate_choice = Dropdown(options={'assembled_virion': 0},layout=Layout(width='150px'))
        self.colormap_dd = Dropdown(options=['viridis', 'jet', 'YlOrRd'],value='YlOrRd',layout=Layout(width='200px'))
        hb2=HBox([self.substrates_toggle,self.substrate_choice,self.colormap_dd], layout=Layout(width='350px', )) # border='1px solid black',))

        self.colormap_fixed_toggle = Checkbox(description='Fix',style = {'description_width': 'initial'}, layout=Layout(width='60px'))
        constWidth2 = '160px'
        self.colormap_min = FloatText(
                    description='Min',
                    value=0,
                    step = 0.1, 
                    layout=Layout(width=constWidth2),)
        self.colormap_max = FloatText(
                    description='Max',
                    value=38,
                    step = 0.1,
                    layout=Layout(width=constWidth2),)
        # hb3=HBox([colormap_fixed_toggle,colormap_min,colormap_max], layout=Layout(width='500px',justify_content='flex-start'))
        hb3=HBox([self.colormap_fixed_toggle,self.colormap_min,self.colormap_max], layout=Layout(justify_content='flex-start'))  # border='1px solid black',

        substrate_vbox=VBox([hb2, hb3], layout=Layout(width='500px'))
        #--------------------------
        self.custom_data_toggle = Checkbox(
                    description='Custom plot',
                    disabled=False,
                    value=False,
            style = {'description_width': 'initial'},
            layout=Layout(width='130px', )  #  border='1px solid black',)
        #           layout=Layout(width=constWidth2),
                )
        self.custom_data_update_button= Button(description='Update', layout=Layout(width='120px', ))  #,style = {'description_width': 'initial'})
        custom_data_vbox1 = VBox([self.custom_data_toggle, self.custom_data_update_button], layout=Layout(justify_content='flex-start',border='1px solid black',))  # width='330px',

        self.custom_data_choice = SelectMultiple(options=['assembled_virion','susceptible','infected', 'dead'],
                    value=['assembled_virion'],rows=3,  layout=Layout(width='160px', ))
        custom_data_hbox = HBox([custom_data_vbox1, self.custom_data_choice])

        #gui=HBox([cells_vbox, substrate_vbox, custom_data_hbox], justify_content='center')  # vs. 'flex-start   , layout=Layout(width='900px'))

        #==========================================================================

        # ------- "observe" functionality (callbacks)
        self.max_frames.observe(self.update_max_frames)

        # self.field_min_max = {'dummy': [0., 1., False]}
        # NOTE: manually setting these for now (vs. parsing them out of data/initial.xml)

        self.mcds_field = Dropdown(
            options={'director signal': 0, 'cargo signal':1},
            value=0,
            disabled = substrates_default_disabled_flag,
            #     description='Field',
           layout=Layout(width=constWidth)
        )
        # print("substrate __init__: self.mcds_field.value=",self.mcds_field.value)
#        self.mcds_field.observe(self.mcds_field_cb)
        self.mcds_field.observe(self.mcds_field_changed_cb)

        self.field_cmap = Dropdown(
            options=['viridis', 'jet', 'YlOrRd'],
            value='YlOrRd',
            #     description='Field',
           layout=Layout(width=constWidth)
        )
#        self.field_cmap.observe(self.plot_substrate)
        self.field_cmap.observe(self.mcds_field_cb)

        self.cmap_fixed_toggle = Checkbox(
            description='Fix',
            disabled=False,
#           layout=Layout(width=constWidth2),
        )
        self.cmap_fixed_toggle.observe(self.mcds_field_cb)

        self.cmap_min = FloatText(
            description='Min',
            value=0,
            step = 0.1,
            disabled=True,
            layout=Layout(width=constWidth2),
        )
        self.cmap_min.observe(self.mcds_field_cb)

        self.cmap_max = FloatText(
            description='Max',
            value=38,
            step = 0.1,
            disabled=True,
            layout=Layout(width=constWidth2),
        )
        self.cmap_max.observe(self.mcds_field_cb)

        def cmap_fixed_toggle_cb(b):
            field_name = self.field_dict[self.mcds_field.value]
            # print(self.cmap_fixed_toggle.value)
            if (self.cmap_fixed_toggle.value):  # toggle on fixed range
                self.cmap_min.disabled = False
                self.cmap_max.disabled = False
                self.field_min_max[field_name][0] = self.cmap_min.value
                self.field_min_max[field_name][1] = self.cmap_max.value
                self.field_min_max[field_name][2] = True
                # self.save_min_max.disabled = False
            else:  # toggle off fixed range
                self.cmap_min.disabled = True
                self.cmap_max.disabled = True
                self.field_min_max[field_name][2] = False
                # self.save_min_max.disabled = True
#            self.mcds_field_cb()
            self.i_plot.update()

        self.cmap_fixed_toggle.observe(cmap_fixed_toggle_cb)

        def cell_edges_toggle_cb(b):
            # self.update()
            if (self.cell_edges_toggle.value):  
                self.show_edge = True
            else:
                self.show_edge = False
            self.i_plot.update()

        self.cell_edges_toggle.observe(cell_edges_toggle_cb)

        def cells_toggle_cb(b):
            # self.update()
            self.i_plot.update()
            if (self.cells_toggle.value):
                self.cell_edges_toggle.disabled = False
                # self.cell_nucleus_toggle.disabled = False
            else:
                self.cell_edges_toggle.disabled = True
                # self.cell_nucleus_toggle.disabled = True

        self.cells_toggle.observe(cells_toggle_cb)

        def substrates_toggle_cb(b):
            if (self.substrates_toggle.value):  # seems bass-ackwards
                self.cmap_fixed_toggle.disabled = False
                self.cmap_min.disabled = False
                self.cmap_max.disabled = False
                self.mcds_field.disabled = False
                self.field_cmap.disabled = False
            else:
                self.cmap_fixed_toggle.disabled = True
                self.cmap_min.disabled = True
                self.cmap_max.disabled = True
                self.mcds_field.disabled = True
                self.field_cmap.disabled = True

        self.substrates_toggle.observe(substrates_toggle_cb)

        #---------------------
        def custom_data_toggle_cb(b):
            # self.update()
            if (self.custom_data_toggle.value):  # seems bass-ackwards
                self.custom_data_choice.disabled = False
            else:
                self.custom_data_choice.disabled = True
            # self.i_plot.update()

        self.custom_data_toggle.observe(custom_data_toggle_cb)

        self.custom_data_update_button.on_click(self.update_custom_data)


        #---------------------
        help_label = Label('select slider: drag or left/right arrows')

        controls_box = HBox([cells_vbox, substrate_vbox, custom_data_hbox], justify_content='center')  # vs. 'flex-start   , layout=Layout(width='900px'))

        if (hublib_flag):
            self.download_button = Download('mcds.zip', style='warning', icon='cloud-download', 
                                                tooltip='Download data', cb=self.download_cb)

            self.download_svg_button = Download('svg.zip', style='warning', icon='cloud-download', 
                                            tooltip='You need to allow pop-ups in your browser', cb=self.download_svg_cb)
            download_row = HBox([self.download_button.w, self.download_svg_button.w, Label("Download all cell plots (browser must allow pop-ups).")])

            # box_layout = Layout(border='0px solid')
            # controls_box = VBox([row1, row2])  # ,width='50%', layout=box_layout)
            # controls_box = HBox([cells_vbox, substrate_vbox, custom_data_hbox], justify_content='center')  # vs. 'flex-start   , layout=Layout(width='900px'))
            self.tab = VBox([controls_box, self.i_plot, download_row])
        else:
            # self.tab = VBox([row1, row2])
            # self.tab = VBox([row1, row2, self.i_plot])
            self.tab = VBox([controls_box, self.i_plot])

    #---------------------------------------------------
    def update_dropdown_fields(self, data_dir):
        # print('update_dropdown_fields called --------')
        self.output_dir = data_dir
        tree = None
        try:
            fname = os.path.join(self.output_dir, "initial.xml")
            tree = ET.parse(fname)
            xml_root = tree.getroot()
        except:
            print("Cannot open ",fname," to read info, e.g., names of substrate fields.")
            return

        xml_root = tree.getroot()
        self.field_min_max = {}
        self.field_dict = {}
        dropdown_options = {}
        uep = xml_root.find('.//variables')
        comment_str = ""
        field_idx = 0
        if (uep):
            for elm in uep.findall('variable'):
                # print("-----> ",elm.attrib['name'])
                field_name = elm.attrib['name']
                self.field_min_max[field_name] = [0., 1., False]
                self.field_dict[field_idx] = field_name
                dropdown_options[field_name] = field_idx

                self.field_min_max[field_name][0] = 0   
                self.field_min_max[field_name][1] = 1

                # self.field_min_max[field_name][0] = field_idx   #rwh: helps debug
                # self.field_min_max[field_name][1] = field_idx+1   
                self.field_min_max[field_name][2] = False
                field_idx += 1

#        constWidth = '180px'
        # print('options=',dropdown_options)
        # print(self.field_min_max)  # debug
        self.mcds_field.value = 0
        self.mcds_field.options = dropdown_options
#         self.mcds_field = Dropdown(
# #            options={'oxygen': 0, 'glucose': 1},
#             options=dropdown_options,
#             value=0,
#             #     description='Field',
#            layout=Layout(width=constWidth)
#         )

    # def update_max_frames_expected(self, value):  # called when beginning an interactive Run
    #     self.max_frames.value = value  # assumes naming scheme: "snapshot%08d.svg"
    #     self.mcds_plot.children[0].max = self.max_frames.value

#------------------------------------------------------------------------------
    def update_params(self, config_tab, user_params_tab):
        # xml_root.find(".//x_min").text = str(self.xmin.value)
        # xml_root.find(".//x_max").text = str(self.xmax.value)
        # xml_root.find(".//dx").text = str(self.xdelta.value)
        # xml_root.find(".//y_min").text = str(self.ymin.value)
        # xml_root.find(".//y_max").text = str(self.ymax.value)
        # xml_root.find(".//dy").text = str(self.ydelta.value)
        # xml_root.find(".//z_min").text = str(self.zmin.value)
        # xml_root.find(".//z_max").text = str(self.zmax.value)
        # xml_root.find(".//dz").text = str(self.zdelta.value)

        self.xmin = config_tab.xmin.value 
        self.xmax = config_tab.xmax.value 
        self.x_range = self.xmax - self.xmin
        self.svg_xrange = self.xmax - self.xmin
        self.ymin = config_tab.ymin.value
        self.ymax = config_tab.ymax.value 
        self.y_range = self.ymax - self.ymin

        self.numx =  math.ceil( (self.xmax - self.xmin) / config_tab.xdelta.value)
        self.numy =  math.ceil( (self.ymax - self.ymin) / config_tab.ydelta.value)

        if (self.x_range > self.y_range):  
            ratio = self.y_range / self.x_range
            self.figsize_width_substrate = self.width_substrate  # allow extra for colormap
            self.figsize_height_substrate = self.height_substrate * ratio
            self.figsize_width_svg = self.width_svg
            self.figsize_height_svg = self.height_svg * ratio
        else:   # x < y
            ratio = self.x_range / self.y_range
            self.figsize_width_substrate = self.width_substrate * ratio 
            self.figsize_height_substrate = self.height_substrate
            self.figsize_width_svg = self.width_svg * ratio
            self.figsize_height_svg = self.height_svg
        print('update_params():  self.figsize_width_svg= ',self.figsize_width_svg)

        self.svg_flag = config_tab.toggle_svg.value
        self.substrates_flag = config_tab.toggle_mcds.value
        # print("substrates: update_params(): svg_flag, toggle=",self.svg_flag,config_tab.toggle_svg.value)        
        # print("substrates: update_params(): self.substrates_flag = ",self.substrates_flag)
        self.svg_delta_t = config_tab.svg_interval.value
        self.substrate_delta_t = config_tab.mcds_interval.value
        self.modulo = int(self.substrate_delta_t / self.svg_delta_t)
        # print("substrates: update_params(): modulo=",self.modulo)        

        if self.customized_output_freq:
#            self.therapy_activation_time = user_params_tab.therapy_activation_time.value   # NOTE: edit for user param name
            # print("substrates: update_params(): therapy_activation_time=",self.therapy_activation_time)
            self.max_svg_frame_pre_therapy = int(self.therapy_activation_time/self.svg_delta_t)
            self.max_substrate_frame_pre_therapy = int(self.therapy_activation_time/self.substrate_delta_t)

#------------------------------------------------------------------------------
#    def update(self, rdir):
#   Called from driver module (e.g., pc4*.py) (among other places?)
    def update(self, rdir=''):
        # with debug_view:
        #     print("substrates: update rdir=", rdir)        
        # print("substrates: update rdir=", rdir)        

        if rdir:
            self.output_dir = rdir

        # print('update(): self.output_dir = ', self.output_dir)

        if self.first_time:
        # if True:
            self.first_time = False
            full_xml_filename = Path(os.path.join(self.output_dir, 'config.xml'))
            # print("substrates: update(), config.xml = ",full_xml_filename)        
            # self.num_svgs = len(glob.glob(os.path.join(self.output_dir, 'snap*.svg')))
            # self.num_substrates = len(glob.glob(os.path.join(self.output_dir, 'output*.xml')))
            # print("substrates: num_svgs,num_substrates =",self.num_svgs,self.num_substrates)        
            # argh - no! If no files created, then denom = -1
            # self.modulo = int((self.num_svgs - 1) / (self.num_substrates - 1))
            # print("substrates: update(): modulo=",self.modulo)        
            if full_xml_filename.is_file():
                tree = ET.parse(full_xml_filename)  # this file cannot be overwritten; part of tool distro
                xml_root = tree.getroot()
                self.svg_delta_t = int(xml_root.find(".//SVG//interval").text)
                self.substrate_delta_t = int(xml_root.find(".//full_data//interval").text)
                # print("substrates: svg,substrate delta_t values=",self.svg_delta_t,self.substrate_delta_t)        
                self.modulo = int(self.substrate_delta_t / self.svg_delta_t)
                # print("substrates: update(): modulo=",self.modulo)        


        # all_files = sorted(glob.glob(os.path.join(self.output_dir, 'output*.xml')))  # if the substrates/MCDS

        all_files = sorted(glob.glob(os.path.join(self.output_dir, 'snap*.svg')))   # if .svg
        if len(all_files) > 0:
            last_file = all_files[-1]
            self.max_frames.value = int(last_file[-12:-4])  # assumes naming scheme: "snapshot%08d.svg"
        else:
            substrate_files = sorted(glob.glob(os.path.join(self.output_dir, 'output*.xml')))
            if len(substrate_files) > 0:
                last_file = substrate_files[-1]
                self.max_frames.value = int(last_file[-12:-4])

    def download_svg_cb(self):
        file_str = os.path.join(self.output_dir, '*.svg')
        # print('zip up all ',file_str)
        with zipfile.ZipFile('svg.zip', 'w') as myzip:
            for f in glob.glob(file_str):
                myzip.write(f, os.path.basename(f))   # 2nd arg avoids full filename path in the archive

    def download_cb(self):
        file_xml = os.path.join(self.output_dir, '*.xml')
        file_mat = os.path.join(self.output_dir, '*.mat')
        # print('zip up all ',file_str)
        with zipfile.ZipFile('mcds.zip', 'w') as myzip:
            for f in glob.glob(file_xml):
                myzip.write(f, os.path.basename(f)) # 2nd arg avoids full filename path in the archive
            for f in glob.glob(file_mat):
                myzip.write(f, os.path.basename(f))

    def update_max_frames(self,_b):
        self.i_plot.children[0].max = self.max_frames.value

    # called if user selected different substrate in dropdown
    def mcds_field_changed_cb(self, b):
        # print("mcds_field_changed_cb: self.mcds_field.value=",self.mcds_field.value)
        if (self.mcds_field.value == None):
            return
        self.field_index = self.mcds_field.value + 4

        field_name = self.field_dict[self.mcds_field.value]
        # print('mcds_field_changed_cb: field_name='+ field_name)
        # print(self.field_min_max[field_name])

        # BEWARE of these triggering the mcds_field_cb() callback! Hence, the "skip_cb"
        self.skip_cb = True
        self.cmap_min.value = self.field_min_max[field_name][0]
        self.cmap_max.value = self.field_min_max[field_name][1]
        self.cmap_fixed_toggle.value = bool(self.field_min_max[field_name][2])
        self.skip_cb = False

        self.i_plot.update()

    # called if user provided different min/max values for colormap, or a different colormap
    def mcds_field_cb(self, b):
        if self.skip_cb:
            return

        self.field_index = self.mcds_field.value + 4

        field_name = self.field_dict[self.mcds_field.value]
        # print('mcds_field_cb: field_name='+ field_name)

        # print('mcds_field_cb: '+ field_name)
        self.field_min_max[field_name][0] = self.cmap_min.value 
        self.field_min_max[field_name][1] = self.cmap_max.value
        self.field_min_max[field_name][2] = self.cmap_fixed_toggle.value

#        self.field_index = self.mcds_field.value + 4
#        print('field_index=',self.field_index)
        self.i_plot.update()

    #------------------------------------------------------------
    def update_custom_data(self,b):
        print('----- update_custom_data')
        self.plot_cell_custom_data()

    #------------------------------------------------------------
    def plot_cell_custom_data(self):
        print('----- plot_cell_custom_data()')
        x = np.linspace(0, 2*np.pi, 400)
        y = np.sin(x**2)
        # self.i_plot.update()
        self.ax1.plot(x, y)

    #------------------------------------------------------------
    # def plot2D_custom_data(self, frame):
    def plot_cell_custom_data_0(self, xname, yname_list, t):
        # global current_idx, axes_max
        global current_frame
        # current_frame = frame
        # fname = "snapshot%08d.svg" % frame
        # full_fname = os.path.join(self.output_dir, fname)
        # print('plot_cell_custom_data: self.output_dir=',self.output_dir)

        cwd = os.getcwd()
        # print('plot_cell_custom_data: cwd=',cwd)
        # if not 'tmpdir' in cwd:
        tdir = os.path.abspath('tmpdir')
        os.chdir(tdir)

        xml_files = glob.glob('output*.xml')
        # xml_files = glob.glob(os.path.join('tmpdir', 'output*.xml'))
        xml_files.sort()
        # print('plot_cell_custom_data: xml_files[0]=',xml_files[0])
        os.chdir(cwd)

        ds_count = len(xml_files)
        # print('plot_cell_custom_data: ds_count=',ds_count)
        # mcds = [pyMCDS(xml_files[i], '.') for i in range(ds_count)]
        mcds = [pyMCDS(xml_files[i], 'tmpdir') for i in range(ds_count)]
        # print('plot_cell_custom_data: mcds[0]=',mcds[0])

        # def cell_data_plot(xname, yname_list, t):
        tname = "time"
        discrete_cells_names = ['virion', 'assembled_virion']
        tval = np.linspace(0, mcds[-1].get_time(), len(xml_files))
        # return

        if xname == tname:
            xval = tval
        elif xname in discrete_cells_names:
            xval = np.array([mcds[i].data['discrete_cells'][xname].sum() for i in range(ds_count)])
        else:
            if xname == 'susceptible_cells':
                xval = np.array([(mcds[i].data['discrete_cells']['assembled_virion'] <= 1).sum() for i in range(ds_count)])
                + np.array([(mcds[i].data['discrete_cells']['cycle_model'] < 6).sum() for i in range(ds_count)])
            elif xname == 'infected_cells':
                xval = np.array([(mcds[i].data['discrete_cells']['assembled_virion'] > 1).sum() for i in range(ds_count)]) \
                + np.array([(mcds[i].data['discrete_cells']['cycle_model'] < 6).sum() for i in range(ds_count)])
            elif xname == 'dead_cells':
                xval = np.array([len(mcds[0].data['discrete_cells']['ID']) - len(mcds[i].data['discrete_cells']['ID']) for i in range(ds_count)]) \
                + np.array([(mcds[i].data['discrete_cells']['cycle_model'] >= 6).sum() for i in range(ds_count)])

        for yname in yname_list:
            if yname in discrete_cells_names:
                yval = np.array([mcds[i].data['discrete_cells'][yname].sum() for i in range(ds_count)])
            else:
                if yname == 'susceptible_cells':
                    yval = np.array([(mcds[i].data['discrete_cells']['assembled_virion'] <= 1).sum() for i in range(ds_count)])
                    + np.array([(mcds[i].data['discrete_cells']['cycle_model'] < 6).sum() for i in range(ds_count)])
                elif yname == 'infected_cells':
                    yval = np.array([(mcds[i].data['discrete_cells']['assembled_virion'] > 1).sum() for i in range(ds_count)])
                    + np.array([(mcds[i].data['discrete_cells']['cycle_model'] < 6).sum() for i in range(ds_count)])
                elif yname == 'dead_cells':
                    yval = np.array([len(mcds[0].data['discrete_cells']['ID']) - len(mcds[i].data['discrete_cells']['ID']) for i in range(ds_count)]) \
                    + np.array([(mcds[i].data['discrete_cells']['cycle_model'] >= 6).sum() for i in range(ds_count)])
            p = self.ax1.plot(xval, yval, label=yname)
            # print('xval=',xval)  # [   0.   60.  120. ...
            # print('yval=',yval)  # [2793 2793 2793 ...
            # print('t=',t)
            if (t >= 0):
                self.ax1.plot(xval[t], yval[t], p[-1].get_color(), marker='o')
            # self.ax1.gca().spines['top'].set_visible(False)
            # self.ax1.gca().spines['right'].set_visible(False)
            # self.ax1.margins(0)

        if xname == tname:
            self.ax1.set_xlabel('time (min)')
        else:
            self.ax1.set_xlabel('total ' * (xname != tname) + xname)
        self.ax1.set_ylabel('total ' + (yname_list[0] if len(yname_list) == 1 else ', '.join(yname_list)))
        # self.ax1.set_legend()
        # self.ax1.tight_layout()
        # self.ax1.show()


    #---------------------------------------------------------------------------
    def circles(self, x, y, s, c='b', vmin=None, vmax=None, **kwargs):
        """
        See https://gist.github.com/syrte/592a062c562cd2a98a83 

        Make a scatter plot of circles. 
        Similar to plt.scatter, but the size of circles are in data scale.
        Parameters
        ----------
        x, y : scalar or array_like, shape (n, )
            Input data
        s : scalar or array_like, shape (n, ) 
            Radius of circles.
        c : color or sequence of color, optional, default : 'b'
            `c` can be a single color format string, or a sequence of color
            specifications of length `N`, or a sequence of `N` numbers to be
            mapped to colors using the `cmap` and `norm` specified via kwargs.
            Note that `c` should not be a single numeric RGB or RGBA sequence 
            because that is indistinguishable from an array of values
            to be colormapped. (If you insist, use `color` instead.)  
            `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
        vmin, vmax : scalar, optional, default: None
            `vmin` and `vmax` are used in conjunction with `norm` to normalize
            luminance data.  If either are `None`, the min and max of the
            color array is used.
        kwargs : `~matplotlib.collections.Collection` properties
            Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
            norm, cmap, transform, etc.
        Returns
        -------
        paths : `~matplotlib.collections.PathCollection`
        Examples
        --------
        a = np.arange(11)
        circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
        plt.colorbar()
        License
        --------
        This code is under [The BSD 3-Clause License]
        (http://opensource.org/licenses/BSD-3-Clause)
        """

        if np.isscalar(c):
            kwargs.setdefault('color', c)
            c = None

        if 'fc' in kwargs:
            kwargs.setdefault('facecolor', kwargs.pop('fc'))
        if 'ec' in kwargs:
            kwargs.setdefault('edgecolor', kwargs.pop('ec'))
        if 'ls' in kwargs:
            kwargs.setdefault('linestyle', kwargs.pop('ls'))
        if 'lw' in kwargs:
            kwargs.setdefault('linewidth', kwargs.pop('lw'))
        # You can set `facecolor` with an array for each patch,
        # while you can only set `facecolors` with a value for all.

        zipped = np.broadcast(x, y, s)
        patches = [Circle((x_, y_), s_)
                for x_, y_, s_ in zipped]
        collection = PatchCollection(patches, **kwargs)
        if c is not None:
            c = np.broadcast_to(c, zipped.shape).ravel()
            collection.set_array(c)
            collection.set_clim(vmin, vmax)

        # ax = plt.gca()
        # ax.add_collection(collection)
        # ax.autoscale_view()
        self.ax0.add_collection(collection)
        self.ax0.autoscale_view()
        # plt.draw_if_interactive()
        if c is not None:
            # plt.sci(collection)
            self.ax0.sci(collection)
        # return collection

    #------------------------------------------------------------
    # def plot_svg(self, frame, rdel=''):
    def plot_svg(self, frame):
        # global current_idx, axes_max
        global current_frame
        current_frame = frame
        fname = "snapshot%08d.svg" % frame
        full_fname = os.path.join(self.output_dir, fname)
        # with debug_view:
            # print("plot_svg:", full_fname) 
        # print("-- plot_svg:", full_fname) 
        if not os.path.isfile(full_fname):
            print("Once output files are generated, click the slider.")   
            return

        xlist = deque()
        ylist = deque()
        rlist = deque()
        rgb_list = deque()

        #  print('\n---- ' + fname + ':')
#        tree = ET.parse(fname)
        tree = ET.parse(full_fname)
        root = tree.getroot()
        #  print('--- root.tag ---')
        #  print(root.tag)
        #  print('--- root.attrib ---')
        #  print(root.attrib)
        #  print('--- child.tag, child.attrib ---')
        numChildren = 0
        for child in root:
            #    print(child.tag, child.attrib)
            #    print("keys=",child.attrib.keys())
            if self.use_defaults and ('width' in child.attrib.keys()):
                self.axes_max = float(child.attrib['width'])
                # print("debug> found width --> axes_max =", axes_max)
            if child.text and "Current time" in child.text:
                svals = child.text.split()
                # remove the ".00" on minutes
                self.title_str += "   cells: " + svals[2] + "d, " + svals[4] + "h, " + svals[7][:-3] + "m"

                # self.cell_time_mins = int(svals[2])*1440 + int(svals[4])*60 + int(svals[7][:-3])
                # self.title_str += "   cells: " + str(self.cell_time_mins) + "m"   # rwh

            # print("width ",child.attrib['width'])
            # print('attrib=',child.attrib)
            # if (child.attrib['id'] == 'tissue'):
            if ('id' in child.attrib.keys()):
                # print('-------- found tissue!!')
                tissue_parent = child
                break

        # print('------ search tissue')
        cells_parent = None

        for child in tissue_parent:
            # print('attrib=',child.attrib)
            if (child.attrib['id'] == 'cells'):
                # print('-------- found cells, setting cells_parent')
                cells_parent = child
                break
            numChildren += 1

        num_cells = 0
        #  print('------ search cells')
        for child in cells_parent:
            #    print(child.tag, child.attrib)
            #    print('attrib=',child.attrib)
            for circle in child:  # two circles in each child: outer + nucleus
                #  circle.attrib={'cx': '1085.59','cy': '1225.24','fill': 'rgb(159,159,96)','r': '6.67717','stroke': 'rgb(159,159,96)','stroke-width': '0.5'}
                #      print('  --- cx,cy=',circle.attrib['cx'],circle.attrib['cy'])
                xval = float(circle.attrib['cx'])

                # map SVG coords into comp domain
                # xval = (xval-self.svg_xmin)/self.svg_xrange * self.x_range + self.xmin
                xval = xval/self.x_range * self.x_range + self.xmin

                s = circle.attrib['fill']
                # print("s=",s)
                # print("type(s)=",type(s))
                if (s[0:3] == "rgb"):  # if an rgb string, e.g. "rgb(175,175,80)" 
                    rgb = list(map(int, s[4:-1].split(",")))  
                    rgb[:] = [x / 255. for x in rgb]
                else:     # otherwise, must be a color name
                    rgb_tuple = mplc.to_rgb(mplc.cnames[s])  # a tuple
                    rgb = [x for x in rgb_tuple]

                # test for bogus x,y locations (rwh TODO: use max of domain?)
                too_large_val = 10000.
                if (np.fabs(xval) > too_large_val):
                    print("bogus xval=", xval)
                    break
                yval = float(circle.attrib['cy'])
                # yval = (yval - self.svg_xmin)/self.svg_xrange * self.y_range + self.ymin
                yval = yval/self.y_range * self.y_range + self.ymin
                if (np.fabs(yval) > too_large_val):
                    print("bogus xval=", xval)
                    break

                rval = float(circle.attrib['r'])
                # if (rgb[0] > rgb[1]):
                #     print(num_cells,rgb, rval)
                xlist.append(xval)
                ylist.append(yval)
                rlist.append(rval)
                rgb_list.append(rgb)

                # For .svg files with cells that *have* a nucleus, there will be a 2nd
                if (not self.show_nucleus):
                #if (not self.show_nucleus):
                    break

            num_cells += 1

            # if num_cells > 3:   # for debugging
            #   print(fname,':  num_cells= ',num_cells," --- debug exit.")
            #   sys.exit(1)
            #   break

            # print(fname,':  num_cells= ',num_cells)

        xvals = np.array(xlist)
        yvals = np.array(ylist)
        rvals = np.array(rlist)
        rgbs = np.array(rgb_list)
        # print("xvals[0:5]=",xvals[0:5])
        # print("rvals[0:5]=",rvals[0:5])
        # print("rvals.min, max=",rvals.min(),rvals.max())

        # rwh - is this where I change size of render window?? (YES - yipeee!)
        #   plt.figure(figsize=(6, 6))
        #   plt.cla()
        # if (self.substrates_toggle.value):
        self.title_str += " (" + str(num_cells) + " agents)"
            # title_str = " (" + str(num_cells) + " agents)"
        # else:
            # mins= round(int(float(root.find(".//current_time").text)))  # TODO: check units = mins
            # hrs = int(mins/60)
            # days = int(hrs/24)
            # title_str = '%dd, %dh, %dm' % (int(days),(hrs%24), mins - (hrs*60))
        # plt.title(self.title_str)
        self.ax0.set_title(self.title_str)

        # plt.xlim(self.xmin, self.xmax)
        # plt.ylim(self.ymin, self.ymax)
        self.ax0.set_xlim(self.xmin, self.xmax)
        self.ax0.set_ylim(self.ymin, self.ymax)

        #   plt.xlim(axes_min,axes_max)
        #   plt.ylim(axes_min,axes_max)
        #   plt.scatter(xvals,yvals, s=rvals*scale_radius, c=rgbs)

        # TODO: make figsize a function of plot_size? What about non-square plots?
        # self.fig = plt.figure(figsize=(9, 9))

#        axx = plt.axes([0, 0.05, 0.9, 0.9])  # left, bottom, width, height
#        axx = fig.gca()
#        print('fig.dpi=',fig.dpi) # = 72

        #   im = ax.imshow(f.reshape(100,100), interpolation='nearest', cmap=cmap, extent=[0,20, 0,20])
        #   ax.xlim(axes_min,axes_max)
        #   ax.ylim(axes_min,axes_max)

        # convert radii to radii in pixels
        # ax1 = self.fig.gca()
        # N = len(xvals)
        # rr_pix = (ax1.transData.transform(np.vstack([rvals, rvals]).T) -
        #             ax1.transData.transform(np.vstack([np.zeros(N), np.zeros(N)]).T))
        # rpix, _ = rr_pix.T

        # markers_size = (144. * rpix / self.fig.dpi)**2   # = (2*rpix / fig.dpi * 72)**2
        # markers_size = markers_size/4000000.
        # print('max=',markers_size.max())

        #rwh - temp fix - Ah, error only occurs when "edges" is toggled on
        if (self.show_edge):
            try:
                # plt.scatter(xvals,yvals, s=markers_size, c=rgbs, edgecolor='black', linewidth=0.5)
                self.circles(xvals,yvals, s=rvals, color=rgbs, edgecolor='black', linewidth=0.5)
                # cell_circles = self.circles(xvals,yvals, s=rvals, color=rgbs, edgecolor='black', linewidth=0.5)
                # plt.sci(cell_circles)
            except (ValueError):
                pass
        else:
            # plt.scatter(xvals,yvals, s=markers_size, c=rgbs)
            self.circles(xvals,yvals, s=rvals, color=rgbs)

        # x = np.linspace(0, 2*np.pi, 100)
        # y = np.sin(x**2)
        # self.i_plot.update()
        # self.ax1.plot(x, y)
        self.plot_cell_custom_data_0("time", ["assembled_virion"], 20)

        # if (self.show_tracks):
        #     for key in self.trackd.keys():
        #         xtracks = self.trackd[key][:,0]
        #         ytracks = self.trackd[key][:,1]
        #         plt.plot(xtracks[0:frame],ytracks[0:frame],  linewidth=5)


    #---------------------------------------------------------------------------
    # assume "frame" is cell frame #, unless Cells is togggled off, then it's the substrate frame #
    # def plot_substrate(self, frame, grid):
    def plot_substrate(self, frame):

        # print("plot_substrate(): frame*self.substrate_delta_t  = ",frame*self.substrate_delta_t)
        # print("plot_substrate(): frame*self.svg_delta_t  = ",frame*self.svg_delta_t)
        self.title_str = ''

        # Recall:
        # self.svg_delta_t = config_tab.svg_interval.value
        # self.substrate_delta_t = config_tab.mcds_interval.value
        # self.modulo = int(self.substrate_delta_t / self.svg_delta_t)
        # self.therapy_activation_time = user_params_tab.therapy_activation_time.value

        # print("plot_substrate(): pre_therapy: max svg, substrate frames = ",max_svg_frame_pre_therapy, max_substrate_frame_pre_therapy)

        # Assume: # .svg files >= # substrate files
#        if (self.cells_toggle.value):

        if (self.substrates_toggle.value):
            # maybe only show 2nd plot if self.custom_data_toggle is True
            self.fig, (self.ax0, self.ax1) = plt.subplots(1, 2, figsize=(30, 12))


            if (self.customized_output_freq and (frame > self.max_svg_frame_pre_therapy)):
                self.substrate_frame = self.max_substrate_frame_pre_therapy + (frame - self.max_svg_frame_pre_therapy)
            else:
                self.substrate_frame = int(frame / self.modulo)

            fname = "output%08d_microenvironment0.mat" % self.substrate_frame
            xml_fname = "output%08d.xml" % self.substrate_frame
            # fullname = output_dir_str + fname

    #        fullname = fname
            full_fname = os.path.join(self.output_dir, fname)
            # print("--- plot_substrate(): full_fname=",full_fname)
            full_xml_fname = os.path.join(self.output_dir, xml_fname)
    #        self.output_dir = '.'

    #        if not os.path.isfile(fullname):
            if not os.path.isfile(full_fname):
                print("Once output files are generated, click the slider.")  # No:  output00000000_microenvironment0.mat
                return

    #        tree = ET.parse(xml_fname)
            tree = ET.parse(full_xml_fname)
            xml_root = tree.getroot()
            mins = round(int(float(xml_root.find(".//current_time").text)))  # TODO: check units = mins
            self.substrate_mins= round(int(float(xml_root.find(".//current_time").text)))  # TODO: check units = mins

            hrs = int(mins/60)
            days = int(hrs/24)
            self.title_str = 'substrate: %dd, %dh, %dm' % (int(days),(hrs%24), mins - (hrs*60))
            # self.title_str = 'substrate: %dm' % (mins )   # rwh

            info_dict = {}
            scipy.io.loadmat(full_fname, info_dict)
            M = info_dict['multiscale_microenvironment']
            f = M[self.field_index, :]   # 4=tumor cells field, 5=blood vessel density, 6=growth substrate

            try:
                xgrid = M[0, :].reshape(self.numy, self.numx)
                ygrid = M[1, :].reshape(self.numy, self.numx)
            except:
                print("substrates.py: mismatched mesh size for reshape: numx,numy=",self.numx, self.numy)
                pass
#                xgrid = M[0, :].reshape(self.numy, self.numx)
#                ygrid = M[1, :].reshape(self.numy, self.numx)

            num_contours = 15
            levels = MaxNLocator(nbins=num_contours).tick_values(self.cmap_min.value, self.cmap_max.value)
            contour_ok = True
            if (self.cmap_fixed_toggle.value):
                try:
                    substrate_plot = self.ax0.contourf(xgrid, ygrid, M[self.field_index, :].reshape(self.numy, self.numx), levels=levels, extend='both', cmap=self.field_cmap.value, fontsize=self.fontsize)
                except:
                    contour_ok = False
                    # print('got error on contourf 1.')
            else:    
                try:
                    substrate_plot = self.ax0.contourf(xgrid, ygrid, M[self.field_index, :].reshape(self.numy,self.numx), num_contours, cmap=self.field_cmap.value)
                except:
                    contour_ok = False
                    # print('got error on contourf 2.')

            if (contour_ok):
                self.ax0.set_title(self.title_str, fontsize=self.fontsize)
                cbar = self.fig.colorbar(substrate_plot, ax=self.ax0)
                cbar.ax.tick_params(labelsize=self.fontsize)

            self.ax0.set_xlim(self.xmin, self.xmax)
            self.ax0.set_ylim(self.ymin, self.ymax)

        # Now plot the cells (possibly on top of the substrate)
        if (self.cells_toggle.value):
            if (not self.substrates_toggle.value):
                # maybe only show 2nd plot if self.custom_data_toggle is True
                # self.fig, (self.ax0, self.ax1) = plt.subplots(1, 2, figsize=(self.figsize_width_svg*2, self.figsize_height_svg))
                # self.fig, (self.ax0, self.ax1) = plt.subplots(1, 2, figsize=(24, 12))
                self.fig, (self.ax0, self.ax1) = plt.subplots(1, 2, figsize=(20, 10))

            self.svg_frame = frame
            # print('plot_svg with frame=',self.svg_frame)
            self.plot_svg(self.svg_frame)


            if not self.custom_data_plotted:
                self.plot_cell_custom_data()
                # self.plot_cell_custom_data('time', ['susceptible_cells', 'infected_cells', 'dead_cells'], 20)
                self.custom_data_plotted = True
