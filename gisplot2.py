#package imports
import hublib.use
import os,sys
import io
# On VHub, xarray and cartopy installed in anaconda2-5.1

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")

import numpy as np
import xarray as xr

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#%matplotlib inline

import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

sys.path.insert(0, os.path.abspath('..'))
import hublib.ui as ui

import ipywidgets as widgets
import math

import requests
from requests.exceptions import HTTPError

from IPython.display import Image 
import time

from collections import Counter

from IPython.display import Javascript

import threading
import cProfile

def load_icesheet_data(file_names):
    #update_progress function is used to show the loading animation while running plotting functions
    update_progress(0)
    #compiling file names into a list
    counter = 0
    global file_vars
    file_vars = []
    for file in file_names:
        file_vars.append(file)
        counter = counter+1
    
    global mapping_var
    if 'AIS' in file_vars[0]:
        mapping_var = 'mapping'
    elif 'GIS' in file_vars[0]:
        mapping_var = 'Polar_Stereographic'
        
    #converting netCDF files into Xarrays and storing them in a list
    counter = 0
    global model_vars
    model_vars = []  
    for f in file_vars:
        #decoded_times = False is needed to process these files
        #As far as I can tell it is only needed for AIS data
        model_vars.append(xr.open_dataset(f,engine='netcdf4',decode_times=False)) 
        counter = counter+1
    
    update_progress(0.1)  
    #collect titles for each of the subplots using the file names
    global ctvs, titles
    ctvs = []
    titles = []
    for f in file_names:
        fbase = os.path.basename(f)
        t = fbase.replace('.nc','')
        lt = t.split('_')
        #try catch to catch error in case the 2D variable name does not exist in the file
        try:
            #ctvs lilst stores the 2D variable names      ##might no longer be needed after adding references
            ctvs.append(lt[0])
        except:
            with output_widget:
                print('Data Variable missing from netCDF file')
            #"return" returns None, which simply exits the function
            return
        title = lt[2]+'_'+lt[3]
        titles.append(title)
    
    #function grabs and stores polar stereographic data from each Xarray 
    counter = 0
    global ctv_vars, ctv_proj_vars
    ctv_vars = []
    ctv_proj_vars = []
    try:
        for m in model_vars:
            ctv_vars.append(m[ctvs[counter]])
            ctv_proj_vars.append(m[mapping_var])
            counter = counter+1
    except:
        with output_widget:
            print('')
        
    update_progress(0.2) 
    
    #sets the stereographic projection based on the user's choice between a model based and standard projection
    global polar_stereographic
    if stereograph_choice.value == 'standard':
        if mapping_var == 'Polar_Stereographic':
            #Set Standard Polar Sterographic Projection definition
            polar_stereographic = ccrs.Stereographic(
                central_latitude=90.0,
                central_longitude=-45.0,
                false_easting=0.0,
                false_northing=0.0,
                true_scale_latitude=70.0,
                globe=ccrs.Globe('WGS84'))
        elif mapping_var == 'mapping':
            polar_stereographic = ccrs.Stereographic(
                central_latitude=-90.0,
                central_longitude=0.0,
                false_easting=0.0,
                false_northing=0.0,
                true_scale_latitude=-71.0,
                globe=ccrs.Globe('WGS84'))
    
    else:
        polar_stereographic = ccrs.Stereographic(
            central_latitude=ctv_proj_vars[0].latitude_of_projection_origin,
            central_longitude=ctv_proj_vars[0].straight_vertical_longitude_from_pole,
            false_easting=ctv_proj_vars[0].false_easting,
            false_northing=ctv_proj_vars[0].false_northing,
            true_scale_latitude=ctv_proj_vars[0].standard_parallel,
            globe=ccrs.Globe('WGS84') )
            
    update_progress(0.3)
    
    #calls the next function in the plotting process
    try:
        check = transform_projection()
        if check=='failed':
            upload_button.clear_output()
            return check
    except:
        with output_widget:
            print('ERROR: Projection Transformation Failed')
        return
        
def transform_projection(): 
    ######################
    # Transform projection
    #setting cartopy map values based on WGS84
    geodetic = ccrs.Geodetic(globe=ccrs.Globe('WGS84'))
    
    yv, xv = np.meshgrid(model_vars[0].y.data, model_vars[0].x.data)
    ll = geodetic.transform_points(src_crs=polar_stereographic, x=xv.flatten(), y=yv.flatten())
    global lons,lats
    lons = ll[:,0]
    lats = ll[:,1]
    
    update_progress(0.4)
    
    counter = 0
    global ctv_mean_vars
    ctv_mean_vars = []
    for l in ctv_vars:
        ctv_mean_vars.append(l.mean(dim='time').data)
        ctv_mean_vars[counter] = ctv_mean_vars[counter].transpose()
        ctv_mean_vars[counter] = ctv_mean_vars[counter].flatten()
        counter = counter+1
        
    update_progress(0.5)
    
    try:
        check = plot_icesheet_data()
        if check == 'failed':
            return check
    except:
        output_widget.clear_output()
        with output_widget:
            print('ERROR: Plotting Failed')
        time.sleep(3)
        return 'failed'
    
    #line not needed just used to profile the following function
#    cProfile.run('plot_icesheet_data()')

def plot_icesheet_data():    
    ####################
    #Plot Transformed Ice Sheet Data
    update_progress(0.6)
    
    #get dimensions for subplotting images
    num_frames = len(ctv_vars)
    if num_frames == 1:
        grid_cols = 1
        grid_rows = 1
        fig_dims = [10,10]
    else:
        grid_cols = 2
        grid_rows = math.ceil(math.sqrt(num_frames))
        fig_dims = [16,8*grid_rows]
    #creates a figure using the dimensions calculated
    plt.figure(figsize=(fig_dims[0],fig_dims[1]))
    
    update_progress(0.7)
    #using all previously gathered data to create the subplots for each dataset
    try:
        counter = 0
        ax_vars = []
        for l in ctv_mean_vars:
            frame = counter+1
            ax_vars.append(plt.subplot(grid_rows,grid_cols,frame, projection=polar_stereographic))
            if mapping_var == 'Polar_Stereographic':
                ax_vars[counter].set_extent([-65, -20, 57, 84]) #not needed for ais plots
#        elif mapping_var == 'mapping':
#            pass
            #ax_vars[counter].set_extent([-180,-160,183,30])#shortest = [-180,-160,180,-10]
            #if this doesn't end up speeding it up just delete
            #ax_vars[counter].set_extent([-65, -20, 57, 84]) #changing these values may speed up the plotting
            ax_vars[counter].coastlines(resolution='10m', zorder=7)
            ax_vars[counter].gridlines(zorder=8)
        #appropriately names the reference data subplot
            if counter==0:
                ax_vars[counter].set_title('Reference Data ('+titles[0]+')', fontsize=20)
            else:
                ax_vars[counter].set_title(titles[counter], fontsize=20)
#        if stereograph_choice.value == 'standard':
            plt.scatter(lons, lats, 1, c=ctv_mean_vars[counter], transform=ccrs.Geodetic(), zorder=0, cmap='viridis')
#        else:
#            plt.scatter(lons_vars[counter], lats_vars[counter], 1, c=ctv_mean_vars[counter], transform=ccrs.Geodetic(), zorder=0, cmap='viridis')
        #sets subplots axis name and units based on the data in the Xarray
            data_table = getattr(model_vars[counter],ctvs[counter])
            name = data_table.attrs['standard_name']
            units = data_table.attrs['units']
            c = plt.colorbar(fraction=0.046, pad=0.04)
            c.set_label('{0} ({1})'.format(name, units), size=16)
            counter = counter+1
        axes = plt.gca()
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    except:
        return 'failed'

    update_progress(0.8)
    #sets title for entire figure
    plt.suptitle((name+' ('+ctvs[0]+')'), fontsize=30) #subtitle name may need to be adjusted later
    plt.subplots_adjust(top=0.88)
    #saves figure for user to download later if they choose
    plt.savefig('GIS_Ice_Sheet_Model_Comparison.png', dpi='figure')
    
    update_progress(0.9)
    #print(ax_vars[0].get_extent())
    output_widget.clear_output(wait = True)
    with output_widget:
        plt.show() 
    with filename_output:
        print('Plotted Files:')
        for f in file_names:
            fbase = os.path.basename(f)
            print(fbase)
def check_files():
    global checked_files, wrong_files
    checked_files = []
    wrong_files = []
    var_names = []
    ex_names = []
    split_names = []
    counter = 0
    #loop to collect 2D variables and experiment names from the provided files
    for f in uploaded_filenames:
        fbase = os.path.basename(f).replace('.nc','')
        split_names.append(fbase.split('_'))
        var_names.append(split_names[counter][0])
        ex_names.append(split_names[counter][4])
        counter = counter+1
    #these two lines find the most commonly occuring 2D variable in the files provided
    variables = Counter(var_names)
    ref_var = variables.most_common(1)[0][0]
    counter = 0
    exs = []
    #loop uses most common 2D variable to find the most common experiment
    for f in uploaded_filenames:
        if ref_var == split_names[counter][0]:
            exs.append(ex_names[counter])
        counter = counter+1
    experiments = Counter(exs)
    ref_ex = experiments.most_common(1)[0][0]
    counter = 0
    #puts the files fitting the reference file into checked files and saves the removed files into wrong files
    for f in uploaded_filenames:
        if ref_var == split_names[counter][0] and ref_ex == split_names[counter][4]:
            checked_files.append(f)
        else:
            wrong_files.append(f)
        counter = counter+1
    ref_file = ''
    #creates file path based on if the user is uploading Greenland or Antarctic Ice Sheet Data
    if 'GIS' in checked_files[0]:
        path = '/data/groups/ghub/tools/reference/gis'
    elif 'AIS' in checked_files[0]:
        path = '/data/groups/ghub/tools/reference/ais'
    #r=root d=directory f=file
    for r,d,f in os.walk(path):
        for file in f:
            rbase = file.replace('.nc','')
            split_ref = rbase.split('_')
            if split_ref[0] == ref_var and split_ref[4] == ref_ex:
                ref_file = os.path.join(r,file)
    #checks in case a rerence file couldn't be found that fits the uploaded data
    if len(ref_file)==0:
        with output_widget:
            print('ERROR: No reference file found')
            return
    checked_files.insert(0,ref_file)
    return 'all good'
    
def uploaded_cb(b):
    download_button(False)
    #d1.w.layout.visibility = 'hidden'
    fm.visible = False
    stereograph_choice.disabled = True
    uploaded_data.disabled = True
    
    global file_names
    file_names = checked_files
    try:
        check = load_icesheet_data(file_names)
    except:
        output_widget.clear_output()
        fm.visible = True
        with output_widget:
            print('Plotting Failed')
            print('File may be formatted incorrectly')
    
    global projection_value
    projection_value = stereograph_choice.value
    
    uploaded_data.disabled = False
    stereograph_choice.disabled = False
    button_output_widget.clear_output()
    fm.visible = True
    download_button(True)
    with trash_output:
        clear_uploads()
    
    if check=='failed':
        download_button(False)
        fm.visible = False
        stereograph_choice.disabled = True
        uploaded_data.disabled = True
        with output_widget:
            print('ERROR: Plotting Failed')
            print('File may be formatted incorrectly')
        fm.visible = True
    
# Called when all files finish uploading
def done_cb(w, name):
    upload_button.clear_output()
    filename_output.clear_output()
    trash_output.clear_output()
    output_widget.clear_output()
    global projection_value
    projection_value = ''
    global uploaded_filenames
    uploaded_filenames = []
    #checks to make sure files are netCDFs
    for n in name:
        fbase = os.path.basename(n)
        if '.nc' in n:
            uploaded_filenames.append(n)
        else:
            with output_widget:
                print('File '+fbase+' Not Uploaded')
                print('Incorrect File Type: must be a netcdf(.nc) file\n')
    #resets upload if no vallid netCDFs are uploaded
    if len(uploaded_filenames)==0:
        w.reset()
        uploaded_data.layout.visibility = 'hidden'
        stereograph_choice.layout.visibility = 'hidden'
        return
    #calls function to check the files and get corresponding reference file
    check = check_files()
    #resets upload if the check files fnction throws an error
    if check == None:
        w.reset()
        uploaded_data.layout.visibility = 'hidden'
        stereograph_choice.layout.visibility = 'hidden'
        return
    uploaded_data.layout.visibility = 'visible'
    stereograph_choice.layout.visibility = 'visible'
    #prints info about what was uploaded to the user
    with output_widget:
        c = 0
        for f in checked_files:
            fbase = os.path.basename(f)
            if c==0:
                print("Reference File: %s" % fbase)
            else:
                print("%s uploaded" % fbase)
            c=c+1
        if len(wrong_files)!=0:
            print('\n\nSome files did not fit the chosen reference file:')
            for f in wrong_files:
                fbase = os.path.basename(f)
                print(fbase)
    reset_upload.layout.visibility = 'hidden'
    # reset clears and re-enables the widget for more uploads
    w.reset()
    display_upload()

def update_progress(progress):
    #generates loading animation throughout plotting process
    title = 'Plotting Uploaded Data:'
    bar_length = 20
    block = int(20.0*progress)
    text = title+" [{0}] {1:.1f}%".format( "#" * block + "-" * (bar_length - block), progress * 100)
    output_widget.clear_output(wait = True)
    with output_widget:
        print(text)
        
def download_button(active):
    ##creates download button when called
    #prevents pop ups when the tool is launched
    if active:
        d1 = ui.Download('GIS_Ice_Sheet_Model_Comparison.png', label = 'Download Plot', 
                  style='success', icon='fa-arrow-circle-down')
        with button_output_widget:
            display(d1)
    else:
        button_output_widget.clear_output()       
        
def button_deactivation(change):
    #deactivates or activates buttons based on when the upload button is clicked
    if len(uploaded_filenames)==0:
        return
    if len(change['old'])==1:
        output_widget.clear_output()
        download_button(False)
        #d1.w.layout.visibility = 'hidden'
        stereograph_choice.disabled = True
        uploaded_data.disabled = True
    elif len(change['old'])>1:
        stereograph_choice.disabled = False
        uploaded_data.disabled = False

def check_permissions():
    #checks to see if user has permission to upload files and read from reference directory
    f = open('output.txt','w+r')
    !getfacl -d /data/groups/ghub/tools/gisplot2 >> output.txt
    for l in f:
        if 'group' in l:
            access = l
    if 'rwx' not in access:
        try:
            !setfacl -d -m g::rwx /data/groups/ghub/tools/gisplot2
        except:
            with output_widget:
                print('Upload Permission Restricted: Contact Admnistrator')
        
def second_plot_activation(change):
    #allows for dynamic adtivation of stereograph_choice widget
    if projection_value != '':
        if projection_value == stereograph_choice.value:
            uploaded_data.disabled = True
        else:
            uploaded_data.disabled = False
    else:
        uploaded_data.disabled = False
        
def clear_uploads():
    #deletes upload directory when needed
    !rm -rf {upload_directory}

def provide_fix(change):
    filename_output.clear_output()
    if len(change['new'])>len(change['old']):
        reset_upload.layout.visibility = 'visible'
    else:
        reset_upload.layout.visibility = 'hidden'

def fix_upload(b):
    try:
        with trash_output:
            fm.reset()
            stereograph_choice.disabled = True
            uploaded_data.disabled = True
    except:
        with output_widget:
            print('Wait until upload finishes.')
         
#import cProfile, pstats, StringIO
#pr = cProfile.Profile()
#pr.enable()

#scripting to run when the tool is launched
#mostly creates widgets
check_permissions()

#gets session number for the user to make a unique upload directory
s = requests.session()
session = str(s).split(' ')
session = session[3].replace('>','')
upload_directory = '/data/groups/ghub/tools/gisplot2/'+session

#creating file upload widget
uploaded_filenames = []
fm = ui.FileUpload('','Please Select Your Files to Upload',
               dir= upload_directory,
               cb=done_cb,
               maxnum=6,
               maxsize='150M')
fm.w.observe(button_deactivation)
fm.w.observe(provide_fix)

#stereograph model is chosen using this widget
stereograph_choice = widgets.RadioButtons(
    description = 'Projection:',
    options = ['model based','standard'],
    disabled = False,
    layout = widgets.Layout(width = '300px',visibility = 'hidden'))
stereograph_choice.observe(second_plot_activation)

#value declared to assist in other operations later
#global projecion_value
projection_value = ''

#button used to plot data after it has been uploaded
uploaded_data = widgets.Button(
    description = 'Plot Data',
    button_style = 'info',
    disabled = False,
    layout = widgets.Layout(visibility = 'hidden'))
uploaded_data.on_click(uploaded_cb)

reset_upload = widgets.Button(
    description = 'Reset Upload',
    button_style = 'danger',
    layout = widgets.Layout(visibility = 'hidden'))
reset_upload.on_click(fix_upload)

#all output widgets
upload_button = widgets.Output()
button_output_widget = widgets.Output()
output_widget = widgets.Output()
filename_output = widgets.Output()
trash_output = widgets.Output()
#pr.disable()
#s = StringIO.StringIO()
#ps = pstats.Stats(pr, stream = s).sort_stats('cumulative')
#ps.print_stats()
#print(s.getvalue())

display(fm)

display(stereograph_choice)

def display_upload():
    with upload_button:
        display(uploaded_data)
        
display(upload_button)

reset_upload

display(button_output_widget)

display(filename_output)

display(output_widget)

#if the user exits the program the upload directory is deleted
import atexit
atexit.register(clear_uploads)

script = '''
    require(["base/js/namespace"],function(Jupyter) {
        Jupyter.notebook.save_checkpoint();
    });
    '''
Javascript(script)
