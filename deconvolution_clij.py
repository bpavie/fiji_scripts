# @File(label="Select a PSF directory", style="file") psf_file
# @File(label="Select an Image directory", style="directory") img_folder
# @File(label="Select an Output directory", style="directory") output_folder
# @String(label="Graphic card", value="RTX", description="Set to empty to pick the first card") gpu_card
# @String(label="File extension", value="nd2", description="File extension") fileExtension
# @Integer(label="Deconvolution Iterations", value=100) iterations 
# @Float(label="Deconvolution Regularization Factor", value=0.0, description="Usually, value should be betwen 0.001 and 0.01") regularizationFactor
# @String(label="Rescale intensity Range",choices={"12bits","16bits", "Min/Max", "None"}, description="Rescale the intensity (12bits:0-4095, 16bits:0-65535, Min/Max:min/max intensity of the channel, should be for visualization only) before saving it in 16bits.") rescale_intensity_method
# @boolean(label="Use in script", value=false) use_in_script
# @OpService ops
# @UIService ui


'''
Deconvolution of multi-channel stacks in batch using the  Richardson Lucy FFT implementation from clijX using the GPU and save it in a dedicated folder

See paper related to fluoresent microscopy and Richardson Lucy Algorithm
Nicolas Dey, Laure Blanc-Féraud, Christophe Zimmer, Pascal Roux, Zvi Kam, et al… 3D Microscopy Deconvolution using Richardson-Lucy Algorithm with Total Variation Regularization.
https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/jemt.20294
Here’s an open version: https://hal.inria.fr/inria-00070726/document

Data organization:
<Images_FOLDER>/a.nd2/tif/...
               /b.nd2/tif/...
<PSF_FOLDER>   psf.nd2/tif/...
Results will be saved in the result folder:
<RESULT_FOLDER>/a.tif

//Benjamin Pavie benjamin.pavie@kuleuven.vib.be
//BioImaging Core Facility Leuven
//2023/06/15
//Please acknowledge us if you use this script in a publication.

Requirements:
clijx
clij2
See https://github.com/clij/clijx
To install, Help > Update... and activate them

Run in a terminal:
Linux:
./ImageJ-linux64 --ij2 --headless --run "/home/my_user/process_deconvolution_in_batch.py" 'psf_file="/path/to/psf.nd2",img_folder="/path/to/img_folder",output_folder="/path/to/output_folder",gpu_card="", iteration=100,regularizationFactor=0.0,rescale_intensity_method="12bits"'
Windows:
ImageJ-win64.exe --ij2 --headless --run "C:\Users\my_user\process_deconvolution_in_batch.py" "psf_file='/path/to/psf.nd2',img_folder='/path/to/img_folder',output_folder='/path/to/output_folder',gpu_card='',iteration=100,regularizationFactor=0.0,rescale_intensity_method='12bits'"

//Cite the following puplication
//Robert Haase, Loic Alain Royer, Peter Steinbach, Deborah Schmidt, Alexandr Dibrov, Uwe Schmidt, Martin Weigert, Nicola Maghelli, Pavel Tomancak, Florian Jug, Eugene W Myers. CLIJ: GPU-accelerated image processing for everyone. Nat Methods 17, 5-6 (2020) doi:10.1038/s41592-019-0650-1

//If the image is to big to be processed, you could use a tiled version, see:
//https://forum.image.sc/t/clij-deconvolution-regularization-fiji-crash/55155/17?u=bpavie
'''

import os
import math
from net.haesleinhuepf.clij.coremem.enums import NativeTypeEnum;
from net.haesleinhuepf.clij import CLIJ;
from net.haesleinhuepf.clij2 import CLIJ2;
from net.haesleinhuepf.clijx import CLIJx;
from net.haesleinhuepf.clijx.plugins import DeconvolveRichardsonLucyFFT;
from ij import IJ
from ij.process import LUT
from ij.plugin import Duplicator
from ij.plugin import RGBStackMerge
from ij.gui import GenericDialog
from ij.process import ImageStatistics
from ij.process import StackStatistics
from loci.formats import ImageReader
from loci.formats import MetadataTools
from loci.plugins.in import ImporterOptions
from loci.plugins import BF
from ome.units import UNITS
from java.awt import Color 

#Filter to list only .jpg files from a directory
def acceptFile(filename, file_extension):
  """ Work only with file_extension. """
  if(len(filename) - 4 == filename.rfind('.'+file_extension)):
    return True
  else:
    return False

#Get the folder path entered by the user
#psf_folder_path=psf_folder.getAbsolutePath()      
#psf_file_list=sorted(os.listdir(psf_folder_path))

img_folder_path=img_folder.getAbsolutePath()  
#img_file_list=sorted(os.listdir(img_folder_path))
img_file_list = sorted([f for f in os.listdir(img_folder_path) if os.path.isfile(os.path.join(img_folder_path, f)) and f.endswith('.'+fileExtension)])

   
options = ImporterOptions()
options.setId(psf_file.getAbsolutePath())
options.setUngroupFiles(True)
options.setColorMode(ImporterOptions.COLOR_MODE_COLORIZED) #Open with Color default from the Original format
options.setOpenAllSeries( False ) #Just open the 1st serie, the other one is just a thumbnail to have a quick overview
imps_psf = BF.openImagePlus(options)
imp_psf = imps_psf[0] 
         
'''
reader = ImageReader()
omeMeta = MetadataTools.createOMEXMLMetadata()
reader.setMetadataStore(omeMeta)
reader.setId(psf_file.getAbsolutePath())
physSizeX = omeMeta.getPixelsPhysicalSizeX(0)
physSizeY = omeMeta.getPixelsPhysicalSizeY(0)
physSizeZ = omeMeta.getPixelsPhysicalSizeZ(0)

print('pixel x/y:'+str(physSizeX.value(UNITS.MICROM)))
print('pixel z:'+str(physSizeZ.value(UNITS.MICROM)))
channels_luts = []
'''

print("Available devices:");
for name in CLIJ.getAvailableDeviceNames():
  print(name);

'''
for nc in range(omeMeta.getChannelCount(0)):
  channel_name = omeMeta.getChannelName(0, nc)
  if channel_name:
    print('Channel name:'+omeMeta.getChannelName(0, nc))
  else:
    print('No channel name detected from the metadata')
  ome_color = omeMeta.getChannelColor(0,nc)
  #print(omeMeta.getChannelColor(0,nc))
  if(ome_color):
    javaColor = Color(ome_color.getRed(), ome_color.getGreen(), ome_color.getBlue(),
        ome_color.getAlpha())
    luts = LUT.createLutFromColor(javaColor)
    channels_luts.append(luts)
  else:
    channels_luts.append(None)
    print('No channel color detected from the metadata')
'''
            
for img_name in img_file_list:
  path = os.path.join(img_folder_path, img_name)
  if os.path.isdir(path):
    # skip directories
    continue
  IJ.log('Processing '+img_name);
  options = ImporterOptions()
  options.setId(os.path.join(img_folder_path, img_name))
  options.setUngroupFiles(True)
  options.setColorMode(ImporterOptions.COLOR_MODE_COLORIZED) #Open with Color default from the Original format
  options.setOpenAllSeries( False ) #Just open the 1st serie, the other one is just a thumbnail to have a quick overview
  imps_img = BF.openImagePlus(options)
  imp_img = imps_img[0]    
  
  reader = ImageReader()
  omeMeta = MetadataTools.createOMEXMLMetadata()
  reader.setMetadataStore(omeMeta)
  reader.setId(psf_file.getAbsolutePath())
  physSizeX = omeMeta.getPixelsPhysicalSizeX(0)
  physSizeY = omeMeta.getPixelsPhysicalSizeY(0)
  physSizeZ = omeMeta.getPixelsPhysicalSizeZ(0)
  
  channels_luts = []
  for nc in range(omeMeta.getChannelCount(0)):
    channel_name = omeMeta.getChannelName(0, nc)
    if channel_name:
      print('Channel name:'+omeMeta.getChannelName(0, nc))
    else:
      print('No channel name detected from the metadata')
    ome_color = omeMeta.getChannelColor(0,nc)
    #print(omeMeta.getChannelColor(0,nc))
    if(ome_color):
      javaColor = Color(ome_color.getRed(), ome_color.getGreen(), ome_color.getBlue(),
          ome_color.getAlpha())
      luts = LUT.createLutFromColor(javaColor)
      channels_luts.append(luts)
    else:
      channels_luts.append(None)
      print('No channel color detected from the metadata')
  
  
  
  n_channel = imp_img.getNChannels()
  n_z = imp_img.   getNSlices()
  n_f = imp_img.getNFrames()
  
  print('Number of channel:'+str(n_channel))
  deconvolution_result_list=[]
  for c in range(1,(n_channel+1)):
    imp_img_c = Duplicator().run(imp_img, c, c, 1, n_z, 1, 1);
    imp_psf_c = Duplicator().run(imp_psf, c, c, 1, n_z, 1, 1);
    
    # initialize a device with a given name
    clij2 = CLIJ2.getInstance(gpu_card)#CLIJ2.getInstance("RTX");
    clijx=CLIJx.getInstance()
    clij2.clear()
    
    gpuPSF = clij2.push(imp_psf_c);
    
    #deconvolved = IJ.createImage("deconvolved", "32-bit", imp_img_c.getWidth(), imp_img_c.getHeight(), imp_img_c.getNSlices());
    
    gpuImg = clij2.push(imp_img_c);
    tempOut = clij2.create(gpuImg.getDimensions(), NativeTypeEnum.Float);    
    DeconvolveRichardsonLucyFFT.deconvolveRichardsonLucyFFT(clij2, gpuImg, gpuPSF, tempOut, iterations, regularizationFactor, False);
    deconvolved = clij2.pull(tempOut);
    
    #Rescale the intensity to 12/16bits range or min/max
    min_i = 0
    max_i = pow(2, 16)-1
    
    if rescale_intensity_method == "12bits":
      max_i = pow(2, 12)-1
    elif rescale_intensity_method == "Min/Max":
      stats2 = StackStatistics(deconvolved)
      min_i = stats2.min
      max_i = stats2.max
    
    IJ.log("Rescale intensity of channel "+str(c)+" to range : "+str(min_i)+"-"+str(max_i))
    deconvolved.setDisplayRange(min_i, max_i, 1)
    
    if rescale_intensity_method != "None":
      IJ.run(deconvolved, "16-bit", "");
    #deconvolved.show()
    
    #deconvolved = ops.run("convert.uint16", deconvolved)
    if(channels_luts[(c-1)]):
      deconvolved.setLut(channels_luts[(c-1)])
    else:
      deconvolved.setLut(imp_img_c.getLuts()[0])
    deconvolution_result_list.append(deconvolved)

  colorImp = RGBStackMerge.mergeChannels(deconvolution_result_list, True) #set to false to remove the original imageplus
  IJ.run(colorImp, "Properties...", "channels="+str(n_channel)+" slices="+str(n_z)+" frames="+str(n_f)+" unit=um pixel_width="+str(physSizeX.value(UNITS.MICROM))+" pixel_height="+str(physSizeY.value(UNITS.MICROM))+" voxel_depth="+str(physSizeZ.value(UNITS.MICROM)) )
  
  IJ.save(colorImp, os.path.join(output_folder.getAbsolutePath(), img_name[0:-4] + "_deconvolved.tiff"))
  IJ.log(img_name+' processed');

if(use_in_script == True):
  print('Deconvolution is Done!')
else:
  gd = GenericDialog('Deconvolution is Done!')
  gd.addMessage('Deconvolution is Done!')
  gd.showDialog()

'''
IJ.run("Close All", "")
'''
