"""Python function to plot SubCortexMesh subcortices
"""
import numpy as np
import copy

from brainspace.plotting import plot_surf as ps
from brainspace.mesh import mesh_creation as mc
import vtk
import numpy as np

def surfplot_subcortical(cdata, lh_celldata, rh_celldata, hemis=['L','R'],size=[350,400], smoothing=0, **qwargs):
    '''
    cdata: array with the shape Vx2xF, where V is the number of vertices, 2 is the number of hemispheres (unless specified), and F is the number of rows/features
    lh_celldata: list containing point and cell data for the left hemisphere, loaded from the R plot_surf() function.
    rh_celldata: list containing point and cell data for the right hemisphere, loaded from the R plot_surf() function.
    kwargs: see https://brainspace.readthedocs.io/en/latest/generated/brainspace.plotting.surface_plotting.plot_surf.html#brainspace.plotting.surface_plotting.plot_surf
    '''
    
    # build left and right surface
    lh = mc.build_polydata(lh_celldata["vertices"], cells=lh_celldata["triangles"])
    rh = mc.build_polydata(rh_celldata["vertices"], cells=rh_celldata["triangles"])
    
    #smoothing option
    if smoothing > 0:
      def smoother(mesh, smoothing):
        s = vtk.vtkWindowedSincPolyDataFilter()
        s.SetInputData(mesh.VTKObject)
        s.SetNumberOfIterations(int(smoothing))
        s.SetPassBand(0.001)
        s.NonManifoldSmoothingOn()
        s.NormalizeCoordinatesOn()
        s.Update()
        mesh.VTKObject.ShallowCopy(s.GetOutput())
        return mesh
    
      lh=smoother(lh, smoothing)
      rh=smoother(rh, smoothing)
    
    #cdata formatting for 2D and 3D input
    cdata = np.reshape(cdata,[cdata.shape[0],len(hemis),-1])
    if len(cdata.shape) == 2: cdata = np.expand_dims(cdata,axis=2)

    #ROI-specific optimised angles/axes
    lh,rh=coord_optimizer(cdata,lh,rh)
    #flipped meshes (will be reangled as dorsal)
    lf = mc.build_polydata(lh.Points.copy(), cells=lh.GetCells2D().copy())
    rf = mc.build_polydata(rh.Points.copy(), cells=rh.GetCells2D().copy())
    
    #keep N vertices across hemis as they vary
    N_lh = lh.points.shape[0] 
    N_rh = rh.points.shape[0]
    
    # set up layout
    surfDict = {'Lh':lh,  'Lf':lf, 'Rh':rh, 'Rf':rf}
    surfList = np.ones((cdata.shape[2],len(hemis)*2),dtype=object)
    arrName = np.ones((cdata.shape[2],len(hemis)*2),dtype=object)
    for h,hemi in enumerate(hemis):
        if hemi=='L':
            surfList[:,[h,h+1]] = np.array([f"{hemi}h",f"{hemi}f"])
            for f in range(cdata.shape[2]):
              lh.append_array(cdata[:N_lh,h,f], name=f'feature{f}', at='point')
              lf.append_array(cdata[:N_lh,h,f], name=f'feature{f}', at='point')
        elif hemi=='R':
            surfList[:,[h*2,h*2+1]] = np.array([f"{hemi}f",f"{hemi}h"])
            for f in range(cdata.shape[2]): 
              rh.append_array(cdata[:N_rh,h,f], name=f'feature{f}', at='point')
              rf.append_array(cdata[:N_rh,h,f], name=f'feature{f}', at='point')
        for f in range(cdata.shape[2]):
            arrName[f,:] = f'feature{f}'  
        
    # extra parameters
    new_qwargs = dict(zoom=1.7, nan_color=(0,0,0,0))
    new_qwargs.update(qwargs)
    new_size=copy.deepcopy(size)
    new_size[0] = new_size[0]*len(hemis)
    new_size[1] = new_size[1]*cdata.shape[2]
    if 'color_bar' in qwargs:
        new_size[0] = new_size[0]+60
        
    #view to slightly differ for Brain-Stem (as no hemispheres)
    #and cerebellum (as inside/coronal view may be less relevant)
    if cdata.shape[0] == 9452 or cdata.shape[0] == 95718:
      view=['ventral','dorsal','medial','lateral']
    elif cdata.shape[0] == 19559:
      view=['ventral','medial','lateral','ventral']
      lf.Points = rotate_points(lf.Points, axis='x', angle_deg=90)
      rf.Points = rotate_points(rf.Points, axis='x', angle_deg=90)
    else: 
      view=['ventral','dorsal','dorsal','ventral']
    
    # plot
    p = ps(surfDict,surfList, array_name=arrName, size=new_size,view=view, **new_qwargs)
    return p

#rotation function to adapt meshes visually
def rotate_points(points, axis='z', angle_deg=90):
    angle = np.deg2rad(angle_deg)
    if axis == 'x':
        R = np.array([[1, 0, 0],
                      [0, np.cos(angle), -np.sin(angle)],
                      [0, np.sin(angle),  np.cos(angle)]])
    elif axis == 'y':
        R = np.array([[ np.cos(angle), 0, np.sin(angle)],
                      [0, 1, 0],
                      [-np.sin(angle), 0, np.cos(angle)]])
    elif axis == 'z':
        R = np.array([[np.cos(angle), -np.sin(angle), 0],
                      [np.sin(angle),  np.cos(angle), 0],
                      [0, 0, 1]])
    else:
        raise ValueError("axis must be 'x', 'y', or 'z'")
    
    return points @ R.T

#ROI-specific visual parameters
#/!\ cdata.shape[0] will have the vertex count of the largest hemisphere, so not always the left
def coord_optimizer(cdata, lh, rh):
  #accumbens nuclei
  if cdata.shape[0] == 1022:
    lh.Points = rotate_points(lh.Points, axis='y', angle_deg=-30)
    lh.Points = rotate_points(lh.Points, axis='x', angle_deg=-25)
    rh.Points[:,2]=rh.Points[:,2]/1.7 #Z changed due to asymmetry
    rh.Points = rotate_points(rh.Points, axis='y', angle_deg=60)
  #amygdalae
  if cdata.shape[0] == 1792:
    lh.Points = rotate_points(lh.Points, axis='z', angle_deg=-10)
    rh.Points = rotate_points(rh.Points, axis='x', angle_deg=30)
  #caudate
  if cdata.shape[0] == 3500:
    lh.Points = rotate_points(lh.Points, axis='y', angle_deg=-30)
    rh.Points = rotate_points(rh.Points, axis='y', angle_deg=30)
    lh.Points = rotate_points(lh.Points, axis='x', angle_deg=10)
    rh.Points = rotate_points(rh.Points, axis='x', angle_deg=10)
  #cerebellum
  if cdata.shape[0] == 19664:
    lh.Points = rotate_points(lh.Points, axis='y', angle_deg=0)
    rh.Points = rotate_points(rh.Points, axis='y', angle_deg=0)
    lh.Points = rotate_points(lh.Points, axis='x', angle_deg=180)
    rh.Points = rotate_points(rh.Points, axis='x', angle_deg=180)
  #hippocampus
  if cdata.shape[0] == 4086:
    lh.Points = rotate_points(lh.Points, axis='x', angle_deg=35)
    rh.Points = rotate_points(rh.Points, axis='x', angle_deg=35)
  #pallidum 
  if cdata.shape[0] == 1600:
    lh.Points = rotate_points(lh.Points, axis='y', angle_deg=30)
    rh.Points = rotate_points(rh.Points, axis='y', angle_deg=-30)
    lh.Points = rotate_points(lh.Points, axis='z', angle_deg=90)
    rh.Points = rotate_points(rh.Points, axis='z', angle_deg=-90)
  #putamen 
  if cdata.shape[0] == 4268:
    lh.Points = rotate_points(lh.Points, axis='y', angle_deg=60)
    rh.Points = rotate_points(rh.Points, axis='y', angle_deg=-60)
    lh.Points = rotate_points(lh.Points, axis='z', angle_deg=90)
    rh.Points = rotate_points(rh.Points, axis='z', angle_deg=-90)
  #Thalamus 
  if cdata.shape[0] == 3936:
    lh.Points = rotate_points(lh.Points, axis='y', angle_deg=60)
    rh.Points = rotate_points(rh.Points, axis='y', angle_deg=-60)
    lh.Points = rotate_points(lh.Points, axis='z', angle_deg=-90)
    rh.Points = rotate_points(rh.Points, axis='z', angle_deg=90)
  #Ventral diencephalon 
  if cdata.shape[0] == 3594:
    lh.Points = rotate_points(lh.Points, axis='x', angle_deg=60)
    rh.Points = rotate_points(rh.Points, axis='x', angle_deg=60)
    lh.Points = rotate_points(lh.Points, axis='z', angle_deg=-50)
    rh.Points = rotate_points(rh.Points, axis='z', angle_deg=50)
  #Brain-Stem 
  if cdata.shape[0] == 9452:
    lh.Points = rotate_points(lh.Points, axis='z', angle_deg=-180)
    rh.Points = rotate_points(rh.Points, axis='x', angle_deg=-100)
  #All aseg
  if cdata.shape[0] == 95718:
    lh.Points = rotate_points(lh.Points, axis='x', angle_deg=-30)
    lh.Points = rotate_points(lh.Points, axis='z', angle_deg=180)
    rh.Points = rotate_points(rh.Points, axis='x', angle_deg=-90)
  
  return lh, rh
