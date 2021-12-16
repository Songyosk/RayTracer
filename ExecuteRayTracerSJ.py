'''
Project A: An optical ray tracer
Aim: To run the python modules and generate results. Please only change the Model name.

Written by Son-Gyo Jung; Group F; Lab Group A9; 04.02.16
'''
import raytracerSJ as rt
reload(rt)

'Insert the name of the model to be considered below:'
#Write one of the following string: 'singlesurface', 'plano-convex', 'convex-plano' or 'optimised'.
Model = 'singlesurface'


####################### Don't touch anything below this line ########################
if Model is 'singlesurface':
    'Parameters for lenses and output plane'
    FirstLensPos = 100
    curv_1 = 0.03
    curv_2 = None
    n1 = 1.0
    n2 = 1.5 #or 1.5168
    apr = 1000
    pos_outputplane = image_dis = 200 #or 197.832817337 when using 1.5168 refractive index
    optimisation = False #True or False
    

elif Model is 'convex-plano':
    'Parameters for lenses and output plane'
    FirstLensPos = 100
    curv_1 = 0.02 #or 0.
    curv_2 = 0. #or -0.02
    n1 = 1.0
    n2 = 1.5168
    apr = 1000.
    pos_outputplane = image_dis = 198.452812504 #or 201.749226006 when curv_2 is -0.02.
    optimisation = False #True or False


elif Model is 'plano-convex':
    FirstLensPos = 100
    curv_1 = 0.
    curv_2 = -0.02
    n1 = 1.0
    n2 = 1.5168
    apr = 1000.
    pos_outputplane = image_dis = 201.749226006 
    optimisation = False #True or False  
    

elif Model is 'optimised':
    'Parameters for lenses and output plane'
    FirstLensPos = 100
    curv_1 = 0.0141749158946 #for optimisation the curvatures can be any random number.
    curv_2 = -0.00532044505389
    n1 = 1.0
    n2 = 1.5168
    apr = 1000.
    pos_outputplane = image_dis = 202.5 
    optimisation = True #True or False



'Parameters for ray bundle'
n = 6
rmax = 5.+5./5. #for a max of x/mm, rmax+rmax/n
m = 6
k = [0.,0.,1.] 


'Implement modules'
r = rt.ray()
s1 = rt.SphericalRefraction([0,0,FirstLensPos], curv_1, n1, n2, apr)
if curv_2 is not None:
    s2 = rt.SphericalRefraction([0,0,FirstLensPos+5], curv_2, n2, n1, apr)
s3 = rt.OutputPlane(pos_outputplane)


s1.propagate_ray(r)
if curv_2 is not None:
    s2.propagate_ray(r)


b = rt.raybundle(n, rmax, m, k) #angle can be given
b.propagate_ray(s1)
if curv_2 is not None:
    b.propagate_ray(s2)

b.propagate_ray(s3)

b.plot2DOutput()
b.plot2D()
b.plot3D()


'Extensions tasks and calculators'  
if optimisation is True:
    rt.RMSvsDiffractionLimit(FirstLensPos, image_dis, n, 20, 1, 30, m, k, 0.0141749158946, -0.00532044505389) #rmax=20, rmin=1, divisions=30 for the graph of RMS and diffraction limit against beam radius
    rt.OptimiseLenses(FirstLensPos, image_dis, n, rmax, m, k)

elif curv_2 is not None and Model is 'convex-plano':
    rt.RMSPlanoConvex(FirstLensPos, n, 10, 1, 30, m, k, curv_1, curv_2)

elif curv_2 is not None and Model is not 'convex-plano':
    rt.RMSPlanoConvex(FirstLensPos, n, 10, 1, 30, m, k, curv_2, curv_1)
    
if curv_2 is not None:
    rt.DiffractionLimit(s1, s2, FirstLensPos, rmax, n, image_dis, optimisation)
    
    if optimisation is False:
        rt.ParaxialFocus(s1, s2)        

else:
    rt.DiffractionLimit(s1, None, FirstLensPos, rmax, n, image_dis, optimisation)
    rt.ParaxialFocus(s1, None) 
      
if Model is not 'optimised':
    b.RMS()            

