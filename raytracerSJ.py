'''
Project A: An optical ray tracer
Aim: A module to model an optical system using the principle of geometric optics with different types of optical surfaces,
and subsequently optimising the curvature of a singlet lens. This module contains codes for plotting all the graphs.

Written by Son-Gyo Jung; Group F; Lab Group A9; 04.02.16
'''

import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *
import scipy.optimize as spo
import genpolarSJ as gp



#Normalisation function and Snell's law are placed at the beginning for convenience.
def normalise(x): 
    """
    Normalises a given 3D vector.
    
    Parameter:
        x - a 3D vector that is an array.
    """
    if x is not None and len(x) == 3:
        return x/np.sqrt(np.dot(x,x)) 

    else:
        raise Exception("Three coordinates are required in the array.")
        return None  


        
def Snell(k1, s_normal, n1, n2): #note: k1 must be a array
    '''
    Refracts the incident rays using Snell's law.
    
    Parameters:
        k1 - the direction of the incident ray.
        s_normal - vector normal to the surface.
        n1,n2 - the refractive indices either side of the surface.
    '''

    n_12 = n1/n2
    k1_hat = normalise(k1)
    s_normal_hat = normalise(s_normal)
    cos_angle1 = fabs(np.dot(k1_hat, s_normal_hat))
    angle1 = acos(cos_angle1)  
    

    if sin(angle1) >= n2/n1:
        print "Total internal reflection."

        return None     

    else:
        k2 = n_12*k1_hat + (n_12*cos_angle1 - np.sqrt(1 - n_12**2 * (1-cos_angle1)**2)) * s_normal_hat #vector form of Snell's law
        k2_hat = normalise(k2)

        return k2_hat
        


class ray:
    """
    Class defining the positions and directions of optical rays.
    
    Parameters:
        p - the position of the ray.
        k - the direction of the ray.
    """

    def __init__(self, p=[0.0, 0.0, 0.0], k=[0.0, 0.0, 1.0]): 
        self.__positions = list()
        self.__directions = list()

        
        if p is None and len(p) != 3:
            raise Exception("p requires exactly 3 arguments.")

        else:
            self.__positions.append(np.array(p))
            
        if k is None and len(k) != 3:
            raise Exception("k requires exactly 3 arguments.")

        else:
            self.__directions.append(normalise(np.array(k))) #append normalised k as an array
            

    def p(self):
        ''' 
        Returns the current position of the ray.
        '''
        return self.__positions[-1]

        
    def k(self):
        ''' 
        Returns the current direction of the ray.
        '''
        return self.__directions[-1] 

        
    def vertices(self):
        '''
        Returns the entire points through which the ray propagated.
        '''
        return self.__positions 

        
    def klist(self):
        '''
        Returns the list of directions of the ray.
        '''
        return self.__directions 

        
    def append(self, new_p=None, new_k=None):
        '''
        Add new position and direction of the ray
        
        Parameters:
            new_p - the new position of the ray.
            new_k - the new direction of the ray.
        '''

        if new_p is not None and len(new_p) == 3:
            self.__positions.append(np.array(new_p))

        else:
            raise Exception("P should only take 3 arguments.")
            
        if new_k is not None and len(new_k) == 3:
            self.__directions.append(np.array(new_k)) 

        else:
            raise Exception("K should only take 3 arguments.")
           


class OpticalElement: #provided in BlackBoard
    '''
    A base class for optical elements being employed; it propagates a ray through the optical element
    '''

    def propagate_ray(self, ray):

        raise NotImplementedError()        



class SphericalRefraction(OpticalElement):
    """
    Initialises the interface of refraction with either a spherical lens or a plane lens, where the latter corresponds to the curvature of zero.
    It returns the location of the intercepts and the new direction calculated by calling the Snell's law function.
    
    Parameters:
        z0 - the intercept of the surface with the z-axis.
        curv - the curvature of the surface.
        n1,n2 - the refractive indices either side of the surface.
        apr - aperture radius which is the maximum extent of the surface from the optical axis.
    """

    def __init__(self, z0=[0, 0, 100], curv=0.03, n1=1.0, n2=1.5168, apr=1000):
        self.__z0 = z0
        self.__curv = curv
        self.__n1 = n1
        self.__n2 = n2
        self.__apr = apr

        
    def intercept(self, ray):
        '''
        Determines the point of interception of the ray with the given lens.
        '''

        if self.__curv != 0: #when the lens is spherical  
            R = 1.0/self.__curv
            self.__centre = np.array([0., 0., R + self.__z0[2]])
            r = ray.p() - self.__centre
            discriminant = np.inner(r, normalise(ray.k()))**2 - (np.inner(r, r) - R**2) #using the equation provided in BlackBoard.
        

            if discriminant < 0:
                print "No real solution as the ray is complex."
                return None #No real solutions so terminate
    

            Lplus = -np.inner(r, normalise(ray.k())) + np.sqrt(discriminant) #Possible distance to the intercept
            Lminus = -np.inner(r, normalise(ray.k())) - np.sqrt(discriminant) #Possible distance to the intercept
            interceptplus = ray.p() + Lplus * normalise(ray.k()) #coordinate of the possible intercept
            interceptminus = ray.p() +Lminus * normalise(ray.k()) #coordinate of the possible intercept
                 

            if min(self.__z0[2],self.__centre[2]) <= interceptminus[2] <= max(self.__z0[2],self.__centre[2]):
                return interceptminus 
               

            if min(self.__z0[2],self.__centre[2]) <= interceptplus[2] <= max(self.__z0[2],self.__centre[2]):
                return interceptplus
                

        else: #when the lens is planar
            number_of_khatz = (self.__z0[2] - ray.p()[2])/normalise(ray.k())[2] #See lab book page 112
            intercept_in_plane = ray.p() + number_of_khatz * normalise(ray.k())
            return intercept_in_plane
                


    def propagate_ray(self,ray): 
        '''
        Determines the new position and direction of the ray after propagating through the lens.
        '''
        intersection = self.intercept(ray)
        if intersection is None:
            print "No intersections with the lens."
            return None

        else: #determining the direction of the surface normal vector (using the convention taught in the vector calculus course)
            if self.__curv == 0:
                new_s_normal_hat = np.array([0., 0., -1.])

            elif self.__curv > 0:
                new_s_normal_hat = normalise(intersection - self.__centre)

            elif self.__curv < 0:
                new_s_normal_hat = normalise(self.__centre - intersection)  

            else:
                raise Exception("Curvature is not defined.")         
                            
            newdir = Snell(normalise(ray.k()), new_s_normal_hat, self.__n1, self.__n2) #calling Snell's law
            newpos = intersection

            ray.append(newpos, newdir)
                      


class OutputPlane(SphericalRefraction): #take SphericalRefraction class as an argument to inherit the methods
    '''
    Terminates the propagation of the ray by creating a plane with a normal vector parallel to the z-axis.
    Note: the output plane is modelled using a planar lens that has infinite aperture radius.
    '''

    def __init__(self, output_z): 
        SphericalRefraction.__init__(self,[0,0,output_z], 0.00, 1., 1., float('inf')) #design output plane using a plane lens with infinite aperture radius 
            


class raybundle(list): 
    '''
    Generates a bundle of collimated rays using the genpolar module. 
    Note: the parameters can be adjusted to vary the beam radius, distribution and the number of rays.
    
    Parameters:
        n - Number of rings (this includes the central ray).
        rmax - Parameter which determines the maximum radius. (Note: the array created starts with zero radius by default.)
        m - The step-size of number of rays in each subsequent ring.
    '''

    def __init__(self, n=6, rmax=6., m=6,k=[0., 0., 1.]):
        for r, t in gp.rtuniform(n, rmax, m):
            self.append(ray([r*cos(t), r*sin(t), 0], normalise(k)))
        

    def propagate_ray(self, lens_element):
        '''
        Propagates the given rays through a single optical element at a time.
        '''
        for r in self:
            lens_element.propagate_ray(r)
            

    def RMS(self):
        '''
        Caculates the Root-Mean-Square of the spot radius at the paraxial focus for the bundle of rays; 
        that is the RMS deviation of the ray positions from the optical axis.
        '''

        rxy_squared = list()

        for r in self:
            rxy_squared.append(r.p()[0]**2 + r.p()[1]**2)   

        RMS = np.sqrt(np.mean(rxy_squared))*1e-3

        print 'RMS is %sm' % (RMS)

        return RMS


                                            
    def plot2D(self):
        '''
        Plots the paths of the rays in two-dimensions.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for r in self:
            x, y, z = zip(*r.vertices())
            ax.plot(z, x, 'b')
        ax.set_title("Ray propagation in z-x plane")
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("x (mm)")
        plt.grid()
        plt.show()
        


    def plot3D(self):
        '''
        Plots the paths of the rays in three-dimensions.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        for r in self:
            x, y, z = zip(*r.vertices())
            ax.plot(z, x, y) #use c = np.random.rand(3,1) to make it pretty
        ax.set_title("Ray propagation in 3D")
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_zlabel("x (mm)")
        plt.grid()
        plt.show()


    
    def plot2DOutput(self):
        '''
        Plots the rays that terminate in the output plane as a scatter graph in the x-y plane.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for r in self:
            x, y, z = zip(*r.vertices())
            ax.plot( x[-1],y[-1], 'bo')
        ax.set_title("Output in the x-y Plane")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
        plt.grid()
        plt.show()



        
#Extension tasks and additional calculators
def ParaxialFocus(fristlens, secondlens): #to calcuate the paraxial focus 
    '''
    Accurately estimates the location of the paraxial focus in the z-direction by utilising a paraxial ray.
    
    Parameters:
        fristlens - the front/frist optical surface defined using the SphericalRefraction class.
        secondlens - the back/second optical surface defined using the SphericalRefraction class.
    '''

    paraxial_ray = ray([0.00001, 0., 0.], [0., 0., 1.,]) #defines the paraxial ray
    fristlens.propagate_ray(paraxial_ray)

    if secondlens is not None:
        secondlens.propagate_ray(paraxial_ray)
    

    paraxial_focus = paraxial_ray.p()[2] - paraxial_ray.p()[0] * paraxial_ray.k()[2]/paraxial_ray.k()[0] # using eq. x = m(z-f), see lab book p.132
    
    print 'Paraxial focus is at %smm' % (paraxial_focus)

    return paraxial_focus
    


def DiffractionLimit(fristlens, secondlens, FirstLensPos, rmax, n, image_dis, optimisation): #to calculate the diffraction limit
    '''
    Accurately estimates the diffraction limit by employing the ParaxialFocus function.
    
    Parameters:
        fristlens, secondlens - first and second optical surface defined using the SphericalRefraction class.
        rmax, n - defined using parameters of raybundle class.      
    '''

    wavelength = 588e-9 
    diameter = np.arange (0.0, rmax, rmax/n)[-1] * 2e-3 
    
    if optimisation is True:
        f_length = image_dis*1e-3 - FirstLensPos*1e-3 - 2.5*1e-3 #for optimising biconvex lens with 5mm thickness

    elif  secondlens is not None: 
        f_length = ParaxialFocus(fristlens, secondlens)*1e-3 - FirstLensPos*1e-3 - 2.5*1e-3 #for a lens of 5mm thickness
    
    else:
        f_length = ParaxialFocus(fristlens, secondlens)*1e-3 - FirstLensPos*1e-3 #for a single refracting surface 
       

    D_limit = wavelength*f_length/diameter

    print 'Diffraction limit is %sm' % (D_limit) 

    return D_limit 
        


def OptimiseLenses(lensposition,image_dis, n, rmax, m, k): #to optimise a singlet lens into its best form for a given image distance.
    '''
    Optimises the curvatures of both lenses by minimising the RMS of spot radius at a given image distance;
    Note: both the positions of the lenses and the output plane is required.
    
    Parameters:
        lensposition - essentially the z-coordinate of z0.
        image_dis - the position of the image (i.e. the position of the output plane).
        n, rmax, m, k - defined using parameters of raybundle class.
    '''    

    def singlet(curvatures):
        '''
        Parameter:
            curvatures - an optimizer which minimises the RMS spot radius using the downhill simplex algorithm.
            [Args: singlet - the objective function to be minimized; trial_solution - the initial guesses.]              
        '''     

        first_lens = SphericalRefraction(z0=[0., 0., lensposition], curv=curvatures[0], n1=1.0, n2=1.5168, apr=float('inf'))
        second_lens = SphericalRefraction(z0=[0., 0., lensposition+5], curv=curvatures[1], n1=1.5168, n2=1.0, apr=float('inf'))
        opp = OutputPlane(image_dis) #output plane positioned at the give image distance
        
        b = raybundle(n, rmax, m, k) #substitute with new values using ExecuteRayTracerSJ.py
        b.propagate_ray(first_lens)
        b.propagate_ray(second_lens)
        b.propagate_ray(opp)
        

        rxy_squared = list()

        for r in b:
            rxy_squared.append(r.p()[0]**2 + r.p()[1]**2)
        RMS = np.sqrt(np.mean(rxy_squared))*1e-3
        print 'RMS is %sm' % (RMS)

        return RMS

        
    trial_solutions = [0.05, -0.05]
    curvatures = spo.fmin(singlet, trial_solutions)

    print 'The best form curvatures for a given image distance of %s mm are: first lens curvature = %s/mm, second lens curvature = %s/mm.' % (image_dis, curvatures[0], curvatures[1])
                  


def RMSvsDiffractionLimit(lensposition, image_dis, n, rmax, rmin, steps, m, k, optcurv_1, optcurv_2): #function to plot RMS spot radius and diffraction limit against increasing beam radius
    rmaxlist = np.linspace(rmin, rmax, steps)
    RMSlist = list() 
    DLlist = list()
   
    first_lens = SphericalRefraction(z0=[0., 0., lensposition], curv=optcurv_1, n1=1.0, n2=1.5168, apr=float('inf'))
    second_lens = SphericalRefraction(z0=[0., 0., lensposition+5], curv=optcurv_2, n1=1.5168, n2=1.0, apr=float('inf'))
    opp = OutputPlane(image_dis) 

    for rmax in rmaxlist:
        b = raybundle(n, rmax, m, k) #substitute with new values using ExecuteRayTracerSJ.py
        b.propagate_ray(first_lens)
        b.propagate_ray(second_lens)
        b.propagate_ray(opp)
        RMSlist.append(b.RMS()) #append all the RMS values
        
        DLlist.append(DiffractionLimit(first_lens, second_lens, lensposition, rmax, n, image_dis, optimisation=True)) #append the diffraction limits calculated using the DiffractionLimit function
        

    plt.figure()
    line_up, = plt.plot(rmaxlist,RMSlist,'g--', label='RMS Spot Radius')
    line_down, = plt.plot(rmaxlist,DLlist,'r-', label='Diffraction Limit')

    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.grid()
    plt.xlabel("Beam radius (mm)")
    plt.ylabel('Limiting factor (m)')
    plt.title('Diffraction Limit and RMS Spot Size with Beam Radius \n for Optimised Biconvex lens', fontsize=13)
    plt.autoscale(tight=True) 
    plt.legend(handles=[line_up, line_down])  
    plt.show()
    


def RMSPlanoConvex(lensposition, n, rmax, rmin, steps, m, k, optcurv_1, optcurv_2): #to plot the behaviour of RMS spot radius with increasing beam size when using the planocovex lens
    rmaxlist = np.linspace(rmin, rmax, steps)
    RMSlist = list()
    RMSlist2 = list() #both orientations considered
   
    first_lens = SphericalRefraction(z0=[0., 0., lensposition], curv=optcurv_1, n1=1.0, n2=1.5168, apr=float('inf'))
    second_lens = SphericalRefraction(z0=[0., 0., lensposition+5], curv=optcurv_2, n1=1.5168, n2=1.0, apr=float('inf'))
    third_lens = SphericalRefraction(z0=[0., 0., lensposition], curv=optcurv_2, n1=1.0, n2=1.5168, apr=float('inf'))
    fourth_lens = SphericalRefraction(z0=[0., 0., lensposition+5], curv=-optcurv_1, n1=1.5168, n2=1.0, apr=float('inf'))
    opp1 = OutputPlane(ParaxialFocus(first_lens,  second_lens)) 
    opp2 = OutputPlane(ParaxialFocus(third_lens, fourth_lens))
    

    for rmax in rmaxlist:
        b = raybundle(n, rmax, m, k) #substitute with new values using ExecuteRayTracerSJ.py
        b.propagate_ray(first_lens)
        b.propagate_ray(second_lens)
        b.propagate_ray(opp1)
        RMSlist.append(b.RMS()) #list of all the RMS values
        

    for rmax in rmaxlist:
        b2 = raybundle(n, rmax, m, k) #substitute with new values using ExecuteRayTracerSJ.py
        b2.propagate_ray(third_lens)
        b2.propagate_ray(fourth_lens)
        b2.propagate_ray(opp2)
        RMSlist2.append(b2.RMS()) #list of all the RMS values
          

    plt.figure()
    line_up, = plt.plot(rmaxlist,RMSlist, 'g-', label='Convex-Plano')
    line_down, = plt.plot(rmaxlist,RMSlist2, 'b--', label='Plano-Convex')
    plt.legend(handles=[line_up, line_down], loc=2)  
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.grid()
    plt.xlabel("Beam radius (mm)")
    plt.ylabel('RMS Spot Size (m)')
    plt.title('RMS spot size with beam radius\n(planoconvex lenses)')
    plt.autoscale(tight=True) 
    plt.show()
    
