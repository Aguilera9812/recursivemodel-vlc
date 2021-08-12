"""
RECURSIVE MODEL CHANNEL FOR VISIBLE LIGHT COMMUNICATION
Juan Felipe Gutierrez
jufgutierrezgo@unal.edu.co

This software includes the following improvements:
- Using of the fast euclidean distance function
- Add a new dimension to the array_points, a wall label
- The array_parameter is computed only with half matrix 
- Was created a general reports about channel impulse reponse.

"""
import numpy as np
import numpy
import math
import os


import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

import fractions
from fractions import Fraction

from fastdist import fastdist

import timeit

from numpy.core.function_base import linspace


c = 3e8 #speed of light in [m/s]
no_walls = 6
tres = 0.2e-9 # time resolution
bins = 300 # bins for power graph 

#Array with normal vectors for each wall.
ew_n = [[0,0,-1],[0,1,0],[1,0,0],[0,-1,0],[-1,0,0],[0,0,1]]

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

cir_path = ROOT_DIR + "/cir/"
report_path = ROOT_DIR + "/report/"


#Function to calculate angls between two vector position
def cos_2points(v1,n1,v2,n2):
    unit_vlos = (v1-v2) / np.linalg.norm(v1-v2)

    cos_phi = np.dot(-1*unit_vlos, n1)
    cos_tetha = np.dot(unit_vlos, n2)   
    
    #print([angle1,angle2])    
    return cos_phi,cos_tetha

def led_pattern(m):
    theta, phi = np.linspace(0, 2 * np.pi, 40), np.linspace(0,np.pi/2, 40)
    THETA, PHI = np.meshgrid(theta, phi)
    R = (m+1)/(2*np.pi)*np.cos(PHI)**m
    X = R * np.sin(PHI) * np.cos(THETA)
    Y = R * np.sin(PHI) * np.sin(THETA)
    Z = R * np.cos(PHI)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    plot = ax.plot_surface(
        X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('jet'),
        linewidth=0, antialiased=False, alpha=0.5)

    plt.show()
    return 0

#Function to calculate the points in each wall from room size specifications.
#x_lim: lenght in x axe.2
#y_lim: lenght in y axe.
#z_lim: lenght in z axe.
# returns 2d-array (3xNc) with [X,Y,Z] coordinates of each points.
def tessellation(x_lim,y_lim,z_lim,scale_factor):
    print("//****** Tessellation *******//")
    x_num = fractions.Fraction(str(x_lim)).numerator
    x_den = fractions.Fraction(str(x_lim)).denominator
    y_num = fractions.Fraction(str(y_lim)).numerator
    y_den = fractions.Fraction(str(y_lim)).denominator    
    z_num = fractions.Fraction(str(z_lim)).numerator
    z_den = fractions.Fraction(str(z_lim)).denominator    

    #print(x_num,x_den,y_num,y_den,z_num,z_den)
    den_lcm = lcm(x_den,y_den,z_den)
        
    n_x = int(den_lcm*x_num/x_den)
    n_y = int(den_lcm*y_num/y_den)
    n_z = int(den_lcm*z_num/z_den)
    #print(n_x,n_y,n_z)
    #print(den_lcm)

    num_gfc = math.gcd(n_z,math.gcd(n_x,n_y))
    #print(num_gfc)

    delta_Lmax = num_gfc / den_lcm 
    print("DeltaL max is: ", delta_Lmax)
    delta_Amax = delta_Lmax**2
    print("DeltaA max is: ", delta_Amax)
    
    #Scaling factor for delta_Lmax (1/2,1/3,1/4....)
    delta_L = delta_Lmax*scale_factor

    #DeltaA is defined from root(2)*DeltaL/2, because is the minimum distance between two points
    #delta_A = (delta_L/20)**2
    
    #DeltaA is defined from time-clearesolution presented in main reference.
    #delta_A = 3.6e-3
    
    #DeltaA defined to fulfill deltaA << (root(2)/2)*deltaL --> the maximun lenght on delta_A must be 10 times less than distane between points 
    delta_A = (delta_L**2)/800

    print("Scale factor for Delta L is: ", scale_factor)
    print("DeltaL[m]: ", delta_L)
    print("DeltaA[m^2]: ", delta_A)


    no_xtick = int(x_lim/delta_L)
    no_ytick = int(y_lim/delta_L)
    no_ztick = int(z_lim/delta_L)

    ew0_points = np.zeros((4,no_xtick*no_ytick))
    ew1_points = np.zeros((4,no_ztick*no_xtick))
    ew2_points = np.zeros((4,no_ztick*no_ytick))
    ew3_points = np.zeros((4,no_ztick*no_xtick))
    ew4_points = np.zeros((4,no_ztick*no_ytick))
    ew5_points = np.zeros((4,no_xtick*no_ytick))

    #Init_index define the index where each point in parameters start, eg. ew_0 start in index 0, and ew_1 in no_xtick*no_ytick index
    init_index = np.zeros(6)
    points = [ew0_points,ew1_points,ew2_points,ew3_points,ew4_points,ew5_points]

    

    for i in range(1,6):
        init_index[i] = int(len(points[i-1][0,:]) + init_index[i-1])
    
    counter_cell = 0

    for j in range(0,no_xtick):
        for i in range(0,no_ytick):
            ew0_points[0,counter_cell] = delta_L/2 + j*delta_L
            ew0_points[1,counter_cell] = delta_L/2 + i*delta_L         
            ew0_points[2,counter_cell] = z_lim
            ew0_points[3,counter_cell] = 0         

            ew5_points[0,counter_cell] = delta_L/2 + j*delta_L
            ew5_points[1,counter_cell] = delta_L/2 + i*delta_L         
            ew5_points[2,counter_cell] = 0
            ew5_points[3,counter_cell] = 5         
            counter_cell += 1

    counter_cell = 0

    for j in range(0,no_ztick):
        for i in range(0,no_xtick):
            ew1_points[0,counter_cell] = x_lim - delta_L/2 - i*delta_L
            ew1_points[1,counter_cell] = 0
            ew1_points[2,counter_cell] = z_lim - delta_L/2 - j*delta_L
            ew1_points[3,counter_cell] = 1        
            
            ew3_points[0,counter_cell] = x_lim - delta_L/2 - i*delta_L
            ew3_points[1,counter_cell] = y_lim
            ew3_points[2,counter_cell] = z_lim - delta_L/2 - j*delta_L
            ew3_points[3,counter_cell] = 3         
            counter_cell += 1

    
    counter_cell = 0

    for j in range(0,no_ztick):
        for i in range(0,no_ytick):
            ew2_points[0,counter_cell] = 0
            ew2_points[1,counter_cell] = delta_L/2 + i*delta_L/2
            ew2_points[2,counter_cell] = z_lim - delta_L/2 - j*delta_L
            ew2_points[3,counter_cell] = 2        
            
            ew4_points[0,counter_cell] = x_lim
            ew4_points[1,counter_cell] = delta_L/2 + i*delta_L/2
            ew4_points[2,counter_cell] = z_lim - delta_L/2 - j*delta_L
            ew4_points[3,counter_cell] = 4        
            counter_cell += 1

    no_points=2*no_xtick*no_ytick + 2*no_ztick*no_xtick + 2*no_ztick*no_ytick
    print("The total number of points is: ",no_points)
    print("//-------- points array created --------------//")
    #print(ew0_points)
    #print(ew5_points)
    return [np.concatenate((ew0_points,ew1_points,ew2_points,ew3_points,ew4_points,ew5_points),axis=1),no_xtick,no_ytick,no_ztick,init_index,delta_A,no_points]

def lcm(a, b,c):
    return abs(a*b*c) // math.gcd(c,math.gcd(a, b))

#Function to create a cross parameters between points.
def make_parameters(array_points,x_lim,y_lim,z_lim,no_xtick,no_ytick,no_ztick):

    no_points = 2*no_xtick*no_ytick + 2*no_ztick*no_xtick + 2*no_ztick*no_ytick
    ew_par = np.zeros((2,no_points,no_points),dtype=np.float16)    

    counter_points = 0

    for ini_point in range(0,no_points):        

        for end_point in range(ini_point+1,no_points):
                               
            if array_points[3,ini_point]==array_points[3,end_point]:
                ew_par[0,ini_point,end_point] = 0
                ew_par[1,ini_point,end_point] = 0
            else:
                #ew_par[0,ini_point,end_point] = math.dist(array_points[:,ini_point],array_points[:,end_point])
                wallinit = int(array_points[3,ini_point])               
                wallend = int(array_points[3,end_point])

                ew_par[0,ini_point,end_point] = fastdist.euclidean(array_points[0:3,ini_point],array_points[0:3,end_point])                 
                ew_par[0,end_point,ini_point]  = ew_par[0,ini_point,end_point]

                ew_par[1,ini_point,end_point],ew_par[1,end_point,ini_point] = cos_2points(array_points[0:3,ini_point],ew_n[wallinit],
                array_points[0:3,end_point],ew_n[wallend])
                
   

    print("//------- parameters array created -----------//")
    #print(h_k[i])   
    #numpy.savetxt("ew_par_dis.csv", ew_par[0,:,:], delimiter=",")  
    #numpy.savetxt("ew_par_cos.csv", ew_par[1,:,:], delimiter=",")  

    return ew_par

#Function to compute the channel impulse respone
# Inputs arguments:
# m: lambertian number to tx emission
# tx_pos: 1d-array with [x,y,z] tx position
# rx_pos: 1d-array with [x,y,z] rx position
# points: List with [x,y,z] cooridinates for every point in each wall
# parameters: List with angle and distance between all points.  
# x_lim,y_lim,z_lim: limits in room dimmensions
# a_r: sensitive area in photodetector
# no_xtick,no_ytick,no_ztick: number of division in each axes.
def h_t(m,tx_pos,rx_pos,points,wall_label,parameters,x_lim,y_lim,z_lim,no_xtick,no_ytick,no_ztick,init_index,a_r,rho,delta_A,k_reflec):
    
    #tx_wall_power = np.zeros(3,2*no_xtick*no_ytick + 2*no_ztick*no_xtick + 2*no_ztick*no_ytick)
    #rx_wall_power = np.zeros(3,2*no_xtick*no_ytick + 2*no_ztick*no_xtick + 2*no_ztick*no_ytick)

    no_cells = len(points[0,:])

    area_factor = (2*x_lim*y_lim + 2*x_lim*z_lim + 2*y_lim*z_lim)/(delta_A*no_cells)

    #define the wall of the tx_pos
    tx_wall = wall_label[tx_pos]
    

    #define the wall of the rx_pos
    rx_wall = wall_label[rx_pos]
    

    for i in range(0,no_cells):
        #print(np.transpose(tx_pos)-points[tx_wall][:,i])        
        if np.allclose(np.transpose(tx_pos),points[:,i]):
            tx_index_point = i
            #print(i)
            break
            

    for i in range(0,no_cells):        
        if np.allclose(np.transpose(rx_pos),points[:,i]):
            rx_index_point = i
            #print(i)
            break

    
    cos_phi = np.zeros((no_cells),dtype=np.float16)
    dis2 = np.zeros((no_cells,no_cells),dtype=np.float16)

    dis2 = np.power(parameters[0,:,:],2)
    
    cos_phi = parameters[1,int(tx_index_point),:]
    tx_power = (m+1)/(2*np.pi)*np.multiply(np.divide(1,dis2[tx_index_point,:],out=np.zeros((no_cells)), where=dis2[tx_index_point,:]!=0),np.power(cos_phi,m))
    rx_wall_factor = a_r*parameters[1,int(rx_index_point),:]

    h0_se = np.zeros((no_cells,2),dtype=np.float32)
    h0_er = np.zeros((no_cells,2),dtype=np.float32)

    
    #Impulse response between source and each discretized wall
    h0_se[:,0] = np.multiply(area_factor*rho*delta_A*tx_power,parameters[1,:,int(tx_index_point)])
    h0_er[:,0] = np.divide(rx_wall_factor,dis2[rx_index_point,:],out=np.zeros((no_cells)), where=dis2[rx_index_point,:]!=0)
    h0_se[:,1] = parameters[0,tx_index_point,:]/c
    h0_er[:,1] = parameters[0,rx_index_point,:]/c

    dP_ij = np.zeros((no_cells,no_cells),np.float32)
    dP_ij = np.divide(rho*delta_A*parameters[1,:,:],dis2,out=np.zeros_like(dP_ij),where=dis2!=0) 
    #dP_ij_1d = dP_ij.flatten()
    #numpy.savetxt("dPij.csv", dP_ij[:,0], delimiter=",")
    

    h_k = []
    hlast_er = []
    

    for i in range(k_reflec+1):
        h_k.append(np.zeros((int(no_cells**i),2),np.float32))
        hlast_er.append(np.zeros((int(no_cells**i),2),np.float32)) 

        if i == 0:           

            h_k[i][0,0] = tx_power[int(rx_index_point)]*rx_wall_factor[int(tx_index_point)]
            h_k[i][0,1] = parameters[0,int(tx_index_point),int(rx_index_point)]/c

            print("//------------- h0-computed ------------------//")            
            numpy.savetxt(cir_path+"h0.csv", h_k[i], delimiter=",")

        elif i==1:
            hlast_er[i][:,0] = h0_er[:,0]     
            hlast_er[i][:,1] = h0_er[:,1]

            h_k[i][:,0] = np.multiply(h0_se[:,0],h0_er[:,0]) 
            h_k[i][:,1] = h0_se[:,1] + h0_er[:,1]

            print("//------------- h1-computed ------------------//")
            numpy.savetxt(cir_path+"h1.csv", h_k[i], delimiter=",")
            

        else:
            count_blocks = 0
            #print(len(hlast_er[i-1][:,0]))
            #print(len(hlast_er[i][:,0]))
            for j in range(len(hlast_er[i-1][:,0])):

                index_dpij = int(j%no_cells)
                
                hlast_er[i][no_cells*j:int(no_cells*(j+1)),0] = hlast_er[i-1][j,0]*dP_ij[index_dpij,:]
                hlast_er[i][no_cells*j:int(no_cells*(j+1)),1] = hlast_er[i-1][j,1] + parameters[0,index_dpij,:]/c                

            len_last = len(hlast_er[i][:,0])

            for l in range(no_cells):
                
                lim_0 = int(l*(no_cells**(i-1)))
                lim_1 = int((l+1)*(no_cells**(i-1)))
                
                #h_k[i][lim_0:lim_1,0] = h0_se[l,0]*[hlast_er[i][m,0] for m in range(l,len_last,no_cells)]
                #h_k[i][lim_0:lim_1,1] = h0_se[l,1] + [hlast_er[i][m,1] for m in range(l,len_last,no_cells)]

                #h_k[i][lim_0:lim_1,0] = h0_se[l,0]*hlast_er[i][lim_0:lim_1,0] 
                #h_k[i][lim_0:lim_1,1] = h0_se[l,1] + hlast_er[i][lim_0:lim_1,1]
                #print(h0_se[l,0])
                #print([hlast_er[i][m,0] for m in range(l,len_last,no_cells)])
                h_k[i][lim_0:lim_1,0] = np.multiply([hlast_er[i][m,0] for m in range(l,len_last,no_cells)],h0_se[l,0])
                h_k[i][lim_0:lim_1,1] = h0_se[l,1] + [hlast_er[i][m,1] for m in range(l,len_last,no_cells)]

            print("//------------- h"+str(i)+"-computed ------------------//")            
            numpy.savetxt(cir_path+"h"+str(i)+".csv", h_k[i], delimiter=",")      
    
      
    return h_k



#Function to create an analysis of the simulation
def create_report(h_k,k_reflec,no_cells):

    print("//------------- Data report ------------------//")
    print("Time resolution [s]:"+str(tres))
    print("Number of Bins:"+str(bins))
    h_power = np.zeros((k_reflec+1))

    hk_aux = []
    
    delay_los = h_k[0][0,1]
    power_data = np.zeros((bins,k_reflec+1))


    for i in range(k_reflec+1):            
        
        hk_aux.append(np.zeros((int(no_cells**i),2)))

        # Compute and print the total power per order reflection
        print("h"+str(i)+"-Response:")                       
        h_power[i] = np.sum(h_k[i][:,0])
        print("Power[w]:",h_power[i])
        if i==0:
            print("Delay[s]:",h_k[i][0,1])

        # Create graphs 
               
        hk_aux[i] = h_k[i]
        hk_aux[i][:,1] = hk_aux[i][:,1] - delay_los
        hk_aux[i][:,1] = np.floor(hk_aux[i][:,1]/tres)

        for j in range(no_cells**i):
           power_data[int(hk_aux[i][j,1]),i] += hk_aux[i][j,0]

                
        time_scale = linspace(0,bins*tres,num=bins)

        fig, (vax) = plt.subplots(1, 1, figsize=(12, 6))
        vax.plot(time_scale,power_data[:,i], 'o',markersize=2)
        vax.vlines(time_scale, [0], power_data[:,i],linewidth=1)

        vax.set_xlabel("time(s) \n Time resolution:"+str(tres)+"s  Bins:"+str(bins),fontsize=15)
        vax.set_ylabel('Power(W)',fontsize=15)
        vax.set_title("Channel Impulse Response h"+str(i)+"(t)",fontsize=20)

        vax.grid(color = 'black', linestyle = '--', linewidth = 0.5)

        print("//-------- h"+str(i)+"-histogram-saved -------------//")            
        numpy.savetxt(report_path+"h"+str(i)+"-histogram.csv", np.transpose([power_data[:,i],time_scale.T]), delimiter=",") 

        fig.savefig(report_path+"h"+str(i)+".png")
        print("Graph created and saved in directory.")
        plt.show()

    
    print("Total-Response:")
    print("Total-Power[W]:"+str(sum(h_power)))    
    print("//---------- total-histogram-saved -------------//")            
    numpy.savetxt(report_path+"total-histogram.csv", np.transpose([np.sum(power_data,axis=1),time_scale.T]), delimiter=",") 

    fig, (vax) = plt.subplots(1, 1, figsize=(12, 6))
    vax.plot(time_scale,np.sum(power_data,axis=1), 'o',markersize=2)
    vax.vlines(time_scale, [0], np.sum(power_data,axis=1),linewidth=1)

    vax.set_xlabel("time(s) \n Time resolution:"+str(tres)+"s  Bins:"+str(bins),fontsize=15)
    vax.set_ylabel('Power(W)',fontsize=15)
    vax.set_title("Total Channel Impulse Response",fontsize=20)

    vax.grid(color = 'black', linestyle = '--', linewidth = 0.5)

    fig.savefig(report_path+"total-histogram.png")
    print("Graph created and saved in directory.")
    plt.show()
    
    
    return 0

####### define input parameters for channel model ###########
#source = {tx_pos,txnormal_vector,lambert_num,power[W]}
#tx_pos: [pos_x,pos_y,pos_z]
#txnormal_vector: [pos_x,pos_y,pos_z]
s = [[1,1,2],[0,0,-1],1,1]

#receiver = {rx_pos,rxnormal_vector,area_receiver[m^2],FOV}
r = [[1,1,0],[0,0,1],1e-4,1]

#envirorment e = {reflectance,scale_factor,size_room,k_reflections}
#size_room: [x_lim,y_lim,z_lim]
e = [0.8,1/5,[2,2,2],3]    

starttime = timeit.default_timer()
#print("The start time is :",starttime)

a = tessellation(e[2][0],e[2][1],e[2][2],e[1])
b = make_parameters(a[0],e[2][0],e[2][1],e[2][2],a[1],a[2],a[3])
c = h_t(s[2],s[0],r[0],a[0][0:3,:],a[0][3,:],b,e[2][0],e[2][1],e[2][2],a[1],a[2],a[3],a[4],r[2],e[0],a[5],e[3])
create_report(c,e[3],a[6])

print("The execution time is :", timeit.default_timer() - starttime)

#print(a[6].shape)
#print(a[8])

#led_pattern(s[2])
