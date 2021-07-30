# RecursiveModel-VLC

The RecursiveModel-VLC is an script to compute the channel impulse response of optical link in the visible range using a recursive model presented in [1]. 

## Installation

Make sure have installed all libraries used in the python script. Once the script execution has finished, the power and delay for every ray  are in 'cir' folder as .csv files. In the 'report' folder are the graph of the channel impulse response for every reflection. 

## Usage

To run this recursive model you only need run the main_model.py.  

This model uses a set of parameters to define features of source, receiver and environment, as follows:

```python

# list to define the source features
s = {tx_pos,txnormal_vector,lambert_num,power[W]}
# example 
s = [[1,1,2],[0,0,-1],1,1]

# list to define the receiver features
r = {rx_pos,rxnormal_vector,area_receiver[m^2],FOV}
# example
r = [[1,1,0],[0,0,1],1e-4,1]

# list to define the environment features
e = {reflectance,scale_factor,size_room,k_reflections}
# example
e = [0.8,1/71,[2,2,2],1]
```

The previous values for source, receiver and environment are the default values. Using this set of parameters, the script compute the channel impulse response based on the follow functions:

```python
# The tessellation function calculates the coordinates of every points in the walls discretization, and returns the array_points.
tessellation(x_lim,y_lim,z_lim,scale_factor):

#This function creates an cross-parmeteres array between points, cosine of output angle and euclidean distance. Returns the cross-parameters array.
make_parameters(array_points,x_lim,y_lim,z_lim,no_xtick,no_ytick,no_ztick):

#This funciton computes the channel impulse response from cross-parameters array, based on number of reflection. Returns a list with the different order response, h0,h1,h2...hk. 
h_t(m,tx_pos,rx_pos,points,parameters,x_lim,y_lim,z_lim,no_xtick,no_ytick,no_ztick,init_index,a_r,rho,delta_A,k_reflec):

#This function creates an analysis of the simulation, the total power for each reflection, total power in the receiver and plots. 
create_report(h_k,k_reflec):
```


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Referencies
[1] Barry, J. R., Kahn, J. M., Krause, W. J., Lee, E. A., & Messerschmitt, D. G. (1993). Simulation of multipath impulse response for indoor wireless optical channels. IEEE journal on selected areas in communications, 11(3), 367-379.
