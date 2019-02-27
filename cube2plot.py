#! /usr/bin/env python

# collection of plotting tools

# cube2plot
# import numpy as np
# from ase.units import Bohr
# from ase.io.cube import read_cube_data
# from ase.io import write

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# from mpl_toolkits.mplot3d import Axes3D
# from skimage import measure


def cube2xyz(cube_data, cube_atoms, outfile_dat="", 
             print_range="", x_coord="", y_coord="", 
             z_coord="", mpl=True, aa=1.0, output=True, 
             p_label=" (a.u.)", show_plot=False, ax=False ):
    """
    Plots a crosssectional landscape or axial profile 
    of an ASE cube data representation 
    and the according atom position projection.
    
    Parameters
    ----------
    cube_data:
        data read from a cube file with ase.io.cube.read_cube_data(...)
    cube_atoms:
        atoms read from a cube file with ase.io.cube.read_cube_data(...)
        While .cube files conventionally use Bohr as the spatial unit,
        ASE uses Angstrom internally. Thus, atom positions and box measures
        are treated as ASE Angstrom within this data structure and converted 
        to atomic units (i.e. Bohr) for plotting in this function.
    outfile_dat:
        data of selection written to this outfile instead of stdout.
        EXPERIMENTAL
    print_range:
        'x', 'y' or 'z'. Prints all discretization points along this axis
        EXPERIMENTAL
        and exits.
    x_coord, y_coord, z_coord:
        specifies position of 2D crossection (if one of these specified)
        or 1D profile (if two of these specified)
    mpl:
        default 'True'. Plots with Matplotlib.
    aa: 
        Scaling factor. Default '1.0' results in Bohr scaling (see above)
    output:
        default 'True'. Whether or not to print some info.
    p_label:
        default ' (a.u.)'. Label for spatial axes on plot.
    show_plot:
        default 'False': Handle to newly created or modified Matplotlib
        axes object is returned, but plot not displayed.
    ax:
        default 'False'. If an axes object is given, the plot is done on these axes
        instead of creating new figure and axes.  

    Return
    ------
    axes:
        New or modified matplotlib.axes.Axes handle
    """
    
    # the adaption of this function from http://larrucea.eu/cube2xyz/ and
    # https://github.com/julenl/molecular_modeling_scripts/blob/master/cube2xyz/cube2xyz-v0.1.py
    # is just ugly. At some point, do proper implementation.
    
    import sys, os
    import numpy as np
    from ase.units import Bohr

    print("Options: \n",
        "{:>20}:{:>20}\n".format("outfile_dat",outfile_dat), 
        "{:>20}:{:>20}\n".format("print_range",print_range),
        "{:>20}:{:>20}\n".format("x_coord",x_coord),
        "{:>20}:{:>20}\n".format("y_coord",y_coord),
        "{:>20}:{:>20}\n".format("z_coord",z_coord),
        "{:>20}:{!s:>20}\n".format("mpl",mpl),
        "{:>20}:{:20.2e}\n".format("aa",aa),
        "{:>20}:{!s:>20}\n".format("output",output),
        "{:>20}:{:>20}\n".format("p_label",p_label) )
    at_coord=[]
    spacing_vec=[]
    nline=0
    values=[]
    fig=0
    # Read cube file and parse all data
    # for line in open(infile_cube,"r"):
    #    nline+=1
    #    if nline==3:
    #        try:
    #            nat=int(line.split()[0]) 
    #            origin=[float(line.split()[1]),
    ##                float(line.split()[2]),float(line.split()[3])]
    #        except:
    #            print("ERROR: non recognized cube format")
    #    elif nline >3 and nline <= 6:
    #        spacing_vec.append(line.split())
    #    elif nline > 6 and nline <= 6+nat:
    #        at_coord.append(line.split())
    #    elif nline > 5:
    #        if nline > 6+nat:
    #            for i in line.split():
    #                values.append(float(i)) 
    nat = cube_atoms.get_number_of_atoms()
    # spacing_vec = cube_atoms.cell / cube_data.shape / Bohr
    # just for now: ugly matrix of shape and spacing
    spacing_vec = np.hstack([np.atleast_2d(cube_data.shape).T,
                             cube_atoms.cell / cube_data.shape / Bohr])
    # at_coord = cube_atoms.get_positions() / Bohr
    at_coord = np.hstack( 
        [ np.atleast_2d(cube_atoms.get_atomic_numbers()).T, 
         np.zeros( ( cube_atoms.get_number_of_atoms(), 1) ),
         cube_atoms.get_positions() / Bohr ] )

    values = cube_data.flatten()
    
    # print(" ")
    def frange(x, y, jump):
        while x < y:
            yield x
            x += jump

    idx=-1

    if print_range != "":
        print(print_range+" range:")
        if print_range == "x":
            for i in range(0,int(spacing_vec[0][0])):
                print(i*float(spacing_vec[0][1])*aa)
        if print_range == "y":
            for i in range(0,int(spacing_vec[1][0])):
                print(i*float(spacing_vec[1][2])*aa)
        if print_range == "z":
            for i in range(0,int(spacing_vec[2][0])):
                print(i*float(spacing_vec[2][3])*aa)
            sys.exit()


    # The filter should work on Bohr coordinates
    filter=""  # Create a filter for the values on a segment or plane
    if x_coord != "":
        filter=filter+ " x > " \
            + str(float(x_coord)-0.9*float(spacing_vec[0][1])) + " and x < " \
            + str(float(x_coord)+0.9*float(spacing_vec[0][1])) + " "
    if y_coord != "":
        if x_coord != "": 
            filter=filter+" and "
        filter=filter+ " y > " \
            + str(float(y_coord)-0.9*float(spacing_vec[1][2]))+" and y < " \
            + str(float(y_coord)+0.9*float(spacing_vec[1][2])) + " "
    if z_coord != "":
        if len(filter)> 3:
            filter=filter+" and "
        filter=filter+ " z > " \
            + str(float(z_coord)-0.9*float(spacing_vec[2][3])) + " and z < " \
            + str(float(z_coord)+0.9*float(spacing_vec[2][3])) + " "
    if filter=="":
        filter="1==1"
        
    print("Filter: {}".format(filter))

     # Set parameter for type of plot
    plttmp=[]
    xyzs=""
    if x_coord:
        plttmp.append(x_coord)
        xyzs=xyzs+"x"
    if y_coord:
        plttmp.append(y_coord)
        xyzs=xyzs+"y"
    if z_coord:
        plttmp.append(z_coord)
        xyzs=xyzs+"z"
    plot_dim=4-len(plttmp)
    # print(" Representing "+str(plot_dim)+"D data...", xyzs) 


      #Print x,y,z,value data to stdout, file or not at all
    data=[]
    
    #print("Spacing vectors: {}".format(spacing_vec))

    if outfile_dat and output:
        #if output file name is provided, print to file instead of STD out
        tmp=sys.stdout
        sys.stdout = open(outfile_dat,'w')
    if not output:
        tmp=sys.stdout
        sys.stdout = open(os.devnull, 'w')

    print("Spacing vectors: {}".format(spacing_vec))
    for i in range(0,int(spacing_vec[0][0])):
        for j in range(0,int(spacing_vec[1][0])):
            for k in range(0,int(spacing_vec[2][0])):
                idx+=1
                x,y,z= i*float(spacing_vec[0][1]), \
                    j*float(spacing_vec[1][2]),k*float(spacing_vec[2][3])
                print("Looking at pos ({},{},{}).".format(x,y,z))
                if eval(filter):
                    print(x/aa,y/aa,z/aa, values[idx])
                    data.append([x/aa,y/aa,z/aa, values[idx]])
    
    print("Filtered {} data points.".format(len(data)))

    if outfile_dat or not output:
        sys.stdout.close()
        sys.stdout=tmp

    ylabel="Cube magnitude"
    axe_labels=["x","y","z"]
    if mpl:
        if plot_dim == 4:
            print("4D plot not implemented")
            sys.exit()

        if plot_dim == 3:
            var_axe1=list("xyz".replace(xyzs[0],""))[0]  
            var_axe2=list("xyz".replace(xyzs[0],""))[1]
            print("Plot in 3D: \n",
                "{:>20}:{:>20}\n".format("var_axe1",var_axe1), 
                "{:>20}:{:>20}\n".format("var_axe2",var_axe2),
                "{:>20}:{:>20}\n".format("axe_labels(var_axe1)",
                                         axe_labels.index(var_axe1)),
                "{:>20}:{:>20}\n".format("axe_labels(var_axe2)",
                                         axe_labels.index(var_axe2)),
                "{:>20}:{:>20}\n".format("len(data)", len(data)))
            
            from mpl_toolkits.mplot3d.axes3d import Axes3D
            from matplotlib import cm
            #X, Y =np.meshgrid(zip(*data)[axe_labels.index(var_axe1)],
            #    zip(*data)[axe_labels.index(var_axe2)])
            
            X,Y=list(zip(*data))[axe_labels.index(var_axe1)],\
                list(zip(*data))[axe_labels.index(var_axe2)]
            Z = list(zip(*data))[3]
            #print "max Z: ", max(Z)
            if not ax:
                fig = plt.figure()
                ax = fig.gca(projection='3d')
            
            # plt.hold(True)
            # surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, 
            #     cmap=cm.coolwarm, linewidth=0, antialiased=False)
            surf = ax.plot_trisurf(X, Y, Z, cmap=cm.jet, linewidth=0.2)
            fig = ax.figure
            fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
            

            ax.set_xlabel(var_axe1 + p_label)
            ax.set_ylabel(var_axe2 + p_label)
            ax.set_zlabel('Cube property')

            coord=[]
            print("at_coord: {}".format(at_coord))
            for c in at_coord: # Format atomic coordinates for plotting
                 coord.append([(float(c[2]))/aa,(float(c[3]))/aa,(float(c[4]))/aa])

            xc=list(zip(*coord))[axe_labels.index(var_axe1)]
            yc=list(zip(*coord))[axe_labels.index(var_axe2)]
            
            spanZ = max(Z) - min(Z)
            zc=[max(Z)+(spanZ*0.5)]*len(coord) #put atoms slightly higher than max val
            # projection of atoms centers above the graph
            ax.scatter(xc,yc,zc,s=150, c='b', marker='o',cmap=None, norm = None, 
                edgecolors='c', lw=3.0, alpha=1) 
           # for i in coord:

        if plot_dim == 2:  # when 2 variables are fixed a 2D plot is produced
            var_axe= "xyz".replace(xyzs[0],"").replace(xyzs[1],"")
            var_idx=axe_labels.index(var_axe)
            plot(list(zip(*data))[var_idx], list(zip(*data))[3])
            plt.xlabel(var_axe + p_label)
            plt.ylabel("Cube property")

        #plt.grid(True) 
        if show_plot:
            show()   
        
        #if fig:
        #    return fig
        if ax:
            return ax

def plotSection2D(cube_data,cube_atoms, position, axis=0, 
                  ax=False, output=False):
    """
    Plots a crosssectional landscape of an ASE cube data representation 
    and the according atom position projections.
    
    Parameters
    ----------
    cube_data:
        data read from a cube file with ase.io.cube.read_cube_data(...)
    cube_atoms:
        atoms read from a cube file with ase.io.cube.read_cube_data(...)
        While .cube files conventionally use Bohr as the spatial unit,
        ASE uses Angstrom internally. Thus, atom positions and box measures
        are treated as ASE Angstrom within this data structure and converted 
        to atomic units (i.e. Bohr) for plotting in this function.
    position:
        Position along specified axis in Bohr
        ATTENTION: Should be Bohr, but apparently still Angstrom.
        Why?
    axis:
        Axis normal to desired crossection, i.e. '0', '1' or '2', or
        equivalently 'x', 'y' or 'z'.
    ax:
        default 'False'. If an axes object is given, the plot is done on these axes
        instead of creating new figure and axes.  
        
    Return
    ------
    axes:
        New or modified matplotlib.axes.Axes handle
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from mpl_toolkits.mplot3d import Axes3D
    if not ax:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        
    if axis==0 or axis=="x":
        tmp_ax = cube2xyz(cube_data,cube_atoms, x_coord=position, 
                          mpl=True, output=output, show_plot=False, ax=ax)
    elif axis==1 or axis=="y":
        tmp_ax = cube2xyz(cube_data,cube_atoms, y_coord=position, 
                          mpl=True, output=output, show_plot=False, ax=ax)
    elif axis==2 or axis=="z":
        tmp_ax = cube2xyz(cube_data,cube_atoms, z_coord=position, 
                          mpl=True, output=output, show_plot=False, ax=ax)
        
    return tmp_ax

def sliceXYZ(cube_data,cube_atoms, sections=6, margin_relative=1e-1,cols=2,width=8,height=5):
    from ase.units import Bohr
    #from cube2xyz import cube2xyz
    #from cube2xyz import plotSection2D
    #cube_data, cube_atoms = read_cube_data(infile_cube)
    
    X = cube_atoms.cell.diagonal()/Bohr
    #X = cube_atoms.cell.diagonal() # JUST FOR NOW
    
    print("Slicing box of {} Bohr...".format(X))
    margin = margin_relative * X
    #cols = 2
    rows = round(sections/cols)
    print("Displaying {:d} sections in {:d} rows and {:d} cols.".format(sections,rows,cols))

    fig = []
    axes = []
    dim_label= ["x",",y","z"]
    for dim in range(X.shape[0]):
        #L = X.shape[3]
        x = np.linspace(margin[dim],X[dim]-margin[dim],sections)

        print("Processing {} Bohr slices in {} direction...".format(x,dim_label[dim]))

        tmp_fig, tmp_axes = plt.subplots(rows,cols,figsize=(cols*width,rows*height),subplot_kw={'projection':'3d'})
        p = 0
        for i in range(rows):
            for j in range(cols):
                if p < sections:
                    print("Calling cube2xyz.plotSection2D at slice {} = {} Bohr...".format(dim_label[dim],x[p]))

                    tmp_ax = plotSection2D(cube_data,cube_atoms, position=x[p], axis=dim, ax=tmp_axes[i,j])
                    tmp_ax.set_title(".cube slice at {} = {:.3e} Bohr".format(dim_label[dim],x[p]))
                    #cube2xyz(infile_cube, y_coord=x[p],output=False, mpl=True,ax=axes[i,j])
                    p += 1
                else:
                    break

        fig.append(tmp_fig)
        axes.append(tmp_axes)
        
    return fig

# Written by Johannes Hörmann, johannes.hoermann@imtek.uni-freiburg.de

# Adapted from Julen Larrucea's cube2xyz
# http://www.larrucea.eu 

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

def isosurfLog(data, sections=6, margin=1e-2,
               cols=2,width=8,height=5,sign=1.0,
               base=10,interval=[]):
    cube_data = data
    
    cube_data = cube_data * sign
    
    data_span = cube_data.max()-cube_data.min()
    print("Data maximum {}.".format(cube_data.max()))
    print("Data minimum {}.".format(cube_data.min()))


    print("Span of data {}.".format(data_span))

    data_offset = cube_data.min() - 1
    print("Shift data by {}.".format(data_offset))
    
    #normalized_data = (cube_data - data_offset)
    cube_log_data = np.log(cube_data - data_offset ) / np.log(base)
    data_log_max = cube_log_data.max() 
    data_log_span = cube_log_data.max() - cube_log_data.min()
    print("Span of logarithms: {}.".format(data_log_span))
    print("Magnitude of logarithm: {}.".format(data_log_max))

    normalized_log_data = cube_log_data / data_log_max
    # margin = data_log_span * margin_relative
    if not interval:
        interval = [ normalized_log_data.min(), normalized_log_data.max() ]
        
    print("Plot isosurfaces of normalized log values in interval {}".format(interval))
    normalized_log_iso_vals = np.linspace(interval[0]+margin,interval[1]-margin,sections)

    iso_vals = (base**(normalized_log_iso_vals*data_log_max) + data_offset)*sign
    #exponents = np.linspace(0+margin,np.log(data_span+offset-margin)/np.log(base),P)
    #levels = base**exponents-offset
    #if sign < 0: # descending order
    #    isopots = cube_data.max() - levels
    #else: # ascending order
    #    isopots = cube_data.min() + levels
     
    rows = round(sections/cols)

    fig, axes = plt.subplots(rows,cols,figsize=(cols*width,rows*height),subplot_kw={'projection':'3d'})
    p = 0
    for i in range(rows):
        for j in range(cols):
            if p < sections:
                #tmp_ax = fig.add_subplot(pos,projection='3d')
                #phi = cube_data.min() + round(stride*(p+0.5))
                verts, faces, normals, values = \
                    measure.marching_cubes(normalized_log_data, normalized_log_iso_vals[p])
                axes[i,j].set_title("isosurface at val = {:.3g}".format(iso_vals[p]))
                axes[i,j].plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],
                            cmap='Spectral', lw=1)
                p += 1
            else: 
                break
        
    return fig


def isosurf(data, sections=6, 
            margin=1e-2, cols=2, width=8, height=5, 
            sign=1.0, base=10, interval=[]):
    #cube_data = data
    
    data = data * sign
    
    data_span = data.max() - data.min()
    print("Data maximum {}.".format(data.max()))
    print("Data minimum {}.".format(data.min()))


    print("Span of data {}.".format(data_span))

    #data_offset = cube_data.min() - 1
    #print("Shift data by {}.".format(data_offset))
    
    #normalized_data = (cube_data - data_offset)
    #cube_log_data = np.log(cube_data - data_offset ) / np.log(base)
    #data_log_max = cube_log_data.max() 
    #data_log_span = cube_log_data.max() - cube_log_data.min()
    #print("Span of logarithms: {}.".format(data_log_span))
    #print("Magnitude of logarithm: {}.".format(data_log_max))

    #normalized_log_data = cube_log_data / data_log_max
    # margin = data_log_span * margin_relative
    if not interval:
        interval = [ data.min(), data.max() ]
        
    #print("Plot isosurfaces of normalized log values in interval {}".format(interval))
    iso_vals = np.linspace(interval[0]+margin,interval[1]-margin,sections)

    #iso_vals = (base**(normalized_log_iso_vals*data_log_max) + data_offset)*sign
    #exponents = np.linspace(0+margin,np.log(data_span+offset-margin)/np.log(base),P)
    #levels = base**exponents-offset
    #if sign < 0: # descending order
    #    isopots = cube_data.max() - levels
    #else: # ascending order
    #    isopots = cube_data.min() + levels
     
    rows = round(sections/cols)

    fig, axes = plt.subplots(rows,cols,
                             figsize=(cols*width,rows*height),
                             subplot_kw={'projection':'3d'})
    p = 0
    for i in range(rows):
        for j in range(cols):
            if p < sections:
                #tmp_ax = fig.add_subplot(pos,projection='3d')
                #phi = cube_data.min() + round(stride*(p+0.5))
                verts, faces, normals, values = \
                    measure.marching_cubes(data, iso_vals[p])
                axes[i,j].set_title("isosurface at val = {:.3g}".format(iso_vals[p]))
                axes[i,j].plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],
                            cmap='Spectral', lw=1)
                p += 1
            else: 
                break
        
    return fig