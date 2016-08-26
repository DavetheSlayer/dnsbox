import h5py
import numpy
from mayavi import mlab


f = h5py.File('vorticity.h5', 'r')
Nx, Ny, Nz = f.get('omegax').shape[:]

Lx = 2.0 * np.pi
Ly = 2.0 * np.pi
Lz = 2.0 * np.pi

x = [i * (Lx/Nx) for i in xrange(- Nx / 2, 1 + Nx / 2)]
y = [i * (Ly/Ny) for i in xrange(- Ny / 2, 1 + Ny / 2)] 
z = [i * (Lz/Nz) for i in xrange(- Nz / 2, 1 + Nz / 2)] 

omegatot = np.sqrt(f.get('omegax')[:, :, :]**2 \
                 + f.get('omegay')[:, :, :]**2 \
                 + f.get('omegaz')[:, :, :]**2)

f.close()

xx = np.zeros((Nx, Ny, Nz), dtype = 'float')
yy = np.zeros((Nx, Ny, Nz), dtype = 'float')
zz = np.zeros((Nx, Ny, Nz), dtype = 'float')

for i in xrange(Nx):
    for j in xrange(Ny):
        for k in xrange(Nz):
            xx[i, j, k] = x[i]
            yy[i, j, k] = y[j]
            zz[i, j, k] = z[k]


src = mlab.pipeline.scalar_field(xx, yy, zz, omegatot, colormap='YlGnBu')

mlab.pipeline.iso_surface(src, contours=[omegatot.min() 
                                         + 0.1 * omegatot.ptp(), ],
                          colormap='YlGnBu', opacity=0.2)

#mlab.pipeline.iso_surface(src, contours=[omegatot.max() 
                                         #- 0.5 * omegatot.ptp(), ],
                          #colormap='YlGnBu', opacity=1.0)

mlab.pipeline.iso_surface(src, contours=[omegatot.max() 
                                         - 0.5 * omegatot.ptp(), ],
                          colormap='YlGnBu', opacity=0.5)

#mlab.pipeline.image_plane_widget(src, plane_orientation='y_axes', 
                                 #slice_index = Ny/2, colormap='YlGnBu', 
                                 #opacity=0.01)

#mlab.pipeline.image_plane_widget(src, plane_orientation='x_axes', 
                                 #slice_index = Nx/2, colormap='YlGnBu', 
                                 #opacity=0.01)

mlab.scalarbar()
mlab.xlabel('x', object=src)
mlab.ylabel('y', object=src)
mlab.zlabel('z', object=src)
