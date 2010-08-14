"""
A module to aid in volumetric reconstruction from light fields
"""

import sys
import os
#sys.path.append("/Users/todicus/Code/Python/MLSL")

if os.uname()[0]=='Darwin': # osx paths
    import mlsl
    import mlsl2
else:	# this might be seriously f'd code, but it seems to work...
    import linux.mlsl
    import linux.mlsl2
    mlsl = linux.mlsl
    mlsl2 = linux.mlsl2


from poly import Polynomial as P

import array
import struct
import math
import time
import thread
import threading
import collections

class FeasibilityError(Exception):
    "Notify barrier method that the current t causes feasibility errors in a Newton step"
    def __init__(self, backtrack):
        Exception.__init__(self, 'feasibility backtrack error: %g' % backtrack)
        self.backtrack = backtrack

class Error(Exception):
    pass

class LFReconGlobals(object):
    """
    A singleton class that handles settings
    """
    def __new__(cls, *args, **kwds):
        instance = cls.__dict__.get('__instance__')
        if instance is not None:
            return instance
        instance = object.__new__(cls)
        cls.__instance__ = instance
        instance.init(*args, **kwds)
        return instance
    
    def init(self, *args, **kwds):
        "Initialize the very first time"
        self.threads = 2
        self.verbose = {}
        self.logfile = sys.stderr
        self.focal_alg = ''

class LFRecon:
    def __init__(self,
                 xlen, ylen, zlen, ulen, vlen,
                 spec_mag, spec_na, spec_n,
                 ulens_pitch, ulens_flen, ulens_n,
                 z_center, z_spacing, radiometry=None):
        """
        Initialise reconstruction operations

        (xlen,ylen,ulen,vlen) are the dimensions (in pixels)
        of the light field
        (xlen,ylen,zlen) are the dimensions in voxels of the
        volume reconstruction
        spec_mag is the magnification of the objective
        spec_na is the NA of the objective
        spec_n is the index of refraction in the specimen
        ulens_pitch is the pitch of the microlens array
                    (each microlens is ulens_pitch x ulens_pitch
                     micron in size)
        ulens_flen is the focal length of the microlens array
                   (specified in micron)
        ulens_n is the index of refraction directly outside the
                   microlens
        z_center is the z offset in micron from the focal plane
                 of the center slice in the volume
                 positive z is "up" from the focal plane towards
                 the objective.  If there are an even number of
                 slices in the volume, then this is the average
                 of the z offset of the two middle slices
        z_spacing is the number of micron between slices in the
                  volume
        radiometry is either None or a mapping (u,v)->float 
                   that specifies the intensity of a uniform
                   fluorescent volume
                   along each of the ray angles.  All the (u,v)
                   ray angles that will be used in the reconstruction
                   must be specified in this mapping
        threads is the maximum number of threads to use when performing
                calculations
        """
        # first, store the parameters
        self.xlen = xlen
        self.ylen = ylen
        self.zlen = zlen
        self.ulen = ulen
        self.vlen = vlen
        self.spec_mag = spec_mag
        self.spec_na = spec_na
        self.spec_n = spec_n
        self.ulens_pitch = ulens_pitch
        self.ulens_flen = ulens_flen
        self.ulens_n = ulens_n
        self.z_center = z_center
        self.z_spacing = z_spacing

        _dummy = array.array('f')
        self._strides = _dummy.itemsize, _dummy.itemsize, _dummy.itemsize*self.xlen

        # calculate z slice values
        self.compute_z_coords()

        # set up radiometry
        self.compute_radiometry(radiometry)

        # set up ray spread functions for each slice
        self.compute_rayspreads()

        if self.focal_alg() == 'psf':
            # set up the psf
            self.compute_psf()
        elif self.focal_alg() == 'sipsf':
            # set up shift-invariant psf
            self.compute_sipsf()

        self.compute_focal_stack = self.compute_focal_stack_fast
        self.compute_lf_projection = self.compute_lf_projection_fast

    def threads(self):
        return LFReconGlobals().threads

    def verbose(self):
        return LFReconGlobals().verbose

    def logfile(self):
        return LFReconGlobals().logfile

    def focal_alg(self):
        return LFReconGlobals().focal_alg

    def log(self, flag, s, encoding='utf8', newline=True):
        if flag in self.verbose():
            t = time.time()
            time_str = time.strftime('%Y-%m-%d %H:%M:%S.', time.localtime(t)) + ('%f' % (t-math.floor(t)))[2:]
            output_str = self.verbose()[flag] + ':' + time_str + ': ' + s
            if newline:
                output_str += '\n'
            output_str = output_str.encode(encoding)
            self.logfile().write(output_str)
            self.logfile().flush()

    # Can alter this method to generate a non-linear series of z-coordinates
    def compute_z_coords(self):
        "Compute the z values for each slice in the volume"
        self.log('i', 'Computing z values for each slice in the volume')
        z_start = self.z_center-0.5*(self.zlen - 1)*self.z_spacing
        self.z_coords = [x*self.z_spacing+z_start for x in range(self.zlen)]

    def compute_radiometry(self, desired_radiometry=None):
        """
        Compute the scaling needed in the projection->light field stage for
        each ray angle
        """
        # the natural radiometry of the rayspread at z=0
        self.radiometry_base = {}
        # what we want the radiometry to be for each ray angle
        self.radiometry_desired = {}
        # the total scale to be applied to every ray angle to
        # make the focal stack operation be relatively energy preserving
        self.radiometry_scale = 0

        self.log('i','Computing radiometry for each ray angle')

        # compute "natural" radiometry
        for v in range(self.vlen):
            for u in range(self.ulen):
                spread = mlsl.rayspread(0, u, v, self.ulen, self.vlen,
                                        self.ulens_pitch, self.ulens_flen, self.ulens_n,
                                        self.spec_mag, self.spec_na, self.spec_n)
                psf = array.array('f',spread[0])
                # compute base radiometry
                base_radiometry = sum(psf)
                if base_radiometry:
                    self.radiometry_base[(u,v)] = base_radiometry
       
        if desired_radiometry is not None:
            # copy over desired radiometry if present
            for (u,v) in self.radiometry_base:
                if (u,v) in desired_radiometry:
                    self.radiometry_desired[(u,v)] = desired_radiometry[(u,v)]
        else:
            self.radiometry_desired.update(self.radiometry_base)
        
        # rescale
        self.radiometry_scale = 1.0 / math.sqrt(sum(self.radiometry_desired.values()) * self.zlen)
        
    def compute_rayspreads(self):
        """
        Compute rayspread functions for each ray angle and z slice
        """
        self.log('i','Computing rayspread functions')
        self.rayspreads = {}
        self.flipped_rayspreads = {}
        for z_index in range(self.zlen):
            z = self.z_coords[z_index]
            self.log('i', 'z = %.3g um' % z)
            for (u,v) in self.radiometry_desired:
                spread = mlsl.rayspread(z, u, v, self.ulen, self.vlen,
                                        self.ulens_pitch, self.ulens_flen, self.ulens_n,
                                        self.spec_mag, self.spec_na, self.spec_n)
                psf = array.array('f',spread[0])
                # rescale
                total = sum(psf)
                if total:
                    scale = self.radiometry_scale * self.radiometry_desired[(u,v)] / total
                    for i in range(len(psf)):
                        psf[i] = psf[i] * scale
                # store
                self.rayspreads[(z_index,u,v)] = (psf, spread[1], spread[2], spread[3])
                # flip
                flipped = mlsl.flip_kernel((psf, spread[1], spread[2], spread[3]),True,True)
                self.flipped_rayspreads[(z_index,u,v)] = flipped

    def compute_psf(self):
        """
        Compute the z-varying PSF
        """
        self.log('i','Computing z-varying PSF')

        psf = {}

        psf_lock = threading.RLock()
        psf_done = threading.Semaphore(0)

        def thread_func(zs):
            try:
                local_psf = {}
                for source_z_index in zs:
                    source_z = self.z_coords[source_z_index]
                    self.log('i', '%.3g um -> * transfer' % source_z)
                    for dest_z_index in range(self.zlen):
                        spread = ('\x00\x00\x00\x00', 0, 1, 0) # start with a zero kernel
                        for (u,v) in self.radiometry_desired:
                            source_spread = self.rayspreads[(source_z_index,u,v)]
                            dest_spread = self.rayspreads[(dest_z_index,u,v)]
                            corr_spread = mlsl.correlate_kernel(source_spread,dest_spread)
                            spread = mlsl.mix_kernel(1.0, spread, 1.0, corr_spread)
                        local_psf[(source_z_index,dest_z_index)] = spread

                psf_lock.acquire()
                try:
                    psf.update(local_psf)
                finally:
                    psf_lock.release()
            finally:
                psf_done.release()

        # run the threads
        threads = min(self.threads(), self.zlen)
        for i in range(threads):
            start = int(round(1.0 * self.zlen * i / threads))
            end = int(round(1.0 * self.zlen * (i+1) / threads))
            thread.start_new_thread(thread_func, (range(start,end),))
        for i in range(threads):
            psf_done.acquire()

        self.psf = psf

    def compute_sipsf(self):
        """
        Compute the spatially invariant PSF
        """
        self.log('i','Computing the spatially invariant PSF')

        sipsf = {}

        z_offsets = range(-self.zlen+1,self.zlen)

        for offset in z_offsets:
            z = offset * self.z_spacing
            self.log('i', '%.3g um -> 0 um transfer' % z)
            accum = ('\x00\x00\x00\x00', 0, 1, 0) # start with a zero kernel
            for (u,v) in self.radiometry_desired:
                spread = mlsl.rayspread(z, u, v, self.ulen, self.vlen,
                                        self.ulens_pitch, self.ulens_flen, self.ulens_n,
                                        self.spec_mag, self.spec_na, self.spec_n)
                psf = array.array('f',spread[0])
                # rescale
                total = sum(psf)
                if total:
                    scale = self.radiometry_scale * self.radiometry_scale * self.radiometry_desired[(u,v)] / total
                    for i in range(len(psf)):
                        psf[i] = psf[i] * scale 
                # accumulate
                accum = mlsl.mix_kernel(1.0, accum, 1.0, (psf, spread[1], spread[2], spread[3]))
            sipsf[offset] = accum

        self.sipsf = sipsf

    def blank_volume(self):
        """
        Return a writeable array that can hold a blank volume, filled with zeros
        """
        return array.array('f', struct.pack('f',0.0)*self.xlen*self.ylen*self.zlen)

    def const_volume(self, value):
        """
        Return a writeable array that can hold a volume, filled with a constant value
        """
        return array.array('f', struct.pack('f',value)*self.xlen*self.ylen*self.zlen)

    def blank_lf(self):
        """
        Return a writeable array that can hold a blank light field, filled with zeros
        """
        return array.array('f', struct.pack('f',0.0)*self.xlen*self.ylen*self.ulen*self.vlen)

    def blank_image(self):
        """
        Return a writeable array that can hold an image, filled with zeros
        """
        return array.array('f', struct.pack('f',0.0)*self.xlen*self.ylen)

    def copy(self, original, dest=None):
        return array.array('f', original.tostring())

    def volume_byte_offset(self, z_index):
        """
        Return the byte offset of a specific slice in the volume
        """
        return z_index * self.xlen * self.ylen * self._strides[1]

    def lf_byte_offset(self, u_index, v_index):
        return ( u_index + self.ulen * v_index ) * self.xlen * self.ylen * self._strides[1]

    def set_voxel(self, vol, x, y, z, value):
        """
        Set a single voxel value
        """
        vol[x+self.xlen*(y+self.ylen*z)] = value

    def get_voxel(self, vol, x, y, z):
        """
        Get the value of a single voxel
        """
        return vol[x+self.xlen*(y+self.ylen*z)] 

    def image_dot(self, A, A_offset, B, B_offset):
        "Takes the element wise dot product of two images"
        return mlsl.image_dot('f', A, A_offset,
                              self._strides[0], self._strides[1], self._strides[2],
                              B, B_offset,
                              self._strides[0], self._strides[1], self._strides[2],
                              1, self.xlen, self.ylen)

    def image_log(self, A, A_offset, B_offset, B=None):
        "1./A -> B"
        return mlsl.image_log('f', A, A_offset,
                                     self._strides[0], self._strides[1], self._strides[2],
                                     B_offset,
                                     self._strides[0], self._strides[1], self._strides[2],
                                     1, self.xlen, self.ylen, B)

    def image_mix(self, alpha, A, A_offset, beta, B, B_offset, C_offset, C=None):
        "alpha*A + beta*B -> C"
        return mlsl.image_mix('f', alpha, A, A_offset,
                              self._strides[0], self._strides[1], self._strides[2],
                              beta, B, B_offset,
                              self._strides[0], self._strides[1], self._strides[2],
                              C_offset,
                              self._strides[0], self._strides[1], self._strides[2],
                              1, self.xlen, self.ylen, C)

    def image_mult(self, A, A_offset, B, B_offset, C_offset, C=None):
        "A.*B -> C"
        return mlsl.image_mult('f', A, A_offset,
                               self._strides[0], self._strides[1], self._strides[2],
                               B, B_offset,
                               self._strides[0], self._strides[1], self._strides[2],
                               C_offset,
                               self._strides[0], self._strides[1], self._strides[2],
                               1, self.xlen, self.ylen, C)

    def image_reciprocal(self, A, A_offset, B_offset, B=None):
        "1./A -> B"
        return mlsl.image_reciprocal('f', A, A_offset,
                                     self._strides[0], self._strides[1], self._strides[2],
                                     B_offset,
                                     self._strides[0], self._strides[1], self._strides[2],
                                     1, self.xlen, self.ylen, B)

    def image_scale(self, alpha, beta, A, A_offset, B_offset, B=None):
        "alpha + beta*A -> B"
        return mlsl.image_scale(alpha, beta, A, 'f', A_offset,
                                self._strides[0], self._strides[1], self._strides[2],
                                'f', B_offset,
                                self._strides[0], self._strides[1], self._strides[2],
                                1, self.xlen, self.ylen, B)

    def volume_dot(self, A, B):
        "Takes the element wise dot product of two volumes"
        return mlsl.image_dot('f', A, 0,
                              self._strides[0], self._strides[1], self._strides[2],
                              B, 0,
                              self._strides[0], self._strides[1], self._strides[2],
                              1, self.xlen, self.ylen*self.zlen)

    def volume_log(self, A, B=None):
        "1./A -> B"
        return mlsl.image_log('f', A, 0,
                              self._strides[0], self._strides[1], self._strides[2],
                              0,
                              self._strides[0], self._strides[1], self._strides[2],
                              1, self.xlen, self.ylen*self.zlen, B)
        
    def volume_mix(self, alpha, A, beta, B, C=None):
        "alpha*A + beta*B -> C"
        return mlsl.image_mix('f', alpha, A, 0,
                              self._strides[0], self._strides[1], self._strides[2],
                              beta, B, 0,
                              self._strides[0], self._strides[1], self._strides[2],
                              0,
                              self._strides[0], self._strides[1], self._strides[2],
                              1, self.xlen, self.ylen*self.zlen, C)

    def volume_mult(self, A, B, C=None):
        "A.*B -> C"
        return mlsl.image_mult('f', A, 0,
                               self._strides[0], self._strides[1], self._strides[2],
                               B, 0,
                               self._strides[0], self._strides[1], self._strides[2],
                               0,
                               self._strides[0], self._strides[1], self._strides[2],
                               1, self.xlen, self.ylen*self.zlen, C)

    def volume_reciprocal(self, A, B=None):
        "1./A -> B"
        return mlsl.image_reciprocal('f', A, 0,
                                     self._strides[0], self._strides[1], self._strides[2],
                                     0,
                                     self._strides[0], self._strides[1], self._strides[2],
                                     1, self.xlen, self.ylen*self.zlen, B)

    def volume_scale(self, alpha, beta, A, B=None):
        "alpha + beta*A -> B"
        return mlsl.image_scale(alpha, beta, A, 'f', 0,
                                self._strides[0], self._strides[1], self._strides[2],
                                'f', 0,
                                self._strides[0], self._strides[1], self._strides[2],
                                1, self.xlen, self.ylen*self.zlen, B)

    def lf_dot(self, A, B):
        "Takes the element wise dot product of two volumes"
        return mlsl.image_dot('f', A, 0,
                              self._strides[0], self._strides[1], self._strides[2],
                              B, 0,
                              self._strides[0], self._strides[1], self._strides[2],
                              1, self.xlen, self.ylen*self.ulen*self.vlen)

    def lf_log(self, A, B=None):
        "1./A -> B"
        return mlsl.image_log('f', A, 0,
                              self._strides[0], self._strides[1], self._strides[2],
                              0,
                              self._strides[0], self._strides[1], self._strides[2],
                              1, self.xlen, self.ylen*self.ulen*self.vlen, B)
        
    def lf_mix(self, alpha, A, beta, B, C=None):
        "alpha*A + beta*B -> C"
        return mlsl.image_mix('f', alpha, A, 0,
                              self._strides[0], self._strides[1], self._strides[2],
                              beta, B, 0,
                              self._strides[0], self._strides[1], self._strides[2],
                              0,
                              self._strides[0], self._strides[1], self._strides[2],
                              1, self.xlen, self.ylen*self.ulen*self.vlen, C)

    def lf_mult(self, A, B, C=None):
        "A.*B -> C"
        return mlsl.image_mult('f', A, 0,
                               self._strides[0], self._strides[1], self._strides[2],
                               B, 0,
                               self._strides[0], self._strides[1], self._strides[2],
                               0,
                               self._strides[0], self._strides[1], self._strides[2],
                               1, self.xlen, self.ylen*self.ulen*self.vlen, C)

    def lf_reciprocal(self, A, B=None):
        "1./A -> B"
        return mlsl.image_reciprocal('f', A, 0,
                                     self._strides[0], self._strides[1], self._strides[2],
                                     0,
                                     self._strides[0], self._strides[1], self._strides[2],
                                     1, self.xlen, self.ylen*self.ulen*self.vlen, B)

    def lf_scale(self, alpha, beta, A, B=None):
        "alpha + beta*A -> B"
        return mlsl.image_scale(alpha, beta, A, 'f', 0,
                                self._strides[0], self._strides[1], self._strides[2],
                                'f', 0,
                                self._strides[0], self._strides[1], self._strides[2],
                                1, self.xlen, self.ylen*self.ulen*self.vlen, B)

    
    def compute_focal_slice(self, lf, volume, z_index):
        """
        Compute a focal slice using the input lf and write it to volume as well as returning it
        The z value is given by self.z_coords[z_index]
        """
        flipped_rayspreads = [(u,v,self.flipped_rayspreads[(z,u,v)]) for (z,u,v) in self.flipped_rayspreads.keys() if z == z_index]
        # create a blank image to integrate into
        accum = self.blank_image()
        temp = self.blank_image()
        for (u,v,flipped_rayspread) in flipped_rayspreads:
            #rayspread = mlsl.flip_kernel(rayspread,True,True)
            temp = mlsl2.lf_blur_subaperture(lf, 1, 'f', 'xyuv',
                                             1, self.xlen, self.ylen, self.ulen, self.vlen, u, v,
                                             flipped_rayspread[0], flipped_rayspread[1], flipped_rayspread[2], flipped_rayspread[3],
                                             'f', 0, struct.calcsize('f'), struct.calcsize('f')*self.xlen, temp, threads=1)
            accum = self.image_mix(1.0, accum, 0,
                                   1.0, temp, 0,
                                   0, accum)
        # now copy the accumulated image into the volume
        self.image_mix(0.0, accum, 0,
                       1.0, accum, 0,
                       self.volume_byte_offset(z_index), volume)
        return accum

    def compute_focal_slices(self, lf, volume, z_indices):
        """
        Compute several focal slices indicated by z_indices using the input lf and write it to volume
        """
        mlsl.compute_focal_stack(lf, volume, 0, 0, self.xlen, self.ylen, self.zlen, self.ulen, self.vlen,
                                 self._strides[1], self._strides[2], self._strides[2]*self.ylen, self._strides[2]*self.ylen*self.ulen,
                                 self._strides[1], self._strides[2], self._strides[2]*self.ylen,
                                 1.0, self.flipped_rayspreads, z_indices)

    def compute_lf_subaperture(self, volume, lf, u, v):
        """
        Compute a subaperture image of a light field from a volume
        """
        rayspreads = [(_z,self.rayspreads[(_z,_u,_v)]) for (_z,_u,_v) in self.rayspreads.keys() if (_u,_v) == (u,v)]
        accum = self.blank_image()
        temp = self.blank_image()
        for (z_index, rayspread) in rayspreads:
            mlsl.correlate2(volume, 'f', self.volume_byte_offset(z_index),
                           0, 0,
                           struct.calcsize('f'), struct.calcsize('f')*self.xlen,
                           self.xlen, self.ylen,
                           'f', 0,
                           0, 0,
                           struct.calcsize('f'), struct.calcsize('f')*self.xlen,
                           self.xlen, self.ylen,
                           self.xlen, self.ylen,
                           rayspread[0], rayspread[1], rayspread[2], rayspread[3],
                           output=temp)
            accum = self.image_mix(1.0, accum, 0,
                                   1.0, temp, 0,
                                   0, accum)
        # now copy the accumulated image into the light field
        self.image_mix(0.0, accum, 0,
                       1.0, accum, 0,
                       self.lf_byte_offset(u,v), lf)
        return accum

    def compute_lf_subapertures(self, volume, lf, uvs):
        """
        Compute several focal slices indicated by z_indices using the input lf and write it to volume
        """
        mlsl.compute_lf_projection(volume, lf, 0, 0, self.xlen, self.ylen, self.zlen, self.ulen, self.vlen,
                                   self._strides[1], self._strides[2], self._strides[2]*self.ylen,
                                   self._strides[1], self._strides[2], self._strides[2]*self.ylen, self._strides[2]*self.ylen*self.ulen,
                                   1.0, self.rayspreads, uvs)
        
    def compute_focal_stack_slow(self, lf, volume):
        self.log('f','Computing focal stack...')

        done = threading.Semaphore(0)

        def thread_func(slice_indices):
            try:
                for z_index in slice_indices:
                    self.compute_focal_slice(lf, volume, z_index)
                    self.log('f', 'Slice #%d (z=%g) done' % (z_index, self.z_coords[z_index]))
            finally:
                done.release()

        threads = min(self.threads(), self.zlen)
        for i in range(threads):
            start = int(round(1.0 * self.zlen * i / threads))
            end = int(round(1.0 * self.zlen * (i+1) / threads))
            thread.start_new_thread(thread_func, (range(start,end),))
        for i in range(threads):
            done.acquire()

    def compute_focal_stack_fast(self, lf, volume):
        self.log('f','Computing focal stack...')

        done = threading.Semaphore(0)

        def thread_func(slice_indices):
            try:
                self.compute_focal_slices(lf, volume, slice_indices)
            finally:
                done.release()

        threads = min(self.threads(), self.zlen)
        for i in range(threads):
            thread.start_new_thread(thread_func, (range(i, self.zlen, threads),))
        for i in range(threads):
            done.acquire()

    def compute_lf_projection_slow(self, volume, lf):
        self.log('p', 'Computing projected light field...')

        done = threading.Semaphore(0)
    
        def thread_func(uvs):
            try:
                for (u,v) in uvs:
                    self.compute_lf_subaperture(volume, lf, u, v)
                    self.log('p','Subaperture (%d,%d) done' % (u, v))
            finally:
                done.release()
                
        full_uvs = sorted(self.radiometry_desired.keys())
        len_full_uvs = len(full_uvs)
        threads = min(self.threads(), len_full_uvs)
        for i in range(threads):
            start = int(round(1.0 * len_full_uvs * i / threads))
            end = int(round(1.0 * len_full_uvs * (i+1) / threads))
            thread.start_new_thread(thread_func, (full_uvs[start:end],))
        for i in range(threads):
            done.acquire()

    def compute_lf_projection_fast(self, volume, lf):
        self.log('p', 'Computing projected light field...')

        done = threading.Semaphore(0)
    
        def thread_func(uvs):
            try:
                self.compute_lf_subapertures(volume, lf, uvs)
            finally:
                done.release()
                
        full_uvs = sorted(self.radiometry_desired.keys())
        len_full_uvs = len(full_uvs)
        threads = min(self.threads(), len_full_uvs)
        for i in range(threads):
            thread.start_new_thread(thread_func, ([full_uvs[x] for x in range(i, len_full_uvs, threads)],))
        for i in range(threads):
            done.acquire()


    def apply_focal_operator_lf(self, input_vol, output_vol, lf=None):
        "Create a focal stack from a fluorescent density volume, using light fields"
        self.log('F', 'Applying focal operator using LF...')

        if lf is None:
            lf = self.blank_lf()

        self.compute_lf_projection(input_vol, lf)
        self.compute_focal_stack(lf, output_vol)

    def apply_focal_operator_psf(self, input_vol, output_vol):
        "Create a focal stack from a fluorescent density volume, using computed psfs"
        self.log('F', 'Applying focal operator using PSF...')

        # set up threads
        done = threading.Semaphore(0)

        def thread_func(slice_indices):
            try:
                # create output image for transfer
                temp = self.blank_image()
                for z_index in slice_indices:
                    # create accumulator
                    accum = self.blank_image()
                    for src_index in range(self.zlen):
                        # get transfer from source to dest
                        transfer = self.psf[(src_index, z_index)]
                        # compute transfer
                        mlsl.correlate2(input_vol, 'f', self.volume_byte_offset(src_index),
                                        0, 0,
                                        struct.calcsize('f'), struct.calcsize('f')*self.xlen,
                                        self.xlen, self.ylen,
                                        'f', 0,
                                        0, 0,
                                        struct.calcsize('f'), struct.calcsize('f')*self.xlen,
                                        self.xlen, self.ylen,
                                        self.xlen, self.ylen,
                                        transfer[0], transfer[1], transfer[2], transfer[3],
                                        output=temp)
                        # accumulate
                        self.image_mix(1.0, accum, 0,
                                       1.0, temp, 0,
                                       0, accum)
                    # copy into volume
                    self.image_mix(0.0, accum, 0,
                                   1.0, accum, 0,
                                   self.volume_byte_offset(z_index), output_vol)
            finally:
                done.release()

        # run the threads
        threads = min(self.threads(), self.zlen)
        for i in range(threads):
            start = int(round(1.0 * self.zlen * i / threads))
            end = int(round(1.0 * self.zlen * (i+1) / threads))
            thread.start_new_thread(thread_func, (range(start,end),))
        for i in range(threads):
            done.acquire()

    def apply_focal_operator_sipsf(self, input_vol, output_vol):
        "Create a focal stack from a fluorescent density volume, using shift invariant psfs"
        self.log('F', 'Applying focal operator using SIPSF...')

        # set up threads
        done = threading.Semaphore(0)

        def thread_func(slice_indices):
            try:
                # create output image for transfer
                temp = self.blank_image()
                for z_index in slice_indices:
                    # create accumulator
                    accum = self.blank_image()
                    for src_index in range(self.zlen):
                        # get transfer from source to dest
                        transfer = self.sipsf[src_index - z_index]
                        # compute transfer
                        mlsl.correlate2(input_vol, 'f', self.volume_byte_offset(src_index),
                                        0, 0,
                                        struct.calcsize('f'), struct.calcsize('f')*self.xlen,
                                        self.xlen, self.ylen,
                                        'f', 0,
                                        0, 0,
                                        struct.calcsize('f'), struct.calcsize('f')*self.xlen,
                                        self.xlen, self.ylen,
                                        self.xlen, self.ylen,
                                        transfer[0], transfer[1], transfer[2], transfer[3],
                                        output=temp)
                        # accumulate
                        self.image_mix(1.0, accum, 0,
                                       1.0, temp, 0,
                                       0, accum)
                    # copy into volume
                    self.image_mix(0.0, accum, 0,
                                   1.0, accum, 0,
                                   self.volume_byte_offset(z_index), output_vol)
            finally:
                done.release()

        # run the threads
        threads = min(self.threads(), self.zlen)
        for i in range(threads):
            start = int(1.0 * self.zlen * i / threads)
            end = int(1.0 * self.zlen * (i+1) / threads)
            thread.start_new_thread(thread_func, (range(start,end),))
        for i in range(threads):
            done.acquire()

    def dump_volume(self, volume):
        "Return an ImageStack TMP representation of a volume"
        return struct.pack('IIII',self.zlen,self.xlen,self.ylen,1) + volume.tostring()

    def dump_lf(self, lf):
        "Return an ImageStack TMP representation of a light field"
        return struct.pack('IIII',self.ulen*self.vlen,self.xlen,self.ylen,1) + lf.tostring()

    def load_volume(self, s):
        "Load an ImageStack TMP representation of a volume"
        zlen,xlen,ylen,chans = struct.unpack('IIII', s[0:16])
        if (zlen,xlen,ylen,chans) != (self.zlen,self.xlen,self.ylen,1):
            raise Error('Not a compatible volume: file was %dx%dx%dx%d, required %dx%dx%dx%d' % (zlen,xlen,ylen,chans,self.zlen,self.xlen,self.ylen,1))
        volume = array.array('f', s[16:])
        return volume

    def load_lf(self, s):
        "Load an ImageStack TMP representation of a light field"
        uvlen,xlen,ylen,chans = struct.unpack('IIII', s[0:16])
        if (uvlen,xlen,ylen,chans) != (self.ulen*self.vlen,self.xlen,self.ylen,1):
            raise Error('Not a compatible light field: file was %dx%dx%dx%d, required %dx%dx%dx%d' % (uvlen,xlen,ylen,chans,self.ulen*self.vlen,self.xlen,self.ylen,1))
        lf = array.array('f', s[16:])
        return lf

    def barrier_recon_1(self, lf, ts, newton_alpha, newton_beta, newton_max, cg_max, newton_e, newton_e2, cg_e, x = None, fs = None):
        """
        Run the barrier method on an input light field to solve for the volume it represents

        lf is the input light field
        ts is a finite sequence of the values of t for barrier method
        newton_alpha is the backtracking line search alpha
        newton_beta is the backtracking line search beta
        newton_max is the maximum number of Newton steps
        cg_max is the maximum number of CG steps
        newton_e is the epsilon for the Newton steps
        cg_e is the epsilon for the CG steps
        x is an initial (positive) guess and will be modified
        fs is the focal stack computed from lf and will be computed if not provided
        """
        self.log('b', 'ts = %s' % repr(ts))
        self.log('b', 'newton_alpha = %g, newton_beta = %g, newton_max = %g, cg_max = %g, newton_e = %g, cg_e = %g' % (newton_alpha, newton_beta, newton_max, cg_max, newton_e, cg_e))
        if x is None:
            const_level = 1.0 / math.sqrt(2*ts[0])
            self.log('b', 'Initializing input to barrier method to be at level %g' % const_level)
            x = self.const_volume(const_level)
        if fs is None:
            fs = self.blank_volume()
            self.compute_focal_stack(lf, fs)
        for t in ts:
            self.log('b', 'Barrier method t = %g, running Newton' % t)
            self.newton_recon(lf, fs, x, t, newton_alpha, newton_beta, newton_max, cg_max, newton_e, newton_e2, cg_e, throw_feas = False)
        self.log('b', 'Barrier method finished')
        return x

    def barrier_recon_2(self, lf, ts, newton_alpha, newton_beta, newton_max, cg_max, newton_e, newton_e2, cg_e, x = None, fs = None):
        """
        Run the barrier method on an input light field to solve for the volume it represents

        lf is the input light field
        ts is a finite sequence of the values of t for barrier method
        newton_alpha is the backtracking line search alpha
        newton_beta is the backtracking line search beta
        newton_max is the maximum number of Newton steps
        cg_max is the maximum number of CG steps
        newton_e is the epsilon for the Newton steps
        cg_e is the epsilon for the CG steps
        x is an initial (positive) guess and will be modified
        fs is the focal stack computed from lf and will be computed if not provided
        """
        self.log('b', 'ts = %s' % repr(ts))
        self.log('b', 'newton_alpha = %g, newton_beta = %g, newton_max = %g, cg_max = %g, newton_e = %g, cg_e = %g' % (newton_alpha, newton_beta, newton_max, cg_max, newton_e, cg_e))
        if x is None:
            const_level = 1.0 / math.sqrt(2*ts[0])
            self.log('b', 'Initializing input to barrier method to be at level %g' % const_level)
            x = self.const_volume(const_level)
        if fs is None:
            fs = self.blank_volume()
            self.compute_focal_stack(lf, fs)
        t = ts[0]
        tqueue = collections.deque(ts)
        while tqueue:
            last_t = t
            t = tqueue.popleft()
            self.log('b', 'Barrier method t = %g, running Newton' % t)
            if t == last_t:
                self.newton_recon(lf, fs, x, t, newton_alpha, newton_beta, newton_max, cg_max, newton_e, newton_e2, cg_e, throw_feas=False)
            else:
                try:
                    self.newton_recon(lf, fs, x, t, newton_alpha, newton_beta, newton_max, cg_max, newton_e, newton_e2, cg_e, throw_feas=True)
                except FeasibilityError, e:
                    # feasibility error, so backtrack state
                    backtrack = e.backtrack
                    self.log('b', 'Feasibility error detected, backtracking by %g' % backtrack)
                    new_t = math.pow(1.0*t/last_t, backtrack) * last_t
                    tqueue.appendleft(t)
                    tqueue.appendleft(new_t)
                    t = last_t
                    
        self.log('b', 'Barrier method finished')
        return x

    def barrier_recon_3(self, lf, t_start, mu, t_max, max_iters, newton_alpha, newton_beta, newton_max, cg_max, newton_e, newton_e2, cg_e, x = None, fs = None):
        """
        Run the barrier method on an input light field to solve for the volume it represents

        lf is the input light field
        ts is a finite sequence of the values of t for barrier method
        newton_alpha is the backtracking line search alpha
        newton_beta is the backtracking line search beta
        newton_max is the maximum number of Newton steps
        cg_max is the maximum number of CG steps
        newton_e is the epsilon for the Newton steps
        cg_e is the epsilon for the CG steps
        x is an initial (positive) guess and will be modified
        fs is the focal stack computed from lf and will be computed if not provided
        """
        self.log('b', 't_start = %g, mu = %g, t_max = %g, max_iters = %d' % (t_start, mu, t_max, max_iters))
        self.log('b', 'newton_alpha = %g, newton_beta = %g, newton_max = %g, cg_max = %g, newton_e = %g, cg_e = %g' % (newton_alpha, newton_beta, newton_max, cg_max, newton_e, cg_e))
        if x is None:
            const_level = 1.0 / math.sqrt(2*t_start)
            self.log('b', 'Initializing input to barrier method to be at level %g' % const_level)
            x = self.const_volume(const_level)
        if fs is None:
            fs = self.blank_volume()
            self.compute_focal_stack(lf, fs)
        t = t_start
        for i in range(max_iters):
            if t > t_max:
                break
            self.log('b', 'Barrier method t = %g, running Newton' % t)
            self.newton_recon(lf, fs, x, t, newton_alpha, newton_beta, newton_max, cg_max, newton_e, newton_e2, cg_e, throw_feas = False)
            t = t * mu
        self.log('b', 'Barrier method finished')
        return x
        
    def compute_merit_function(self, x, lf, t, temp_lf = None, temp_vol = None):
        "Compute the value of the merit function"
        if temp_lf is None:
            temp_lf = self.blank_lf()
        if temp_vol is None:
            temp_vol = self.blank_volume()
        self.compute_lf_projection(x, temp_lf)
        self.lf_mix(1.0, lf, -1.0, temp_lf, temp_lf)
        self.volume_log(x, temp_vol)
        return self.lf_dot(temp_lf, temp_lf) - sum(temp_vol)/t

    def compute_merit_function_2(self, x, lf, t):
        Pv = self.blank_lf()
        PTPv = self.blank_volume()
        PTl = self.blank_volume()
        logx = self.blank_volume()
        self.compute_lf_projection(x, Pv)
        self.compute_focal_stack(lf, PTl)
        self.compute_focal_stack(Pv, PTPv)
        self.volume_log(x, logx)
        return self.volume_dot(x, PTPv) - self.lf_dot(lf, Pv) - self.volume_dot(x, PTl) + self.lf_dot(lf, lf) - sum(logx)/t

    def compute_merit_function_3(self, x, lf, fs, t, focal_alg='psf'):
        PTPv = self.blank_volume()
        logx = self.blank_volume()
        if focal_alg=='psf':
            self.apply_focal_operator_psf(x, PTPv)
        elif focal_alg=='sipsf':
            self.apply_focal_operator_sipsf(x, PTPv)
        else:
            self.apply_focal_operator_lf(x, PTPv)
        self.volume_log(x, logx)
        return self.volume_dot(x, PTPv) - 2*self.volume_dot(x, fs) + self.lf_dot(lf, lf) - sum(logx)/t
           
    def newton_recon(self, lf, fs, x, t, newton_alpha, newton_beta, newton_max, cg_max, newton_e, newton_e2, cg_e, throw_feas=False):
        """
        Run Newton's method on an input light field with a fixed barrier term t and input x
        """
        goal = self.blank_volume()
        logx = self.blank_volume()
        x_nt = self.blank_volume()
        invx = self.blank_volume()
        y = self.blank_volume()
        y2 = self.blank_volume()
        gradient = self.blank_volume()
        lf_error = self.blank_lf()
        temp_x = self.blank_volume()
        ones = self.const_volume(1.0)
        for i in range(newton_max):
            self.log('n', 'Newton iteration %d' % i)
            self.log('n', 'x magnitude squared = %g' % self.volume_dot(x,x))
            self.log('n', 'writing x')
            open('newton_%g_%03d.tmp' % (t,i),'wb').write(self.dump_volume(x))
            self.log('n', 'Computing y(x)')
            # compute y(x)
            self.volume_reciprocal(x, invx)
            self.volume_scale(0.0, 1.0/t, invx, y)
            self.log('n', 'y magnitude squared = %g' % self.volume_dot(y,y))
            # compute y2(x)
            self.log('n', 'Computing y2(x)')
            self.volume_mult(invx, invx, y2)
            self.volume_scale(0.0, 1.0/t, y2, y2)
            self.log('n', 'y2 magnitude squared = %g' % self.volume_dot(y2,y2))
            self.log('n', 'y2 average value = %g' % (sum(y2)/(self.xlen*self.ylen*self.zlen)))

            # compute gradient
            focal_alg = self.focal_alg()

            if focal_alg == 'lf':
                self.log('n','Applying focal operator using LF...')
                self.apply_focal_operator_lf(x, gradient, lf_error)
                self.volume_mix(-1.0, fs, 1.0, gradient, gradient)
            elif focal_alg == 'psf':
                self.log('n','Applying focal operator using PSF...')
                self.apply_focal_operator_psf(x, gradient)
                self.volume_mix(-1.0, fs, 1.0, gradient, gradient)
            elif focal_alg == 'sipsf':
                self.log('n','Applying focal operator using SIPSF...')
                self.apply_focal_operator_sipsf(x, gradient)
                self.volume_mix(-1.0, fs, 1.0, gradient, gradient)
            else:
                self.log('n','Computing gradient (projection from current volume into lightfield)')
                self.compute_lf_projection(x, lf_error)
                self.log('n', 'light field magnitude squared: %g' % self.lf_dot(lf_error, lf_error))
                self.lf_mix(-1.0, lf, 1.0, lf_error, lf_error)
                self.log('n', 'light field error magnitude squared: %g' % self.lf_dot(lf_error, lf_error))
                self.log('n', 'Computing gradient (focal stack of light field error)')
                self.compute_focal_stack(lf_error, gradient)
            
            self.log('n', 'light field error focal stack magnitude squared: %g' % self.volume_dot(gradient, gradient))
            self.volume_mix(-1.0, y, 2.0, gradient, gradient)
            self.log('n', 'gradient energy = %g' % self.volume_dot(gradient, gradient))
            self.log('n', 'writing light field error')
            open('lferror_%g_%03d.tmp' % (t,i),'wb').write(self.dump_lf(lf_error))
            self.log('n', 'writing gradient')
            open('grad_%g_%03d.tmp' % (t,i),'wb').write(self.dump_volume(gradient))

            if self.volume_dot(gradient,gradient) <= 2*newton_e2:
                self.log('n', 'gradient energy low enough, terminating newton iteration')
                break

            # compute the newton step
            self.log('n', 'Running CG')
            x_nt = self.cg_recon(y2, gradient, cg_max, cg_e)
            self.volume_scale(0.0, -1.0, x_nt, x_nt)
            self.log('n', 'writing x_nt')
            open('nt_%g_%03d.tmp' % (t,i),'wb').write(self.dump_volume(x_nt))
            # compute decrement (lambda^2) from Newton step
            decrement = -self.volume_dot(x_nt, gradient)
            self.log('n', 'Newton decrement: %g' % decrement)
            # check stopping criterion
            if decrement <= 2*newton_e:
                self.log('n', 'decrement low enough, terminating newton iteration')
                break
            
            # BACKTRACKING LINE SEARCH
            self.log('n', 'Backtracking line search..')

            # compute the current value of the optimization function

            if focal_alg == 'lf':
                cur_val = self.compute_merit_function_3(x, lf, fs, t, focal_alg=focal_alg)
            elif focal_alg == 'psf':
                cur_val = self.compute_merit_function_3(x, lf, fs, t, focal_alg=focal_alg)
            elif focal_alg == 'sipsf':
                cur_val = self.compute_merit_function_3(x, lf, fs, t, focal_alg=focal_alg)
            else:
                self.volume_log(x, logx)
                cur_val = self.lf_dot(lf_error, lf_error) - self.volume_dot(logx, ones)/t

            
            self.log('n', 'current merit function value = %g' % cur_val)
            self.log('n', 'target a full step = %g' % (cur_val - newton_alpha * decrement))

            search_amt = 1.0

            # decrease search amount until we're feasible
            self.volume_mix(1.0, x, search_amt, x_nt, temp_x)
            feas_err = False
            while min(temp_x) < 0.0:
                feas_err = True
                search_amt = search_amt * newton_beta
                self.log('n', 'search_amt <- %g (feasibility)' % search_amt)
                self.volume_mix(1.0, x, search_amt, x_nt, temp_x)

            if feas_err and throw_feas:
                raise FeasibilityError(search_amt)

            # compute f(x+search_amt*x_nt)
            self.log('n', 'backtrack error calculation')
    
            self.volume_mix(1.0, x, search_amt, x_nt, temp_x)
            if focal_alg == 'lf':
                new_val = self.compute_merit_function_3(temp_x, lf, fs, t, focal_alg=focal_alg)
            elif focal_alg == 'psf':
                new_val = self.compute_merit_function_3(temp_x, lf, fs, t, focal_alg=focal_alg)
            elif focal_alg == 'sipsf':
                new_val = self.compute_merit_function_3(temp_x, lf, fs, t, focal_alg=focal_alg)
            else:
                new_val = self.compute_merit_function(temp_x, lf, t, temp_lf=lf_error, temp_vol=logx)
            self.log('n', 'new_val = %g' % new_val)
            # loop until condition is satisfied for backtracking line search
            while new_val >= cur_val - newton_alpha * search_amt * decrement:
                search_amt = search_amt * newton_beta
                self.log('n', 'search_amt <- %g (backtracking)' % search_amt)
                self.log('n', 'backtrack error calculation')
                self.volume_mix(1.0, x, search_amt, x_nt, temp_x)
                last_new_val = new_val
                if focal_alg == 'lf':
                    new_val = self.compute_merit_function_3(temp_x, lf, fs, t, focal_alg=focal_alg)
                elif focal_alg == 'psf':
                    new_val = self.compute_merit_function_3(temp_x, lf, fs, t, focal_alg=focal_alg)
                elif focal_alg == 'sipsf':
                    new_val = self.compute_merit_function_3(temp_x, lf, fs, t, focal_alg=focal_alg)
                else:
                    new_val = self.compute_merit_function(temp_x, lf, t, temp_lf=lf_error, temp_vol=logx)
                self.log('n', 'new_val = %g' % new_val)
                if new_val > last_new_val:
                    search_amt = search_amt / newton_beta
                    break
            
            # update (x = temp_x)
            self.volume_scale(0.0, 1.0, temp_x, x)
        # return answer
        return x

    def cg_recon(self, y2, gradient, cg_max, cg_e):
        """
        Run CG to calculate the Newton increment
        """
        w = self.blank_volume()
        temp_lf = self.blank_lf()
        temp_vol = self.blank_volume()
        x = self.blank_volume()
        r = self.blank_volume()
        p = self.blank_volume()
        Mr = self.blank_volume()

        f=open('krylov.m','ab')
        t = time.time()
        f.write('%% CG Iteration at %s\n' % (time.strftime('%Y-%m-%d %H:%M:%S.', time.localtime(t)) + ('%f' % (t-math.floor(t)))[2:]))
        f.flush()

        # calculate approximation to Jacobi preconditioner
        M = self.blank_volume()
        self.volume_scale(2.0, 1.0, y2, M)
        self.volume_reciprocal(M, M)

        self.log('c', 'preconditioning min: %g, max %g, condition %g' % (min(M),max(M),max(M)/min(M)))

        # x_0 = 0
        x_coeffs = P([])
        # r_0 = b
        self.volume_scale(0.0, 1.0, gradient, r)
        r_coeffs = P([1.0])

        self.log('C', 'x = %s' % str(x_coeffs))
        self.log('c', 'r = %s' % str(r_coeffs))

        f.write('r{1} = [' + ','.join([('%g' % r_coeffs[i]) for i in range(len(r_coeffs))]) + '];\n')
        f.flush()

        # Mr = Mr
        self.volume_mult(M, r, Mr)

        rho = self.volume_dot(r,Mr)
        bnorm = rho
        last_rho = 0.0
        for i in range(cg_max):
            self.log('c', 'cg iter %d' % i)
            self.log('c', 'rho = %g' % rho)
            self.log('c', 'bnorm = %g' % bnorm)
            if rho <= cg_e * cg_e * bnorm:
                break
            if i:
                # p = Mr + (rho/last_rho) * p
                self.volume_mix(1.0, Mr, rho/last_rho, p, p)
                p_coeffs = r_coeffs + (rho/last_rho)*p_coeffs
            else:
                # p_0 = Mr
                self.volume_mult(M, r, p)
                p_coeffs = P(r_coeffs)
            self.log('C', 'p = %s' % str(p_coeffs))
            # w := Ap
            focal_alg = self.focal_alg()
            if focal_alg == 'psf':
                self.apply_focal_operator_psf(p,w)
            elif focal_alg == 'sipsf':
                self.apply_focal_operator_sipsf(p,w)
            elif focal_alg == 'lf':
                self.apply_focal_operator_lf(p,w,temp_lf)
            else:
                self.compute_lf_projection(p, temp_lf)
                self.compute_focal_stack(temp_lf, w)
            self.volume_mult(y2, p, temp_vol)
            self.volume_mix(1.0, temp_vol, 2.0, w, w)
            w_coeffs = p_coeffs << 1
            # alpha := rho / (p^Tw)
            alpha = rho / self.volume_dot(p, w)
            # x := x + alpha * p
            self.volume_mix(1.0, x, alpha, p, x)
            x_coeffs = x_coeffs + alpha*p_coeffs
            # r := r - alpha * w
            self.volume_mix(1.0, r, -alpha, w, r)
            self.log('c', 'alpha = %g' % alpha)
            r_coeffs = r_coeffs - alpha*w_coeffs
            f.write('r{'+str(i+2)+'} = [' + ','.join([('%g' % r_coeffs[i]) for i in range(len(r_coeffs))]) + '];\n')
            f.flush()
            last_rho = rho
            # Mr = Mr
            self.volume_mult(M, r, Mr)
            # rho_k := r^T Mr
            rho = self.volume_dot(r,Mr)
            self.log('C', 'x = %s' % str(x_coeffs))
            self.log('c', 'r = %s' % str(r_coeffs))
        f.close()
        return x
