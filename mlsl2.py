"""
Multithreaded extensions to mlsl
"""

import mlsl

import thread
import threading
import array
import struct
import random
import math

class Error(Exception):
    pass

def grayify(img, img_type, src_col_offset, src_row_offset,
            src_col_stride, src_row_stride, src_chan_stride,
            dst_col_offset, dst_row_offset,
            dst_col_stride, dst_row_stride, cols, rows, chans,
            weights, gamma, output=None, threads=2):

    if threads < 1:
        raise Error('the argument threads must be positive')
    if threads < 2:
        return mlsl.grayify(img, img_type, src_col_offset, src_row_offset,
                            src_col_stride, src_row_stride, src_chan_stride,
                            dst_col_offset, dst_row_offset,
                            dst_col_stride, dst_row_stride,
                            cols, rows, chans,
                            weights, gamma, output)

    if threads > rows:
        threads = rows

    # create output buffer
    output_is_array = False
    if output is None:
        output_size = (struct.calcsize('f') +
                       dst_col_stride*(dst_col_offset+cols-1) +
                       dst_row_stride*(dst_row_offset+rows-1))
        output = array.array('f', output_size*'\x00')
        output_is_array = True

    done = threading.Semaphore(0)

    def thread_func(start, end):
        try:
            mlsl.grayify(img, img_type, src_col_offset,
                         src_row_offset + start,
                         src_col_stride, src_row_stride,
                         src_chan_stride, dst_col_offset,
                         dst_row_offset + start,
                         dst_col_stride, dst_row_stride,
                         cols, end-start, chans,
                         weights, gamma, output)
        finally:
            done.release()

    starts = [int(round(1.0*x*rows/threads)) for x in range(threads)]
    starts.append(cols)

    for j in range(len(starts)-1):
        start = starts[j]
        end = starts[j+1]
        thread.start_new_thread(thread_func, (start,end))

    for j in range(len(starts)-1):
        done.acquire()

    if output_is_array:
        output = output.tostring()

    return output

def correlate2s(src, src_img_type, src_byte_offset,
               src_col_offset, src_row_offset, src_col_stride,
               src_row_stride, src_col_total, src_row_total,
               dst_img_type, dst_byte_offset, dst_col_offset,
               dst_row_offset, dst_col_stride, dst_row_stride,
               dst_col_total, dst_row_total, cols, rows, kernel_h,
               kernel_h_offset, kernel_v, kernel_v_offset, min_val=0,
               max_val=0, output=None, threads=2):

    if threads < 1:
        raise Error('the argument threads must be positive')
    if threads < 2:
        return mlsl.correlate2s(src, src_img_type, src_byte_offset,
                               src_col_offset, src_row_offset, src_col_stride,
                               src_row_stride, src_col_total, src_row_total,
                               dst_img_type, dst_byte_offset, dst_col_offset,
                               dst_row_offset, dst_col_stride, dst_row_stride,
                               dst_col_total, dst_row_total, cols, rows, kernel_h,
                               kernel_h_offset, kernel_v, kernel_v_offset, min_val,
                               max_val, output)

    if threads > min(cols,rows):
        threads = min(cols,rows)

    # create intermediate buffer
    temp_col_stride = struct.calcsize('f')
    temp_row_stride = temp_col_stride * (cols|1)
    temp_size = temp_row_stride*rows

    temp = array.array('f', '\x00'*temp_size)

    # create output buffer if there isn't one

    output_size = (dst_byte_offset
                   + max(dst_col_stride*(dst_col_offset+cols-1),0)
                   + max(dst_row_stride*(dst_row_offset+rows-1),0)
                   + struct.calcsize(dst_img_type))

    output_is_array = False
    if output is None:
        output = array.array(dst_img_type, '\x00'*output_size)
        output_is_array = True

    done_h = threading.Semaphore(0)
    done_v = threading.Semaphore(0)

    def horiz_thread(start, end):
        "Thread for horizontal filtering"
        try:
            mlsl.correlate2s(src, src_img_type, src_byte_offset,
                            src_col_offset, src_row_offset+start,
                            src_col_stride, src_row_stride,
                            src_col_total, src_row_total,
                            'f', 0,
                            0, start,
                            temp_col_stride, temp_row_stride,
                            cols, rows,
                            cols, end-start,
                            kernel_h, kernel_h_offset, None, 0,
                            0.0, 0.0, output=temp)
        finally:
            done_h.release()

    def vert_thread(start, end):
        "Thread for vertical filtering"
        try:
            mlsl.correlate2s(temp, 'f', 0,
                            start, 0,
                            temp_col_stride, temp_row_stride,
                            cols, rows,
                            dst_img_type, dst_byte_offset,
                            dst_col_offset+start, dst_row_offset,
                            dst_col_stride, dst_row_stride,
                            dst_col_total, dst_row_total,
                            end-start, rows,
                            None, 0, kernel_v, kernel_v_offset,
                            min_val, max_val, output=output)
        finally:
            done_v.release()

    # horizontal filtering
    for j in range(threads):
        start = int(round(1.0*j*rows/threads))
        end = int(round(1.0*(j+1)*rows/threads))
        thread.start_new_thread(horiz_thread, (start,end))

    for j in range(threads):
        done_h.acquire()

    # vertical filtering
    for j in range(threads):
        start = int(round(1.0*j*cols/threads))
        end = int(round(1.0*(j+1)*cols/threads))
        thread.start_new_thread(vert_thread, (start,end))

    for j in range(threads):
        done_v.acquire()

    # return the result
    if output_is_array:
        return output.tostring()
    else:
        return output

def correlate2(src, src_img_type, src_byte_offset,
              src_col_offset, src_row_offset,
              src_col_stride, src_row_stride,
              src_col_total, src_row_total,
              dst_img_type, dst_byte_offset,
              dst_col_offset, dst_row_offset,
              dst_col_stride, dst_row_stride,
              dst_col_total, dst_row_total,
              cols, rows,
              kernel, kernel_h_offset,
              kernel_width, kernel_v_offset,
              min_val=0, max_val=0, output=None, threads=2):
    if threads < 1:
        raise Error('the argument threads must be positive')
    if threads < 2:
        return mlsl.correlate2(src, src_img_type, src_byte_offset,
                              src_col_offset, src_row_offset, src_col_stride,
                              src_row_stride, src_col_total, src_row_total,
                              dst_img_type, dst_byte_offset, dst_col_offset,
                              dst_row_offset, dst_col_stride, dst_row_stride,
                              dst_col_total, dst_row_total, cols, rows,
                              kernel, kernel_h_offset, kernel_width,
                              kernel_v_offset, min_val, max_val, output=output)

    if threads > rows:
        threads = rows

    # create output buffer if there isn't one

    output_size = (dst_byte_offset
                   + max(dst_col_stride*(dst_col_offset+cols-1),0)
                   + max(dst_row_stride*(dst_row_offset+rows-1),0)
                   + struct.calcsize(dst_img_type))

    output_is_array = False
    if output is None:
        output = array.array(dst_img_type, '\x00'*output_size)
        output_is_array = True

    done_sema = threading.Semaphore(0)

    def thread_func(start, end):
        try:
            mlsl.correlate2(src, src_img_type, src_byte_offset,
                           src_col_offset, src_row_offset+start, src_col_stride,
                           src_row_stride, src_col_total, src_row_total,
                           dst_img_type, dst_byte_offset, dst_col_offset,
                           dst_row_offset+start, dst_col_stride, dst_row_stride,
                           dst_col_total, dst_row_total, cols, end-start,
                           kernel, kernel_h_offset, kernel_width,
                           kernel_v_offset, min_val, max_val, output)
        finally:
            done_sema.release()
    
    for i in range(threads):
        start = int(round(1.0*i/threads*rows))
        end = int(round(1.0*(i+1)/threads*rows))
        thread.start_new_thread(thread_func, (start,end))

    for i in range(threads):
        done_sema.acquire()

    if output_is_array:
        return output.tostring()
    else:
        return output


def maxima(img, img_type, col_offset, row_offset,
           col_stride, row_stride, col_total, row_total,
           cols, rows, distance, threads=2):

    if threads < 1:
        raise Error('the argument threads must be positive')
    if threads < 2:
        return mlsl.maxima(img, img_type, col_offset, row_offset,
                           col_stride, row_stride, col_total,
                           row_total, cols, rows, distance)


    output = []
    lock = threading.RLock()
    done = threading.Semaphore(0)

    if threads > rows:
        threads = rows;

    def thread_func(start, end):
        try:
            local_maxima = mlsl.maxima(img, img_type, col_offset, start,
                                       col_stride, row_stride, col_total,
                                       row_total, cols, end-start, distance)
            lock.acquire()
            try:
                output.extend(local_maxima)
            finally:
                lock.release()
        finally:
            done.release()

    for j in range(threads):
        start = row_offset + int(round(1.0*j*rows/threads))
        end = row_offset + int(round(1.0*(j+1)*rows/threads))
        thread.start_new_thread(thread_func, (start,end))

    for j in range(threads):
        done.acquire()

    return output

def edges(vertices, min_length, max_length,
          start_vertex=0, end_vertex=None, threads=2):

    if end_vertex is None:
        end_vertex = len(vertices)

    if threads < 1:
        raise Error('the argument threads must be positive')
    if threads < 2:
        return mlsl.edges(vertices, min_length, max_length,
                          start_vertex, end_vertex)

    output = []
    lock = threading.RLock()
    done = threading.Semaphore(0)

    def thread_func(start, end):
        try:
            out_edges = mlsl.edges(vertices, min_length, max_length,
                                   start, end)
            lock.acquire()
            try:
                output.extend(out_edges)
            finally:
                lock.release()
        finally:
            done.release()

    starts = []
    ends = []

    begin_edges = start_vertex - 1
    end_edges = end_vertex - 2
    total_edges = (end_vertex-start_vertex) * (end_vertex+start_vertex-3) / 2
    edges_per_thread = total_edges / threads + 1

    cur_start = start_vertex
    for j in range(threads):
        if cur_start >= end_vertex:
            break
        bucket = 0
        for end in range(cur_start,end_vertex+1):
            if bucket >= edges_per_thread:
                starts.append(cur_start)
                ends.append(end)
                cur_start = end
                break
            bucket += (end-1)
        starts.append(cur_start)
        ends.append(end_vertex)
        break

    if ends[-1] != end_vertex:
        raise Error("Unable to allocate jobs")

    for j in range(len(starts)):
        thread.start_new_thread(thread_func, (starts[j],ends[j]))

    for j in range(len(starts)):
        done.acquire()

    return output

def latticise2(vertices, edges, seed_vertex, seed_i, seed_j):
    result = mlsl.latticise(vertices, edges, seed_vertex, seed_i, seed_j, True)
    if not result[2]:
        result = mlsl.latticise(vertices, edges, seed_vertex, seed_i, seed_j, False)
    return result

def latticise(vertices, edges, seed_vertex=None, seed_i=0, seed_j=0, axis_aligned=None, max_iters=-1):

    if seed_vertex is not None:
        return latticise2(vertices, edges, seed_vertex, seed_i, seed_j)

    best_result = None

    unlabeled = set([i for i in range(len(vertices)) if len(vertices[i]) == 2])

    maxlabel = 0

    iters = max_iters

    while unlabeled and iters != 0:

        # choose a random vertex to seed

        seed_vertex = random.choice(tuple(unlabeled))

        # run the indexing operation

        result = latticise2(vertices, edges, seed_vertex, seed_i, seed_j)

        # figure out how many vertices were labeled

        output_vertices = result[0]
        labeled = set([i for i in range(len(output_vertices)) if len(output_vertices[i]) > 2])

        # keep best one
        
        if len(labeled) > maxlabel:
            maxlabel = len(labeled)
            best_result = result

        # remove labeled vertices from unlabeled set

        prev_len = len(unlabeled)

        unlabeled.difference_update(labeled)

        # stop if we've found something good enough or if we can't possibly find anything better
        
        next_len = len(unlabeled)
        if next_len < maxlabel or maxlabel > len(vertices)/2 or next_len == prev_len:
            break

        iters -= 1

    return best_result

def lf_gaussians(img, img_type, x_stride, y_stride, u_stride, v_stride,
                 x_len, y_len, u_len, v_len,
                 ulens_pitch, ulens_flen, ulens_n,
                 spec_mag, spec_na, spec_n, min_value, min_corr=0.25,
                 u_center=None, v_center=None,
                 u_width=None, v_width=None, threads=2):

    # get geometry of scene
    
    (pitch, slopes) = mlsl.geometry(u_len, v_len, ulens_pitch, ulens_flen, ulens_n,
                                    spec_mag, spec_na, spec_n, u_center, v_center,
                                    u_width, v_width)

    # get a list of useful u,v coordinates

    uv_list = slopes.keys()

    uv_number = len(uv_list)

    # get center
    x_center = 0.5*(x_len-1)
    y_center = 0.5*(y_len-1)

    # output

    output = []
    lock = threading.RLock()
    done = threading.Semaphore(0)

    # set up the threads

    def thread_func(start, end):
        thread_output = []
        for i in range(start,end):
            try:
                u,v = uv_list[i]
                img_offset = u_stride*u + v_stride*v
                u_f, v_f, u_s2, v_s2 = slopes[(u,v)]
                result = mlsl.gaussians(img, img_type, x_stride, y_stride,
                                        0, 0, x_len, y_len, x_len, y_len,
                                        min_value, min_corr, img_offset)
                for (x_i, y_i, x_f, y_f, x_s2, y_s2, x_r, y_r) in result:
                    x_f -= x_center
                    y_f -= y_center
                    x_f *= pitch
                    y_f *= -pitch
                    x_s2 *= (pitch*pitch)
                    y_s2 *= (pitch*pitch)
                    thread_output.append((x_i,y_i,u,v,
                                          x_f,y_f,u_f,-v_f,
                                          x_s2,y_s2,u_s2,v_s2))
            except Exception, e:
                print e

        # append the output to the main result
        lock.acquire()
        output.extend(thread_output)
        lock.release()
        done.release()
    
    # start the threads
    threads = min(threads, uv_number)
    for i in range(threads):
        start = int(round(1.0*i*uv_number/threads))
        end = int(round(1.0*(i+1)*uv_number/threads))
        thread.start_new_thread(thread_func,(start,end))

    # wait for threads to finish
    for i in range(threads):
        done.acquire()

    return output
    
    
def lf_points(rays, points=None, min_bundle_size=50, max_distance=1.0, iterations=10,
              threads=2):
    """
    Given a set of rays, find likely points that have generated these rays

    rays must be either a list of tuples
        (x_i,y_i,u_i,v_i,x_f,y_f,u_f,v_f,x_s2,y_s2,u_s2,v_s2)
        or a list of tuples
        (x,y,u,v,x_s2,y_s2,u_s2,v_s2)
    """
    if len(rays[0]) == 12:
        rays = [(e,f,g,h,i,j,k,l) for (a,b,c,d,e,f,g,h,i,j,k,l) in rays]
    elif len(rays[0]) != 8:
        raise Error('input rays not specified correctly')

    # get a list of candidate points if needed
    if not points:
        points, remainder = mlsl.extract_points(rays, 50, 2.0*max_distance)

    # quit if we don't have guess
    if not points:
        return ([], [], rays[:])

    # split up to refine
    num_points = len(points)
    if threads > num_points: threads = num_points

    result_points = []
    result_bundles = []
    result_lock = threading.RLock()
    result_done = threading.Semaphore(0)

    def thread_func(input):
        try:
            my_points = []
            my_bundles = []
            for (x,y,z) in input:
                bundle = mlsl.find_rays(rays, (x,y,z), max_distance)
                for i in range(iterations):
                    x,y,z = mlsl.triangulate_simple(bundle)
                    bundle = mlsl.find_rays(rays, (x,y,z), max_distance)
                if not bundle:
                    continue
                my_points.append((x,y,z))
                my_bundles.append(bundle)

            result_lock.acquire()
            try:
                for point,bundle in zip(my_points,my_bundles):
                    if point not in result_points:
                        result_points.append(point)
                        result_bundles.append(bundle)
            finally:
                result_lock.release()
        finally:
            result_done.release()

    for i in range(threads):
        start_index = int(round(1.0*i*num_points/threads))
        end_index = int(round(1.0*(i+1)*num_points/threads))
        thread.start_new_thread(thread_func, (points[start_index:end_index],))

    for i in range(threads):
        result_done.acquire()

    def found(elem,contlist):
        for cont in contlist:
            if elem in cont:
                return True
        return False

    remainder = [x for x in rays if not found(x,result_bundles)]

    return result_points, result_bundles, remainder
    
def lf_to_lenslets(lf_data, lf_align, lf_format, lf_layout,
                   coffset, clen, xlen, ylen, ulen, vlen,
                   img_type, byte_offset, chan_stride, col_stride, row_stride,
                   chan_offset, col_offset, row_offset,
                   chans, cols, rows,
                   chan_total, col_total, row_total,
                   rtm, rect, image=None, threads=2):
    # create output
    if image is None:
        img_size = byte_offset + max(0,row_stride * (row_total-1)) + max(col_stride * (col_total-1),0) + max(chan_stride * (chan_total-1),0) + struct.calcsize(img_type)
        image = array.array(img_type, '\x00'*img_size)
        convert_back = True
    else:
        convert_back = False

    # create done signal
    result_done = threading.Semaphore(0)

    # create thread function
    def thread_func(row_start, row_end):
        try:
            mlsl.lf_to_lenslets(lf_data, lf_align, lf_format, lf_layout,
                                coffset, clen, xlen, ylen, ulen, vlen,
                                img_type, byte_offset, chan_stride, col_stride, row_stride,
                                chan_offset, col_offset, row_start,
                                chans, cols, row_end-row_start,
                                chan_total, col_total, row_total,
                                rtm, rect, image)
        finally:
            result_done.release()

    # assign threads
    if threads > rows:
        threads = rows

    for i in range(threads):
        row_start = row_offset + int(round(1.0*i*rows/threads ))
        row_end = row_offset + int(round(1.0*(i+1)*rows/threads ))
        thread.start_new_thread(thread_func, (row_start, row_end))

    # wait for threads to finish
    for i in range(threads):
        result_done.acquire()

    # done
    if convert_back:
        return image.tostring()
    else:
        return image

def lf_blur3d(in_data, in_align, in_format, in_layout,
              out_align, out_format, out_layout,
              clen, xlen, ylen, ulen, vlen,
              var_x, var_y, var_z, pitch, slopes, out_data=None,
              threads=2):
    """
    Blur a light field such that the generated focal stack is blurred by
    the given variances (in microns^2).  The light field geometry
    is given by pitch and slopes (output from mlsl.geometry)

    """

    # calculate the light field strides
    (byte_offset,
     c_stride,
     x_stride, y_stride,
     u_stride, v_stride,
     full_size) = mlsl.strides(in_format, in_layout, in_align,
                               clen, xlen, ylen, ulen, vlen)

    (out_byte_offset,
     out_c_stride,
     out_x_stride, out_y_stride,
     out_u_stride, out_v_stride,
     out_full_size) = mlsl.strides(out_format, out_layout,
                                   out_align,
                                   clen, xlen, ylen, ulen, vlen)

    out_pixel_size = struct.calcsize(out_format)

    if out_full_size % out_pixel_size:
        out_full_size += out_pixel_size - (out_full_size % out_pixel_size)

    # allocate an output buffer if it's not there
    make_string = False
    if out_data is None:
        make_string = True
        out_data = array.array(out_format, '\x00'*out_full_size)

    all_jobs = []

    # do some precalculations
    var_x_norm = var_x / (2*math.pi)
    var_y_norm = var_y / (2*math.pi)
    var_z_norm = var_z / (2*math.pi)

    for (u,v) in slopes:
        # get the actual slope for a u,v
        u_f, v_f, u_s2, v_s2 = slopes[(u,v)]

        print u_f, v_f

        # calculate new covariance matrix
        var_h_norm = var_x_norm + var_z_norm * u_f * u_f
        var_v_norm = var_y_norm + var_z_norm * v_f * v_f
        var_c_norm = u_f * v_f * var_z_norm

        # rescale variances
        var_h = 2*math.pi * var_h_norm / (pitch*pitch)
        var_v = 2*math.pi * var_v_norm / (pitch*pitch)
        var_c = 2*math.pi * var_c_norm / (pitch*pitch)

        print var_h, var_v, var_c

        # create Gaussian kernel
        kernel = mlsl.gaussian2(var_h, var_v, var_c)

        # iterate over color channels
        for c in range(clen):
            src_byte_offset = byte_offset + c*c_stride + u*u_stride + v*v_stride
            dst_byte_offset = out_byte_offset + c*out_c_stride + u*out_u_stride + v*out_v_stride
            # append to job list
            all_jobs.append((src_byte_offset, dst_byte_offset, kernel))

    if len(all_jobs) < threads:
        threads = len(all_jobs)

    # threads

    done = threading.Semaphore(0)

    def work_thread(jobs):
        try:
            for (src_byte_offset, dst_byte_offset, kernel) in jobs:
                k, k_h_offset, k_width, k_v_offset = kernel
                mlsl.correlate2(in_data, in_format, src_byte_offset,
                               0, 0, x_stride, y_stride,
                               xlen, ylen,
                               out_format, dst_byte_offset,
                               0, 0, out_x_stride, out_y_stride,
                               xlen, ylen,
                               xlen, ylen,
                               k, k_h_offset, k_width, k_v_offset, output=out_data)
        finally:
            done.release()

    for i in range(threads):
        job_start = int(round(1.0*len(all_jobs)*i/threads))
        job_end = int(round(1.0*len(all_jobs)*(i+1)/threads))
        job_segment = all_jobs[job_start:job_end]
        thread.start_new_thread(work_thread,(job_segment,))

    for i in range(threads):
        done.acquire()

    if make_string:
        out_data = out_data.tostring()

    return out_data

        
def lf_blur_subaperture(in_data, in_align, in_format, in_layout,
                        clen, xlen, ylen, ulen, vlen, u, v,
                        k, k_h_offset, k_width, k_v_offset,
                        dst_img_type, dst_byte_offset,
                        dst_col_stride, dst_row_stride,
                        output=None, threads=2):
    """
    Return a subaperture image (u,v) from given light field correlated
    with the given kernel
    """

    # calculate light field strides
    (lf_byte_offset,
     c_stride,
     x_stride, y_stride,
     u_stride, v_stride,
     full_size) = mlsl.strides(in_format, in_layout, in_align,
                               clen, xlen, ylen, ulen, vlen)

    if u < 0 or u >= ulen or v < 0 or v >= vlen:
        raise Error('u or v out of bounds')

    lf_byte_offset += u_stride * u + v_stride * v
    
    if threads < 1:
        raise Error('the argument threads must be positive')
    if threads < 2:
        return mlsl.correlate2(in_data, in_format, lf_byte_offset,
                              0, 0, x_stride, y_stride,
                              xlen, ylen,
                              dst_img_type, dst_byte_offset,
                              0, 0, dst_col_stride, dst_row_stride,
                              xlen, ylen,
                              xlen, ylen,
                              k, k_h_offset, k_width, k_v_offset, output=output)

    if threads > ylen:
        threads = ylen

    # create a shared output buffer if needed

    if output is None:
        max_col_reach = max(0, (dst_col_stride-1)*xlen)
        max_row_reach = max(0, (dst_row_stride-1)*ylen)
        max_reach = dst_byte_offset + max_col_reach + max_row_reach + struct.calcsize(dst_img_type)
        output = array.array(dst_img_type, '\x00'*max_reach)
        is_array = True
    else:
        is_array = False

    # create semaphores for task completion
    completed = threading.Semaphore(0)

    def thread_func(start,end):
        try:
            mlsl.correlate2(in_data, in_format, lf_byte_offset,
                           0, start, x_stride, y_stride,
                           xlen, ylen,
                           dst_img_type, dst_byte_offset,
                           0, start, dst_col_stride, dst_row_stride,
                           xlen, ylen,
                           xlen, end-start,
                           k, k_h_offset, k_width, k_v_offset, output=output)
        finally:
            completed.release()

    # start tasks

    for i in range(threads):
        start = int(round(1.0*i/threads*ylen))
        end = int(round(1.0*(i+1)/threads*ylen))
        thread.start_new_thread(thread_func, (start,end))

    # wait for task completion

    for i in range(threads):
        completed.acquire()

    if is_array:
        return output.tostring()
    else:
        return output

