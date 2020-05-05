import time

iteration_number = 3
b = 2

def fun_c():
    import ctypes
    cdll = ctypes.cdll.LoadLibrary(r"D:\SourceCode\FlashDTA_sourcecode\test\STALite.dll")
    time_start = time.time()

    network_compu = cdll.network_assignment
    # add_fun.argtypes=[ctypes.c_float,ctypes.c_float]
    network_compu.restype = ctypes.c_double
    c = network_compu(iteration_number,b)


    time_end = time.time()
    print('STALite')
    print('output: {}, total time: {}'.format(c, time_end-time_start))


if __name__ == "__main__":
    fun_c()
    







