import array as binarray

def read_array(fin, astype='B', num=1):
    '''
        Reads an element from a binary stream.
    '''
    var = binarray.array(astype)
    var.fromfile(fin, num)
    if num == 1:
        return var[0]
    return var.tolist()

def read_byte(fin, num=1):
    '''
        Reads bytes from a binary stream.  Returns a byte stream.
    '''
    var = binarray.array('B')
    var.fromfile(fin, num)
    return var.tostring()

def read_integer(fin, num=1):
    '''
        Reads integers from a binary stream
    '''
    return read_array(fin, astype='i', num=num)

def read_float(fin, num=1):
    '''
        Reads floating point values from a binary stream
    '''
    return read_array(fin, astype='f', num=num)

def read_char(fin, num=1):
    '''
        Reads signed char values from a binary stream
    '''
    return read_array(fin, astype='b', num=num)

def read_shortint(fin, num=1):
    '''
        Reads halfword signed integer values from a binary stream
    '''
    return read_array(fin, astype='h', num=num)

def read_shortuint(fin, num=1):
    '''
        Reads halfword unsigned integer values from a binary stream
    '''
    return read_array(fin, astype='H', num=num)

def read_cints(fin, num=1):
    '''
        Reads integer complex values from a binary stream
    '''
    values = read_array(fin, astype='i', num=2*num)
    return [complex(x[0],x[1]) for x in zip(*[values[::2],values[1::2]])]

def read_fints(fin, num=1):
    '''
        Reads floating point complex values from a binary stream
    '''
    values = read_array(fin, astype='f', num=2*num)
    return [complex(x[0],x[1]) for x in zip(*[values[::2],values[1::2]])]

def write_array(fout, values, astype='B'):
    '''
        Writes values to a binary stream.
    '''
    var = binarray.array(astype)
    try:
        var.extend(values)
    except TypeError:
        var.append(values)
    var.tofile(fout)

def write_integer(fout, values):
    '''
        Writes integer values to a binary stream.
    '''
    write_array(fout, values, 'i')

def write_float(fout, values):
    '''
        Writes floating point values to a binary stream.
    '''
    write_array(fout, values, 'f')

def write_byte(fout, values):
    '''
        Writes bytes to a binary stream.
    '''
    write_array(fout, values, 'B')

def write_char(fout, values):
    '''
        Writes signed char values to a binary stream.
    '''
    write_array(fout, values, 'b')

def write_shortint(fout, values):
    '''
        Writes halfword unsigned integer values to a binary stream.
    '''
    write_array(fout, values, 'h')

def write_shortuint(fout, values):
    '''
        Writes halfword unsigned integer values to a binary stream.
    '''
    write_array(fout, values, 'H')

def write_cint(fout, values):
    '''
        Writes integer complex values to a binary stream.
    '''
    for value in values:
        write_integer(fout, (value.real, value.imag))

def write_fint(fout, values):
    '''
        Writes floating point complex values to a binary stream.
    '''
    for value in values:
        write_float(fout, (value.real, value.imag))
