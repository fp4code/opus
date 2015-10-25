# -*- coding: utf-8 -*-
# (C) Fabrice Pardo, MiNaO, LPN, CNRS: version 0.0.2

import struct
import numpy as np

OPUS_NAMES = {
    'ABP':'Absolute Peak Pos in Laser*2',
    'ACC':'Accessory',
    'AG2':'Actual Signal Gain 2nd Channel',
    'AN1':'Analog Signal 1',
    'AN2':'Analog Signal 2',
    'APF':'Apodization Function',
    'APT':'Aperture Setting',
    'AQM':'Acquisition Mode',
    'ARG':'Actual Signal Gain',
    'ARS':'Number of Background Scans',
    'ASG':'Actual Ref. Signal Gain',
    'ASS':'Actual Signal Gain',
    'BBW':'Number of Bad BW Scans',
    'BFW':'Number of Bad FW Scans',
    'BLD':'Building',    
    'BMS':'Beamsplitter Setting',
    'CHN':'Measurement Channel',
    'CNM':'Operator Name',
    'COR':'Correlation Test Mode',
    'CPY':'Company',
    'CSF':'Y - Scaling Factor',
    'DAQ':'Data Aquisition Status',
    'DAT':'Date of Measurement',
    'DEL':'Delay Before Measurement',
    'DLY':'Stabilization Delay',
    'DPF':'Data Point Format',
    'DPM':'Department',
    'DTC':'Detector Setting',
    'DUR':'Scan time (sec)',
    'DXU':'X Units',
    'EXP':'Experiment',
    'FOC':'Focal Length',
    'FXV':'Frequency of First Point',
    'GBW':'Number of Good BW Scans',
    'GFW':'Number of Good FW Scans',
    'HFL':'High Folding Limit',
    'HFQ':'End Frequency Limit for File',
    'HFW':'Wanted High Frequency Limit',
    'HPF':'High Pass Filter',
    'HUM':'Humidity Interferometer',
    'INS':'Instrument Type',
    'IST':'Instrument Status',
    'LCT':'Location',
    'LFL':'Low Folding Limit',
    'LFQ':'Start Frequency Limit for File',
    'LFW':'Wanted Low Frequency Limit',
    'LPF':'Low Pass Filter',
    'LWN':'Laser Wavenumber',
    'LXV':'Frequency of Last Point',
    'MVD':'Max. Velocity Deviation',
    'MXY':'Y - Maximum',
    'MNY':'Y - Minimum',
    'NLI':'Non Linearity Correction',
    'NPT':'Number of Data Points',
    'NSN':'Scan Number',
    'NSS':'Sample Scans',
    'OPF':'Optical Filter Setting',
    'P2A':'Peak Amplitude 2nd Channel',
    'P2K':'Backward Peak Location 2nd Channel',
    'P2L':'Peak Location 2nd Channel',
    'P2R':'Backward Peak Amplitude 2nd Channel',
    'PGN':'Preamplifier Gain',
    'PHR':'Phase Resolution',
    'PHZ':'Phase Correction Mode',
    'PKA':'Peak Amplitude',
    'PKL':'Peak Location',    
    'PLF':'Result Spectrum',
    'PRA':'Backward Peak Amplitude',
    'PRL':'Backward Peak Location',
    'PRS':'Pressure Interferometer (hPa)',
    'RCH':'Background Measurement Channel',
    'RDX':'Extended Ready Check',
    'RDY':'Ready Check',    
    'RES':'Resolution',
    'RG2':'Signal Gain, Background 2nd Channel',
    'RGN':'Signal Gain, Background',
    'RSN':'Running Sample Number',
    'SFM':'Sample Form',
    'SG2':'Signal Gain, Sample 2nd Channel',
    'SGN':'Signal Gain, Sample',
    'SNM':'Sample Name',
    'SON':'External Synchronisation',
    'SPZ':'Stored Phase Mode',
    'SRC':'Source Setting',
    'SRN':'Instrument Serial Number',    
    'SRT':'Start time (sec)',
    'SSM':'Sample Spacing Multiplicator',
    'SSP':'Sample Spacing Divisor',
    'TDL':'To do list',
    'TPX':'Total Points X',
    'TSC':'Scanner Temperature',
    'TIM':'Time of Measurement',   
    'VEL':'Scanner Velocity',    
    'VSN':'Firmware version',
    'XPP':'Experiment Path',
    'ZFF':'Zero Filling Factor',    
}

OPUS_TYPES = {
    'B3':'Blackman-Harris 3-Term',
    'DD':'Double Sided,Forward-Backward',
    'ML':'Mertz',
    'NO':'No',
    'PNT':'Points',
    'RFL':'Reflectance',
    'WN':'Wavenumber cm-1',
    '-1':'Automatic',
    '0':'Open or OFF',
    '1':'ON',
    '2':'Open',
    '3':'3',
    '8':'8',
    '16':'16',
    '1.6':'1.60',
    '2.5':'2.50',
    '5.0':'5.00',
}

OPUS_TYPES_SPECIFIC = {
    'HPF,0':'Open',
    'RDX,0':'OFF',
}

OPUS_BLOCS = {
    (0, 0, 0, 0):    'garbage blocs',
    (0, 0, 104, 64): 'Datafile History',
    (7, 4, 0, 0):    'ScSm',
    (23, 4, 0, 0):   'Data Parameters ScSm',
    (7, 8, 0, 0):    'IgSm',
    (23, 8, 0, 0):   'Data Parameters IgSm',
    (7, 12, 0, 0):   'PhSm',
    (23, 12, 0, 0):  'Data Parameters PhSm',
    (32, 0, 0, 0):   'Instrument Parameters',
    (40, 0, 0, 0):   'Instrument Parameters Rf',
    (48, 0, 0, 0):   'Acquisition Parameters',
    (64, 0, 0, 0):   'FT - Parameters',
    (96, 0, 0, 0):   'Optics Parameters',
    (104, 0, 0, 0):  'Optics Parameters Rf',
    (160, 0, 0, 0):  'Sample Parameters',
}

def unpack(s,d):
    a = struct.unpack(s,d)
    assert(len(a) == 1)
    return a[0]

def get_ibloci(bindata, i):
    return struct.unpack('<12B', bindata[12*i:12*(i+1)])

def get_ibloc(bindata, i):
    return struct.unpack('<4B2i', bindata[12*i:12*(i+1)])

def get_float_array(b):
    assert(len(b)%4 == 0)
    # print(len(b)//4)
    return(np.array(struct.unpack('<'+str(len(b)//4) + 'f', b),'f'))

def get_params(b):
    inext = 0
    imax = len(b)
    dtypes = {}
    dvalues = {}
    key_order = []
    while(inext < imax):
        bname = unpack('4s', b[inext:inext+4]).rstrip('\x00')
        btype = unpack('<H', b[inext+4:inext+6])
        blen = unpack('<H', b[inext+6:inext+8])
        ilast = inext+8+2*blen
        bdata = b[inext+8:ilast]
        if btype == 0:
            if bname == 'END':
                break
            bknown =  unpack('<i', bdata[0:4])
            bunknown = bdata[4:]
        elif btype == 1:
            bknown =  unpack('<d', bdata[0:8])
            bunknown = bdata[8:]
        elif btype == 2:
            i0 = bdata.find('\x00')
            bknown = bdata[0:i0]
            bunknown = bdata[i0:]
        elif btype == 3:
            i0 = bdata.find('\x00')
            bd = bdata[0:i0].rstrip('\x00')
            if OPUS_TYPES_SPECIFIC.has_key(bname + ',' + bd):
                bknown = OPUS_TYPES_SPECIFIC[bname + ',' + bd]
            elif OPUS_TYPES.has_key(bd):
                bknown = OPUS_TYPES[bd]
            else:
                bknown = '$$$ UNKNOWN ' + bd
            bunknown = bdata[i0:]
        elif btype == 4:
            i0 = bdata.find('\x00')
            bknown = bdata[0:i0]
            bunknown = bdata[i0:]
        else:
            bknown = ''
            bunknown = bdata
        if OPUS_NAMES.has_key(bname):
            humankey =  OPUS_NAMES[bname]
        else:
            humankey = "$$$ Unknown " + bname
        # print(bname, btype, blen, bknown, bunknown)
        # print(btype, blen, humankey, bknown)
        assert (not dtypes.has_key(bname))
        dtypes[bname] = btype
        dvalues[bname] = bknown
        key_order.append(bname)
        # print(humankey + ' | ' + str(bknown))
        inext = ilast
    reste = b[inext:]
    return (dvalues,dtypes,key_order,reste)


def opus_read(filename):
    """
    Return blocs from an opus file.

    calling:
        blocs, dict_params, dict_data, dict_unclassified = opus_read(filename)

    parameters:
        filename - an OPUS binary filename

    returns:
        blocs - list of opus blocs
        dict_params - dictionary of parameter bloc indexes
        dict_data - dictionary of data bloc indexes
        dict_unclassified - dictionary of unclassified (yet) bloc indexes

    blocs[0]: (h0, h1, bloc_indexes list)
    parameter bloc: (key-value dictionary, key-type dictionary, key_order, extra bloc bytes)
    data bloc: numpy float-32 array
        
    h0: 12 bytes tuple
    h1: 12 bytes tuples
    bloc_index: (t0, t1, t2, t3, length, index)

    (t0, t1, t2, t3): 4 bytes tuple indexing dict_params and dict_data
    length: bloc length
    index: bloc position in Opus file
    """
    
    with open(filename, 'rb') as file:
        bindata = file.read()
    h_unknown0 = get_ibloci(bindata, 0)
    h_unknown1 = get_ibloci(bindata, 1)
    h_pointers = []
    # print(h_unknown1)
    ib = 1
    (u0, u1, u2, u3, l, p) = get_ibloc(bindata, ib+1)
    h_pointers.append((u0, u1, u2, u3, l, p))
    pn = p+4*l
    nb = len(bindata)
    rblocs = []
    btypes = []
    while(pn < nb):
        rblocs.append(bindata[p:pn])
        btypes.append((u0, u1, u2, u3))
        ib += 1
        (u0, u1, u2, u3, l, p) = get_ibloc(bindata, ib+1)
        h_pointers.append((u0, u1, u2, u3, l, p))
        assert(p == pn)
        pn = p+4*l
    assert (pn == nb)
    rblocs.append(bindata[p:pn])
    btypes.append((u0, u1, u2, u3))

    lasti = len(rblocs)*3*4
    assert (rblocs[0][lasti:] == '\x00' * (len(rblocs[0]) - lasti))

    blocs = []
    params_dict = {}
    data_dict = {}
    unclassified_dict = {}

    blocs.append((h_unknown0, h_unknown1, h_pointers))
    
    for ib in range(1,len(rblocs)):
        b = rblocs[ib]
        # print('')
        # print('bloc ' + str(ib) + ':')
        # print(btypes[ib], ib, len(rblocs[ib]))
        if b[-8:-5] == 'END':
            # print('params')
            dvalues,dtypes,key_order,reste = get_params(rblocs[ib])
            blocs.append((dvalues,dtypes,key_order,reste))
            if btypes[ib] != (0,0,0,0):            
                assert(not params_dict.has_key(btypes[ib]))
                params_dict[btypes[ib]] = ib
            elif params_dict.has_key(btypes[ib]):
                params_dict[btypes[ib]].append(ib)
            else:
                params_dict[btypes[ib]] = [ib] 
        elif btypes[ib][0] == 7:
            # print('data')
            blocs.append(get_float_array(rblocs[ib]))
            assert(not data_dict.has_key(btypes[ib]))
            data_dict[btypes[ib]] = ib
        else:
            # print('unclassified')
            blocs.append(rblocs[ib])
            unclassified_dict[btypes[ib]] = ib
    return(blocs, params_dict, data_dict, unclassified_dict)

OPUS_NAMES_WIDTH = max(map(len, OPUS_NAMES.values()))

def opus_print(blocs, params_dict, data_dict, unclassified_dict):
    for ik in [
            (23, 4, 0, 0),
            (7, 4, 0, 0),
            (7, 12, 0, 0),
            (23, 12, 0, 0),
            (23, 8, 0, 0),
            (7, 8, 0, 0),
            (48, 0, 0, 0),
            (96, 0, 0, 0),
            (40, 0, 0, 0),
            (104, 0, 0, 0),
            (64, 0, 0, 0),
            (160, 0, 0, 0),
            (32, 0, 0, 0),
            (0, 0, 104, 64),]:
        print('')
        print('*' * len(OPUS_BLOCS[ik]))
        print(OPUS_BLOCS[ik])
        print('*' * len(OPUS_BLOCS[ik]))
        if params_dict.has_key(ik):
            b = blocs[params_dict[ik]]
            for k in b[2]:
                print(format(OPUS_NAMES[k], str(OPUS_NAMES_WIDTH) + 's') + ' | ' +  str(b[0][k]))
        elif data_dict.has_key(ik):
            ib = data_dict[ik]
            print('len(blocs[' + str(ib) + ']) = ' + str(len(blocs[ib])))
        else:
            b = blocs[unclassified_dict[ik]].split('\x00')
            for line in b:
                if len(line) != 0:
                    print(line)
    ik = (0,0,0,0)
    for ib in params_dict[ik]:
        print('')
        title = OPUS_BLOCS[ik] + ': ' + str(ib)
        print('*' * len(title))
        print(title)
        print('*' * len(title))
        b = blocs[ib]
        for k in b[2]:
            print(format(OPUS_NAMES[k], str(OPUS_NAMES_WIDTH) + 's') + ' | ' +  str(b[0][k]))

filename = '/home/fab/Z/Fabrice/2015-10-21_murs/15h14_GaAs_40_deg_38_deg_optim_detect_0.5mm_v7_air.0'
blocs, params_dict, data_dict, unclassified_dict = opus_read(filename)
opus_print(blocs, params_dict, data_dict, unclassified_dict)
