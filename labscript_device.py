
import math
import numpy as np
from labscript import Device, LabscriptError, set_passed_properties
from labscript import config    
class HDAWG(Device):
    """A labscript_device for the HDAWG such that it can be controlled via python.
    """
    description = 'HDAWG via python'
    def __init__(self, name):
        Device.__init__(self, name, None, None)
        self.name = name
        self.BLACS_connection = 'HDAWG' # This is needed to create the blacs_tab and will appear as [conn: self.BLACS_connection]
        self.IQLoad = []
        self.dim = []
        
    def generate_code(self, hdf5_file):
        # Code wich is executed by runmanager to create/srote date in h5 file.
        if self.IQLoad != []:
            # create group with the name of device and its data
            grp = hdf5_file.create_group('/devices/'+self.name)
            grp.create_dataset('IQLoad', data = self.IQLoad)

    # takes list of pulses and loads IQSeq to h5 file
    def ks_pulse(self, pulse_dict):
        pulse_program = np.zeros(7)
        pulse_program[0] = pulse_dict['ks_frequency1']
        pulse_program[1] = pulse_dict['ks_frequency2']
        pulse_program[2] = pulse_dict['ks_amplitude1']
        pulse_program[3] = pulse_dict['ks_amplitude2']
        pulse_program[4] = pulse_dict['ks_phase1']
        pulse_program[5] = pulse_dict['ks_phase2']
        pulse_program[6] = pulse_dict['dt']
        return np.array(pulse_program)
    
    # takes list of pulses and loads IQSeq to h5 file
    def pulse_seq(self, pulse_list):
        # start with empyt list
        IQSeq = []
        # append the init pulse (first in sequence)
        IQSeq.append(self.ks_pulse(pulse_list[0]))
        # if we have more to do, check what to do
        if len(pulse_list) > 1:
            for i in range(0, len(pulse_list)-1):
                # for a similar pulse, add the time to last entry of IQSeq
                if pulse_list[i]['ks_frequency1'] == pulse_list[i+1]['ks_frequency1'] and\
                    pulse_list[i]['ks_frequency2'] == pulse_list[i+1]['ks_frequency2'] and\
                    pulse_list[i]['ks_amplitude1'] == pulse_list[i+1]['ks_amplitude1'] and\
                    pulse_list[i]['ks_amplitude2'] == pulse_list[i+1]['ks_amplitude2'] and\
                    pulse_list[i]['ks_phase1'] == pulse_list[i+1]['ks_phase1'] and\
                    pulse_list[i]['ks_phase2'] == pulse_list[i+1]['ks_phase2']:
                    
                    IQSeq[-1][6] += pulse_list[i+1]['dt']
                # for a new pulse: append to IQSeq
                else:
                    IQSeq.append(self.ks_pulse(pulse_list[i+1]))
        
        # finally, load IQSeq to h5-file
        return self.IQLoadAndPlay(IQSeq)

    def IQLoadAndPlay(self, IQSeq):
        self.IQLoad.append(IQSeq)
        # self.dim.append(len(pulse_program))
            
  
        