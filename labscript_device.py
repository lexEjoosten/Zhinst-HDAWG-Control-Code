
import math
import numpy as np
from labscript import Device, LabscriptError, set_passed_properties
from labscript import config    
class VFG150(Device):
    """A labscript_device for the VFG150 such that it can be controlled via python.
    """
    description = 'VFG150 via python'
    
    def __init__(self, name):
        Device.__init__(self, name, None, None)
        self.name = name
        self.BLACS_connection = 'VFG150' # This is needed to create the blacs_tab and will appear as [conn: self.BLACS_connection]
        self.S = []
        
    def generate_code(self, hdf5_file):
        # Code wich is executed by runmanager to create/store data in h5 file.
        # create group with the name of device
        grp = hdf5_file.create_group('/devices/'+self.name)
        # Write data to new dataset
        if self.S != []:
            grp.create_dataset('S', data = self.S)

    def freq(self,frequency):
        self.S +=f_freq(frequency)
    
    def ampl(self,amplitude):
        self.S +=f_ampl(amplitude)
        
    def ampl_slope(self,slope):
        self.S +=f_ampl_slope(slope)
        
    def aux(self,port, state):
        self.S +=f_aux(port, state)
    
    def interval(self,time):
        self.S +=f_interval(time)
        
    def phase(self,phi):
        self.S +=f_phase(phi)
    
    def reset_timebase(self):
        self.S +=f_reset_timebase()
        
    def set_phase_continuous(self):
        self.S +=f_set_phase_continuous()
        
    def set_phase_coherent(self):
        self.S +=f_set_phase_coherent()
    
    def trigger(self):
        self.S +=f_trigger()
    
    def exponentialFrequencyRamp(self,initialFrequency, finalFrequency, dropTime, pulseDuration, maxPhaseError = np.pi/2, amplitude = 1):
        self.S +=f_exponentialFrequencyRamp(initialFrequency, finalFrequency, dropTime, pulseDuration, maxPhaseError, amplitude)
       
    def vfg_pulse(self,old_pulse_dict, new_pulse_dict):
        self.S +=f_vfg_pulse(old_pulse_dict, new_pulse_dict)
    
    def first_vfg_pulse(self,pulse_dict):
        self.S +=f_first_vfg_pulse(pulse_dict)
    
    def pulse_seq(self,pulse_list):
        self.S +=f_pulse_seq(pulse_list)
            
            
  
        