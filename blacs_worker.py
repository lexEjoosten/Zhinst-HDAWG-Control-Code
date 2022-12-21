from blacs.tab_base_classes import Worker
import labscript_devices.HDAWG.HDAWG_funcs as fn

class HDAWGworker(Worker):
    def init(self):
        global properties; import labscript_utils.properties
        global h5_lock, h5py; import labscript_utils.h5_lock, h5py
        global np; import numpy as np
        global math; import math
        global threading; import threading
        global textwrap; import textwrap
        global zhcore; import zhinst.core
        global zhutils; import zhinst.utils

        self.device_id = "dev8542"  # Device serial number available on its rear panel.
        self.interface = "1GbE"  # For Ethernet connection.
        #self.interface = "USB" # For all instruments connected to the host computer via USB.

        server_host = "localhost"
        server_port = 8004
        api_level = 6  # Maximum API level supported for all instruments.



        (self.daq, self.device, _) = zhinst.utils.create_api_session(
        self.from_childdevice_id, api_level, server_host=server_host, server_port=server_port
        )
        # Create an API session to the Data Server.
        self.daq = zhinst.core.ziDAQServer(server_host, server_port, api_level)
        # Establish a connection between Data Server and Device.
        self.daq.connectDevice(self.device_id, self.interface)


        zhinst.utils.api_server_version_check(self.daq)

        # Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
        zhinst.utils.disable_everything(self.daq, self.device)


        # Create an instance of the AWG Module
        awgModule = self.daq.awgModule()        
        awgModule.set("device", self.device)    
        awgModule.execute()


        exp_setting = [
            [f"/{self.device_id}/sigouts/*/on", 1], # turn on the outputs   x
            [f"/{self.device_id}/sigouts/*/range", 1], #set the maximum voltage to 1Vpp
            [f"/{self.device_id}/awgs/0/outputs/*/amplitude", 1], #set the amplitude to 1, (this means that a value of 1 in sequencer corresponds to 0.5V on the output, as this is 1 Vpp for a sine function)
            [f"/{self.device_id}/awgs/0/time", 0], 
            [f"/{self.device_id}/awgs/*/enable", 0],
            [f"/{self.device_id}/awgs/0/userregs/0", 0],
            [f"/{self.device_id}/awgs/1/time", 0],
            [f"/{self.device_id}/awgs/1/userregs/0", 0],
            [f"/{self.device_id}/triggers/out/0/source", 4], #Marker output sent to the marker port at the front of the device
            [f"/{self.device_id}/triggers/out/1/source", 6], #Marker output sent to the marker port at the front of the device
            [f"/{self.device_id}/triggers/out/2/source", 4], #Marker output sent to the marker port at the front of the device
            [f"/{self.device_id}/triggers/out/3/source", 6], #Marker output sent to the marker port at the front of the device
            [f"/{self.device_id}/awgs/0/outputs/*/modulation/mode", 0], # Channels 1 and 2 are not modulated
            [f"/{self.device_id}/awgs/1/outputs/0/modulation/mode", 1], # channels 3 and 4 are modulated by sine generator 2
            [f"/{self.device_id}/awgs/1/outputs/1/modulation/mode", 2], # channels 3 and 4 are modulated by sine generator 2
            [f"/{self.device_id}/system/awg/oscillatorcontrol",1], # oscillators are controlled by the API not, GUI
            [f"/{self.device_id}/system/awg/channelgrouping", 2] # Groups channels in 1x4 configuration
        ]
        self.daq.set(exp_setting)
        self.daq.sync()
        

    def program_AWG(self):
        with h5py.File(self.h5_file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+self.device_name]
            # read out the dimension of the matrix
            # dim = group['dim'][()]
    
            # generate the structure of pulse_program
            # pulse_program = np.zeros((int(dim),7))
            # load the data from h5 file
            ### read_direct: old version
            # group['IQLoad'].read_direct(pulse_program)
            # print(pulse_program) ### testing
            
            ### testing: just load the data
            pulse_program = group['IQLoad'][()]
            print(pulse_program)
            # send functions to iq_modulator
            
            ### we need better input control
            ### how is it saved to h5 file?
            self.PlaySequence(pulse_program.reshape(np.shape(pulse_program)[1], 7))

    def program_manual(self):
        return None
    
    def transition_to_buffered(self, device_name, h5file, initial_values, fresh):
        self.h5_file = h5file
        self.device_name = device_name
        with h5py.File(h5file, 'r') as hdf5_file:
            group = hdf5_file['/devices/'+self.device_name]
            # if there is no data assigned to iq_modulator: do nothing
            if 'IQLoad' in list(group.keys()):        
                # start the programming 0.1s into the experiment
                thread = threading.Timer(0.1, self.program_AWG)
                thread.start()
            else:
                self.IQModulatorOutputOff()
                print('no IQLoad found for this run')

        return initial_values

    def abort_buffered(self):
        seq=fn.Sequence(While=False)
        seq.Upload()
        return self.transition_to_manual(True)
            
    def transition_to_manual(self, abort = False):    
        return True
        
    def shutdown(self):
        print('shutdown')
        exp_setting = [
            [f"/{self.device_id}/system/shutdown", 1],]
        self.daq.set(exp_setting)
        self.daq.sync()
        return True
    
    def Playsequence(self,Pulses):

        SEQ=fn.Sequence()
        Mark=[0,0,0,0]

        A1=0
        A2=0

        F1=0
        F2=0
        
        P1=0
        P2=0

        for P in Pulses:
            
            
            pulse_dict['ks_frequency1'] = P[0]
            pulse_dict = P[1]
            pulse_dict['ks_amplitude1'] = P[2]
            pulse_dict['ks_amplitude2'] = P[3] 
            pulse_dict['ks_phase1'] = P[4]
            pulse_dict['ks_phase2'] = P[5]
            pulse_dict['dt'] = P[6]
            
            t=0
            dt=pulse_dict["dt"]
            if dt==None or dt==0:
                raise(Exception("Pulses need to have non-zero duration"))
            #AUX=[pulse_dict["aux1"],pulse_dict["aux2"],pulse_dict["aux3"],pulse_dict["aux4"]]
            #for m in range(AUX):
            #    if AUX[m]!=None:
            #        Mark[m]=AUX[m]
            

            if pulse_dict["ks_amplitude1"]!=None:
                A1=pulse_dict["ks_amplitude1"]
            if pulse_dict["ks_amplitude2"]!=None:
                A2=pulse_dict["ks_amplitude2"]

            if pulse_dict["ks_frequency1"]!=None:
                F1=pulse_dict["ks_frequency1"]
            if pulse_dict["ks_frequency2"]!=None:
                F2=pulse_dict["ks_frequency2"]

            if pulse_dict["ks_phase1"]!=None:
                P1=pulse_dict["ks_phase1"]
            if pulse_dict["ks_phase2"]!=None:
                P2=pulse_dict["ks_phase2"]

            if A1!=0 and A2!=0:
                phase11=P1
                phase21=P2
                if dt>1e-4:
                    print("long duration two-tones can cause unpredictable memory issues for the HDAWG")

                SEQ.addtwotone(dur=dt,phase11=phase11,phase21=phase21,amp1=A1,amp2=A2,freq1=F1,freq2=F2,phasetime=t,Mark=Mark)
            
            elif A1!=0:
                phase1=P1
                SEQ.singletone(dur=dt,amp=A1,freq=F1,phase1=phase1,phasetime=t,Mark=Mark)

            elif A2!=0:
                phase1=P2
                SEQ.singletone(dur=dt,amp=A1,freq=F2,phase1=phase1,phasetime=t,Mark=Mark)
            
            else:
                SEQ.addwait(dur=dt,Mark=Mark)
            t+=dt
        SEQ.Upload(self.device,self.daq)