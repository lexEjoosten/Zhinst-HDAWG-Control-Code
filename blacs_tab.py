from blacs.device_base_class import DeviceTab

class HDAWGTab(DeviceTab):
    def initialise_GUI(self):
        pass
    
    def initialise_workers(self):
        self.create_worker(
            'main_worker',
            'labscript_devices.HDAWG.blacs_worker.HDAWGWorker'
        )
        self.primary_worker = 'main_worker'    
    
    
