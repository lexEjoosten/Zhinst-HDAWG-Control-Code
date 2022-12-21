
"""
Sets which BLACS_tab belongs to each labscript device.
"""

import labscript_devices

labscript_devices.register_classes(
    'HDAWG',
    BLACS_tab='labscript_devices.HDAWG.blacs_tab.HDAWGTab',
    runviewer_parser=None
)