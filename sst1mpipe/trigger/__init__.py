"""Software emulation of the SST-1M camera trigger (patch7 and TDSCAN), see emulator.py."""
from sst1mpipe.trigger.emulator import TriggerEmulator, TriggerPatch7Emulator, TriggerTDSCANEmulator

__all__ = ["TriggerEmulator", "TriggerPatch7Emulator", "TriggerTDSCANEmulator"]
