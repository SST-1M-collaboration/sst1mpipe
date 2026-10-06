"""Software emulation of the SST-1M camera trigger (patch7 and TDSCAN)."""
from sst1mpipe.trigger.emulator import TriggerEmulator, as_sst1m_event

__all__ = ["TriggerEmulator", "as_sst1m_event"]
