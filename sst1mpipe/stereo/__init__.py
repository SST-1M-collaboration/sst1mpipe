"""Stereo events of the SST-1M telescopes: combination of the coincident events of the telescopes"""
from .combiner import (
    SlidingWindowCombiner,
    StereoEventCombiner,
    SWATEventIdCombiner,
    WhiteRabbitClosestCombiner,
    nearest_index,
)

__all__ = [
    "SlidingWindowCombiner",
    "StereoEventCombiner",
    "SWATEventIdCombiner",
    "WhiteRabbitClosestCombiner",
    "nearest_index",
]
