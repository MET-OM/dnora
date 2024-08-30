from abc import ABC, abstractmethod


def DataProcessor(ABC):
    @abstractmethod
    def __call__(self, obj):
        """aaa"""

    @abstractmethod
    def __str__(self):
        """Describes how the spectral values as processed"""
        pass
