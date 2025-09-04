import warnings
from functools import wraps

# DEPRECIATE = ["metno", "era5", "cmems", "noaa", "nchmf"]
DEPRECIATE = []


def deprecated_class_call(data_source: str, package: str, data_type: str):
    """
    Decorator to wrap the __call__ method of a class and issue a deprecation warning.
    """

    def decorator(cls):
        classname = cls.__name__
        # Check if the class has a __call__ method
        if hasattr(cls, "__call__"):
            # Wrap the existing __call__ method
            original_call = cls.__call__

            @wraps(original_call)
            def wrapped_call(self, *args, **kwargs):
                message = (
                    f"{data_source} products have been moved to package dnora-{package} ('pip install dnora-{package}'). The 'dnora.read.{data_type}.{package}.{classname}' class is deprecated and will be removed in a future release. "
                    f"Please use 'dnora_{package}.{data_type}.{classname}' instead."
                )
                if package in DEPRECIATE:
                    warnings.warn(
                        message,
                        DeprecationWarning,
                        stacklevel=2,  # Points to the user's code
                    )
                return original_call(self, *args, **kwargs)

            # Replace the original __call__ with the wrapped version
            cls.__call__ = wrapped_call

        return cls

    return decorator
