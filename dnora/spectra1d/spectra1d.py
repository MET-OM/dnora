from copy import copy

# Import objects
from .process import spectral_processor_for_convention_change
from dnora.spectral_conventions import SpectralConvention, convert_2d_to_1d

# Import abstract classes and needed instances of them
from dnora import msg

from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import add_time, add_frequency, add_datavar

from .process import SpectralProcessor
import geo_parameters as gp


@add_datavar(name="spec", coord_group="all", default_value=0.0)
@add_datavar(name="dirm", coord_group="all", default_value=0.0)
@add_datavar(name="spr", coord_group="all", default_value=0.0)
@add_frequency(grid_coord=False)
@add_time(grid_coord=True)
class Spectra1D(PointSkeleton):
    meta_dict = {"spec": gp.wave.Ef, "dirm": gp.wave.Dirm, "spr": gp.wave.Spr}

    def process_spectra(
        self, spectral_processors: list[SpectralProcessor] | None = None
    ):
        """Process all the individual spectra of the spectra object.

        E.g. change convention form Meteorological to Oceanic, interpolate spectra to, or multiply everything with a constant.
        """

        if spectral_processors is None:
            msg.info("No SpectralProcessor provided. Doing Nothing.")
            return

        if not isinstance(spectral_processors, list):
            spectral_processors = [spectral_processors]

        for processor in spectral_processors:
            msg.process(f"Processing spectra with {type(processor).__name__}")
            print(processor)
            old_convention = processor._convention_in()
            if old_convention is not None:
                if old_convention != self.convention():
                    msg.warning(
                        f"Spectral convention ({self.convention()}) doesn't match that expected by the processor ({old_convention})!"
                    )

            new_spec, new_dirs, new_freq, new_spr, new_inds = processor(
                self.spec(), self.dirm(), self.freq(), self.spr(), self.inds()
            )

            new_inds = list(new_inds)

            if not new_inds:
                self.ds_manager.set_new_ds(None)
                msg.warning(f"No spectra left after processing. Removing all data.")
                return

            del_inds = list(set(self.inds()) - set(new_inds))
            if del_inds:
                lon, lat = self.lonlat(native=True)
                msg.info(f"Removing the following points:")
                for ind in list(set(self.inds()) - set(new_inds)):
                    msg.plain(
                        f"ind = {ind}, {self.x_str} = {lon[ind]}, {self.y_str} = {lat[ind]}"
                    )

            metadata = self.metadata()

            self._init_structure(
                x=self.x(strict=True),
                y=self.y(strict=True),
                lon=self.lon(strict=True),
                lat=self.lat(strict=True),
                time=self.time(),
                freq=new_freq,
            )
            self.set_spec(new_spec)
            self.set_dirm(new_dirs)
            self.set_spr(new_spr)
            self.set_metadata(metadata)  # Global attributes

            # Set new convention if the processor changed it
            new_convention = processor._convention_out()
            if new_convention is not None:
                self._mark_convention(new_convention)
                msg.info(
                    f"Changing convention from {old_convention} >>> {new_convention}"
                )

        return

    def set_convention(self, convention: SpectralConvention) -> None:
        """Processes spectra to new directional spectral convention and updates metadata."""
        spectral_processor = spectral_processor_for_convention_change(
            current_convention=self.convention(),
            wanted_convention=convert_2d_to_1d(convention),
        )

        if spectral_processor is None:
            msg.info(
                f"Convention ({self.convention()}) already equals wanted convention ({convention})."
            )
            return

        self.process_spectra(spectral_processor)

    def _mark_convention(self, convention: SpectralConvention) -> None:
        """Marks new convention in metadata etc. but does nothing to the spectra"""
        self._convention = convention
        self.set_metadata({"spectral_convention": self.convention().value}, append=True)
        print(f"Spectral convention is now: {self.convention()}")

    def convention(self) -> SpectralConvention | None:
        """Returns the convention (OCEAN/MET/MATH) of the spectra"""
        if not hasattr(self, "_convention"):
            return None
        return copy(self._convention)

    # def process_spectra(self, spectral_processors: List[SpectralProcessor] = None):
    #     """Process all the individual spectra of the spectra object.

    #     E.g. change convention form WW3 to Oceanic, interpolate spectra to
    #     new frequency grid, or multiply everything with a constant.
    #     """

    #     if spectral_processors is None:
    #         msg.info("No SpectralProcessor provided. Doing Nothing.")
    #         return

    #     if not isinstance(spectral_processors, list):
    #         spectral_processors = [spectral_processors]

    #     convention_warning = False

    #     for processor in spectral_processors:
    #         msg.process(f"Processing spectra with {type(processor).__name__}")
    #         # self._history.append(copy(processor))
    #         old_convention = processor._convention_in()
    #         if old_convention is not None:
    #             if old_convention != self.convention():
    #                 msg.warning(
    #                     f"Spectral convention ({self.convention()}) doesn't match that expected by the processor ({old_convention})!"
    #                 )
    #                 convention_warning = True

    #         new_spec, new_dirs, new_freq, new_spr = processor(
    #             self.spec(), self.mdir(), self.freq(), self.spr()
    #         )
    #         metadata = self.metadata()

    #         self._init_structure(
    #             x=self.x(strict=True),
    #             y=self.y(strict=True),
    #             lon=self.lon(strict=True),
    #             lat=self.lat(strict=True),
    #             time=self.time(),
    #             freq=new_freq,
    #         )
    #         self.set_spec(new_spec)
    #         self.set_mdir(new_dirs)
    #         self.set_spr(new_spr)
    #         self.set_metadata(metadata)  # Global attributes

    #         # Set new convention if the processor changed it
    #         new_convention = processor._convention_out()
    #         if new_convention is not None:
    #             self._set_convention(new_convention, process=False)
    #             if convention_warning:
    #                 msg.warning(
    #                     f"Convention variable set to {new_convention}, but this might be wrong..."
    #                 )
    #             else:
    #                 msg.info(
    #                     f"Changing convention from {old_convention} >>> {new_convention}"
    #                 )

    #         print(processor)

    #         msg.blank()
    #     return

    # def _set_convention(
    #     self, convention: SpectralConvention, process: bool = True
    # ) -> None:
    #     """Sets a new spectral directional convention. To not touch spectra, use process=False."""
    #     if isinstance(convention, str):
    #         convention = SpectralConvention[convention.upper()]

    #     spectral_processor = spectral_processor_for_convention_change(
    #         current_convention=self.convention(),
    #         wanted_convention=convert_2d_to_1d(convention),
    #     )

    #     if spectral_processor is None:
    #         msg.info(
    #             f"Convention ({self.convention()}) already equals wanted convention ({convention})."
    #         )
    #         return

    #     if process:
    #         self.process_spectra(spectral_processor)
    #     else:
    #         self._convention = convention
    #         self.set_metadata(
    #             {"spectral_convention": self.convention().value}, append=True
    #         )
    #         print(f"Spectral convention is now: {self.convention()}")

    # def convention(self):
    #     """Returns the convention (OCEAN/MET/MATH) of the spectra"""
    #     if not hasattr(self, "_convention"):
    #         return None
    #     return copy(self._convention)

    # def __str__(self) -> str:
    #     """Prints status of spectra."""

    #     msg.header(self, f"Status of spectra {self.name}")
    #     if self.x() is not None:
    #         msg.plain(f"Contains data ({len(self.x())} points) for {self.start_time()} - {self.end_time()}")
    #         msg.plain(f"Data covers: lon: {min(self.lon())} - {max(self.lon())}, lat: {min(self.lat())} - {max(self.lat())}")
    #     if len(self._history) > 0:
    #         msg.blank()
    #         msg.plain("Object has the following history:")
    #         for obj in self._history:
    #             msg.process(f"{obj.__class__.__bases__[0].__name__}: {type(obj).__name__}")
    #     #msg.print_line()
    #     #msg.plain("The Boundary is for the following Grid:")
    #     #print(self.grid())

    #     msg.print_line()

    #     return ''
