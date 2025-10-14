from dnora import msg

from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import add_time, add_frequency, add_datavar

from dnora.process.spectra import SpectralProcessor
import geo_parameters as gp
from typing import Union, Optional


@add_datavar(name=gp.wave.Ef("spec"), coord_group="all", default_value=0.0)
@add_datavar(name=gp.wave.Dirm("dirm"), coord_group="all", default_value=0.0)
@add_datavar(name=gp.wave.Spr("spr"), coord_group="all", default_value=0.0)
@add_frequency(grid_coord=False)
@add_time(grid_coord=True)
class Spectra1D(PointSkeleton):
    def process(self, spectral_processors: Optional[list[SpectralProcessor]] = None):
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
            new_spec, new_dirs, new_freq, new_inds, new_time, new_spr = processor(
                self.spec(squeeze=False),
                self.dirm(squeeze=False),
                self.freq(),
                self.inds(),
                self.time(),
                self.spr(squeeze=False),
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
                        f"ind = {ind}, {self.core.x_str} = {lon[ind]}, {self.core.y_str} = {lat[ind]}"
                    )

            metadata = self.meta.get()

            self._init_structure(
                x=self.x(strict=True, inds=new_inds),
                y=self.y(strict=True, inds=new_inds),
                lon=self.lon(strict=True, inds=new_inds),
                lat=self.lat(strict=True, inds=new_inds),
                time=new_time,
                freq=new_freq,
            )
            self.set_spec(new_spec)
            self.set_dirm(new_dirs)
            self.set_spr(new_spr)
            self.meta.set(metadata)  # Global attributes

        return
