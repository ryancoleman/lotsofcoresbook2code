from ase.io.eps import EPS


class PNG(EPS):
    def write_header(self):
        from matplotlib.backends.backend_agg import RendererAgg

        try:
            from matplotlib.transforms import Value
        except ImportError:
            dpi = 72
        else:
            dpi = Value(72)

        self.renderer = RendererAgg(self.w, self.h, dpi)

        #self.gc = GraphicsContextBase()
        #self.gc.set_linewidth(2)

    def write_trailer(self):
        renderer = self.renderer
        if hasattr(renderer._renderer, 'write_png'):
            # Old version of matplotlib:
            renderer._renderer.write_png(self.filename)
        else:
            from matplotlib import _png
            # buffer_rgba does not accept arguments from version 1.2.0
            # https://github.com/matplotlib/matplotlib/commit/f4fee350f9fbc639853bee76472d8089a10b40bd
            import matplotlib
            if matplotlib.__version__ < '1.2.0':
                x = renderer._renderer.buffer_rgba(0, 0)
                _png.write_png(renderer._renderer.buffer_rgba(0, 0),
                               renderer.width, renderer.height,
                               self.filename, 72)
            else:
                x = renderer._renderer.buffer_rgba()
                _png.write_png(renderer._renderer.buffer_rgba(),
                               renderer.width, renderer.height,
                               self.filename, 72)

                
def write_png(filename, atoms, **parameters):
    if isinstance(atoms, list):
        if len(atoms) > 1:
            raise RuntimeError("Don't know how to save more than " +
                               "one image to PNG image!")
        else:
            atoms = atoms[0]
    PNG(atoms, **parameters).write(filename)
