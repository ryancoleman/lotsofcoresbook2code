from pylab import *
def twiny(ay=None):
    """
    Make a second axes overlay ay (or the current axes if ay is None)
    sharing the yaxis.  The ticks for ay2 will be placed on the top,
    and the ay2 instance is returned.  See examples/two_scales.py
    """
    if ay is None:
        ay=gca()


    ay2 = gcf().add_axes(ay.get_position(), sharey=ay, frameon=False)
    ay2.xaxis.tick_top()
    ay2.xaxis.set_label_position('top')
    ay.xaxis.tick_bottom()
    draw_if_interactive()
    return ay2
