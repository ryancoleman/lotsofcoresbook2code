# -*- coding: utf-8 -*-

from ase.neb import get_NEB_plot


def NudgedElasticBand(images):
    fig = get_NEB_plot(images)
    fig.show()
