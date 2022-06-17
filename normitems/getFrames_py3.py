import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


def figure_handler(fig=None,
                   figout=None,
                   figext=None,
                   dpi=None,
                   pad_inches=0.1,
                   transparent=True):
    """ Handle the saving, drawing, and formating of figure objects.

    Parameters
    ----------
    fig: Figure object
        target figure.
    figout: str
        Name of the figure.
    figext: str
        format of the output figure, choose among 'pdf', 'eps', 'eye', 'jpg',
        and 'png', or it could be any combination of those four, like 'pdfeps'
        or 'epseys', etc.
    dpi: int
        dpi of output image file, especially useful for 'png' files (not
        necessary for vector images like 'pdf' and 'eps').
    pad_inches: float
        padding.
    transparent: bool
        Set true if the transparent axes background is needed.

    Returns
    -------
    printed: bool
        Printed or not.
    """
    printed = False
    if figext is None:
        plt.show()
        printed = True
    elif isinstance(fig, plt.Figure):
        if figout is None:
            for i in range(100):
                # figout = "fig"+str(i)
                figout = "fig" + str(i).zfill(2)
                if (len(glob(figout + ".*")) == 0):
                    break
                if (i == 99):
                    print("too many fig* files in current dir!")
                    printed = False
                    return (printed)
        if "pdf" in figext:
            fig.savefig(figout + ".pdf", format="pdf", dpi=dpi)
            printed = True
        if "jpg" in figext:
            fig.savefig(figout + ".jpg", format="jpg", dpi=dpi)
            printed = True
        if "png" in figext:
            fig.savefig(figout + ".png",
                        transparent=transparent,
                        format="png",
                        dpi=dpi)
            printed = True
        if "eps" in figext:
            fig.savefig(figout + ".eps",
                        papertype="a4",
                        format="eps",
                        bbox_inches='tight',
                        pad_inches=pad_inches,
                        dpi=dpi)
            printed = True
        if "eye" in figext:
            plt.show()
            printed = True
        if (printed is False):
            print("savefig failed, currently only pdf, png, and eps allowed")

    return (printed)


def getFrames(nframe,
              aspectratio="1x1",
              ncol=None,
              wid=10.0,
              x0=0.09,
              y0=0.09,
              x1=0.95,
              y1=0.95,
              xgut=0.1,
              ygut=0.1,
              vsplit=None,
              vsgut=0.0,
              fig=None,
              **axeskwargs):
    """ return `nframe` axes with given geometry and dimension.

    Returns
    -------
    fig: Figure
        Figure
    axes: list
        List of Axes
    ncol: integer
        number of columns
    nrow: integer
        number of rows
    """
    w2h = float(aspectratio.split("x")[0]) / float(aspectratio.split("x")[1])
    if ncol is None:
        if w2h > 1.0 or nframe == 1:
            # wide   -> vertical arrangement
            ncol = 1
        else:
            # narrow -> horizontal arrangement
            ncol = nframe

    nextra = np.mod(nframe, ncol)
    if nextra == 0:
        nrow = np.int( (nframe - nextra) / ncol )
    else:
        nrow = np.int( (nframe - nextra) / ncol + 1 )

    if isinstance(xgut, list):
        # xgut.insert(0, 0)
        xgut.append(0)
        xgut = np.atleast_1d(xgut)

    else:
        xgut = np.atleast_1d( xgut * ncol )

    if isinstance(ygut, list):
        ygut.append( 0 )
        ygut = np.atleast_1d( ygut )

    else:
        ygut = np.atleast_1d( ygut * nrow )

    # width of cell in ax fraction
    w = (x1 - x0 - np.sum(xgut)) / ncol

    # height of cell in ax fraction
    h = (y1 - y0 - np.sum(ygut)) / nrow

    # wid * w  divided by hig * h should be equal to w2h
    hig = wid * w / (w2h * h)

    if fig is None:
        fig = plt.figure(figsize=(wid, hig))

    # left to right, top to bottom
    axes = []
    if vsplit is not None:
        t2b = float(vsplit.split("x")[0]) / float(vsplit.split("x")[1])

    iframe = 0
    for i in range(nrow):
        for j in range(ncol):
            if iframe >= nframe:
                # make sure there is no redundant axes
                break

            rect = [
                x0 + j * w + np.sum(xgut[0:j]),
                y1 - (i + 1) * h - np.sum(ygut[0:i]), w, h
            ]
            if vsplit is None:
                ax = fig.add_axes(rect, **axeskwargs)
                axes.append(ax)

            else:
                _x0, _y0, _w, _h = rect
                htop = (_h - vsgut) * (t2b / (t2b + 1.0))
                hbot = _h - htop - vsgut
                rect_bot = [_x0, _y0, _w, hbot]
                rect_top = [_x0, _y0 + hbot + vsgut, _w, htop]
                ax_top = fig.add_axes(rect_top, **axeskwargs)
                ax_bot = fig.add_axes(rect_bot, **axeskwargs)
                axes.append([ax_top, ax_bot])
            iframe += 1

    return (fig, axes, ncol, nrow)


def getHandyFrames(geometry='1x1-1x1', wid=8, **kwargs):
    """handy frames."""
    geo = geometry.split("-")
    print (geo)

    if len(geo) == 2:
        colrow, aspectratio = geo
        vsplit = None

    elif len(geo) == 3:
        colrow, aspectratio, vsplit = geo

    if colrow == '1x1' and aspectratio == '1x1':
        _kwargs = {
            'x0': 0.12,
            'x1': 1.08,
            'y0': 0.01,
            'y1': 0.98,
        }
    else:
        _kwargs = {}

    mykwargs = _kwargs.copy()
    mykwargs.update(kwargs)

    ncol, nrow = [int(r) for r in colrow.split('x')]
    nframe = ncol * nrow

    return (getFrames(nframe,
                      aspectratio=aspectratio,
                      ncol=ncol,
                      wid=wid,
                      vsplit=vsplit,
                      **mykwargs))


def get_cax(fig, ax, position="right", width=0.02, pad=0.05):
    """ Make axis for colorbars by expanding from existing axis.

    should be deprecated soon and switch to inset_axes in the future.
    https://matplotlib.org/3.1.1/gallery/axes_grid1/demo_colorbar_with_inset_locator.html
    """
    bbox = ax.get_position()
    l, b, w, h = bbox.x0, bbox.y0, bbox.width, bbox.height

    if position == 'right':
        pos = [l + w + pad, b, width, h]

    elif position == 'left':
        pos = [l - pad - width, b, width, h]

    elif position == 'top':
        pos = [l, b+h+pad, w, width]

    elif position == 'bottom':
        pos = [l, b-pad-width, w, width]

    cax = fig.add_axes(pos)

    return(cax)


if __name__ == "__main__":

    fig, axes, ncol, nrow = getHandyFrames('2x2-1x1', wid=10)

    # cax0 = get_cax(fig, axes[1], pad = 0.01)
    # cax1 = get_cax(fig, axes[3], pad = 0.01)

    figure_handler(fig, figout=None, figext=None)

