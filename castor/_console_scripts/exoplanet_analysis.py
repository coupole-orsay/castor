#!/usr/bin/env python3

import argparse
import os

from tqdm import tqdm
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.ioff()

from castor import files_handling, preparation, alignment, photometry

def get_parsed_args():
    parser = argparse.ArgumentParser(
        description='Exoplanet analysis tool.')
    parser.add_argument(
        'name', type=str,
        help='Name of the exoplanet')
    parser.add_argument(
        '--sci-path', type=str,
        help='Directory containing the science FITS. Default: {name}/sci')
    parser.add_argument(
        '--sci-dark-path', type=str,
        help='Directory containing the science dark FITS. Default: {name}/sci_dark')
    parser.add_argument(
        '--flat-path', type=str,
        help='Directory containing the flat FITS. Default: {name}/flat')
    parser.add_argument(
        '--flat-dark-path', type=str,
        help='Directory containing the flat dark FITS. Default: {name}/flat_dark')
    parser.add_argument(
        '--output-path', type=str,
        help='Directory where output files are saved. Default: {name}/')
    parser.add_argument(
        '--sep-threshold', type=float, default=20,
        help=('Source extraction threshold, passed to sep.extract(). '
              'https://sep.readthedocs.io/en/latest/api/sep.extract.html '
              'Default: 20. '))
    parser.add_argument(
        '--running-average', type=str, default='none',
        help=('The type of running average applied to the transit light curve. '
              'Options: none, boxcar, triangle. '
              'Default: none '))
    parser.add_argument(
        '--running-average-width', type=float, default=60,
        help=('Time width of the running average, in seconds. '
              'Default: 60 '))

    args = parser.parse_args()

    if args.sci_path is None:
        args.sci_path = os.path.join(args.name, 'sci')
    if args.flat_path is None:
        args.flat_path = os.path.join(args.name, 'flat')
    if args.sci_dark_path is None:
        args.sci_dark_path = os.path.join(args.name, 'sci_dark')
    if args.flat_dark_path is None:
        args.flat_dark_path = os.path.join(args.name, 'flat_dark')
    if args.output_path is None:
        args.output_path = args.name

    return args

# Plot tools

def source_label(i):
    if i == 0:
        return 'exo'
    else:
        return 'ref {}'.format(i)

def short_source_label(i):
    if i == 0:
        return 'exo'
    else:
        return str(i)

class ImageSeriesPlot():
    def __init__(self, images, default_image=0):
        ''' Represent a series of images of the same field of view.

        Parameters
        ==========
        images : 3D arrray of size MxNxt
            The series of images.
        default_image : float < t (default: 0)
            The image of the series to use by default in plots where a single
            image is needed.
        '''
        self.images = images
        self.default_image = default_image

    def coverage(self):
        coverage = np.sum(np.isfinite(self.images), axis=0)
        coverage = coverage / len(self.images)
        return coverage

    def plot_fov(self, ax, frame='default', coverage_transparency=True,
                 coverage_contours=True, **kwargs):
        ''' Plot the field of view 

        Parameters
        ==========
        ax : matplotlib axes
            Axes in which to plot the field of view.
        frame : None, int, or func (default: 'default')
            The frame to display.
            If None, use self.default_image.
            If int, use the corresponding image of the series.
            If func, it must take as parameter the cube of images, and return a
            single frame. (e.g. functools.partial(np.mean, axis=-1))
        coverage_transparency : bool (default: True)
            If True, the transparency of a given pixel is proportional to its
            coverage (i.e. how many times it is finite in the series of images).
        coverage_contours : bool (default: True)
            If True, show the 90%, 95%, 99%, and 100% coverage contours.
        **kwargs : passed to plt.imshow()
        '''

        coverage = self.coverage()

        if frame == 'default':
            image = self.images[self.default_image]
        elif isinstance(frame, int):
            image = self.images[frame]
        elif callable(frame):
            image = frame(self.images)
        else:
            raise ValueError('Invalid frame type.')

        vmin, vmax = np.nanpercentile(image, [0.25, 99.75])
        default_norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        if coverage_transparency:
            norm = kwargs.pop('norm', default_norm)
            cmap = mpl.cm.get_cmap(kwargs.pop('cmap', None))
            norm(image[np.isfinite(image)])
            image = norm(image)
            image = cmap(image)
            image[:, :, -1] = coverage**2

        im = ax.imshow(image, **kwargs)

        if coverage_contours:
            clevels = np.array([.9, .95, .99, 1]) - 1e-10
            cs = ax.contour(coverage, levels=clevels,
                            colors='white', linewidths=.5)
            ax.clabel(cs, fmt='%1.2f', fontsize='smaller')

class Source():
    def __init__(self, parent, radius, annulus_radii=None,
                 color_active='#6bc341', color_inactive='#b6e2a2'):
        self.parent = parent
        self.is_active = False
        self.radius = radius
        self.annulus_radii = annulus_radii
        self.color_active = color_active
        self.color_inactive = color_inactive

        self.circ = plt.Circle(
            (None, None), self.radius,
            fill=False, color=self.color_inactive)
        self.parent.im_ax.add_patch(self.circ)

        if self.annulus_radii is not None:
            inner_radius, outer_radius = self.annulus_radii
            inner_circ = plt.Circle(
                (None, None), inner_radius,
                fill=False, color=self.color_inactive)
            outer_circ = plt.Circle(
                (None, None), outer_radius,
                fill=False, color=self.color_inactive)
            self.parent.im_ax.add_patch(inner_circ)
            self.parent.im_ax.add_patch(outer_circ)
            self.annulus_circs = (inner_circ, outer_circ)
        else:
            self.annulus_circs = None

    def set_coordinates(self, x, y):
        self.circ.center = x, y
        if self.annulus_circs is not None:
            for circ in self.annulus_circs:
                circ.center = x, y
        self.parent.fig.canvas.draw()

    def get_coordinates(self):
        return self.circ.center

    def circ_contains(self, event):
        return self.circ.contains(event)[0]

    def annulus_contains(self, event):
        try:
            return self.annulus_circs[1].contains(event)[0]
        except TypeError:
            raise ValueError('This source contains no annulus.')

    def contains(self, event):
        if self.annulus_circs is not None:
            return self.annulus_contains(event)
        else:
            return self.circ_contains(event)

    def set_active(self):
        self.is_active = True
        self.circ.set_color(self.color_active)
        if self.annulus_circs is not None:
            for circ in self.annulus_circs:
                circ.set_color(self.color_active)

    def set_inactive(self):
        self.is_active = False
        self.circ.set_color(self.color_inactive)
        if self.annulus_circs is not None:
            for circ in self.annulus_circs:
                circ.set_color(self.color_inactive)

    def remove(self):
        self.circ.remove()
        if self.annulus_circs is not None:
            for circ in self.annulus_circs:
                circ.remove()

    def plot_label(self, label, **kwargs):
        x, y = self.circ.center
        if self.annulus_radii is not None:
            r = self.annulus_radii[-1]
        else:
            r = self.radius
        if 'verticalalignment' not in kwargs:
            kwargs.update(verticalalignment='center')
        return self.parent.im_ax.text(x + r * 1.5, y, label, **kwargs)

class MultipleSourcesSelector():
    def __init__(self, images_plot, circ_radius=15, annulus_radii=None,
            mandatory_sources=1):
        self.images_plot = images_plot
        self.circ_radius = circ_radius
        self.annulus_radii = annulus_radii
        self.mandatory_sources = mandatory_sources

    def _inactivate_all_sources(self):
        for source in self.sources:
            source.set_inactive()
        self.active_source = None

    def _activate_source(self, source):
        self._inactivate_all_sources()
        source.set_active()
        self.active_source = source

    def _add_source_callback(self, event):
        source = Source(self, self.circ_radius,
                        annulus_radii=self.annulus_radii)
        source.set_coordinates(event.xdata, event.ydata)
        self.sources.append(source)
        self._activate_source(source)

    def _remove_source_callback(self, event):
        if len(self.sources) > self.mandatory_sources:
            source = self.active_source
            source.remove()
            self.sources.remove(source)
            new_source = self.sources[-1]
            self._activate_source(new_source)
            self.fig.canvas.draw()

    def _button_press_callback(self, event):
        if event.inaxes is not self.im_ax:
            return

        if event.button == 1:
            if self.active_source is not None:
                self.active_source.set_coordinates(event.xdata, event.ydata)

        elif event.button == 3:
            sources = []
            for source in self.sources:
                if source.contains(event):
                    sources.append(source)
            if len(sources) == 1:
                self._activate_source(sources[0])
            if len(sources) == 0:
                self._add_source_callback(event)

        self.fig.canvas.draw()

    def _build_figure(self):
        buttons = [
            ('Add', '#88d353', '#77cc3a', self._add_source_callback),
            ('Remove', '#d35355', '#cc3a3c', self._remove_source_callback),
            ('Done', '#d2a550', '#cc9a3a', self._terminate),
            ]

        self.fig = plt.figure(0, clear=True)
        n_buttons = len(buttons)
        gs = mpl.gridspec.GridSpec(2, n_buttons,
            height_ratios=[1, .1], figure=self.fig)

        # background image
        self.im_ax = plt.subplot(gs[0, :])
        self.images_plot.plot_fov(self.im_ax, cmap='gray')

        # buttons
        self.button_axes = [plt.subplot(gs[1, i]) for i in range(n_buttons)]
        self.buttons = []
        for i, (name, c, hc, cb) in enumerate(buttons):
            button_ax = plt.subplot(gs[1, i])
            button = mpl.widgets.Button(
                button_ax, name,
                color=c, hovercolor=hc)
            button.on_clicked(cb)
            self.buttons.append(button)

        # sources
        self.sources = []
        x = np.mean(self.im_ax.get_xlim())
        y = np.mean(self.im_ax.get_ylim())
        for i in range(self.mandatory_sources):
            source = Source(
                self,
                self.circ_radius,
                annulus_radii=self.annulus_radii,
                color_active='#9941c3',
                color_inactive='#cda2e2')
            source.set_coordinates(x, y)
            self._activate_source(source)
            self.sources.append(source)

        self.callbacks_cid = []
        cid = self.fig.canvas.mpl_connect(
            'button_press_event', self._button_press_callback)
        self.callbacks_cid.append(cid)

        self.fig.suptitle(
            'Left click: move source — Right click: select source')

    def get_coordinates(self):
        return [src.get_coordinates() for src in self.sources]

    def run(self):
        self._build_figure()
        plt.show()
        self._clean()
        return self.get_coordinates()

    def _terminate(self, event):
        plt.close()

    def _clean(self):
        for cid in self.callbacks_cid:
            self.fig.canvas.mpl_disconnect(cid)

class SourcesMap():
    def __init__(self, images_plot, coordinates,
                 circ_radius=15, annulus_radii=None):
        self.images_plot = images_plot
        self.coordinates = coordinates
        self.circ_radius = circ_radius
        self.annulus_radii = annulus_radii

    def plot(self):
        self.fig = plt.figure(0, clear=True)

        # background image
        self.im_ax = plt.gca()
        self.images_plot.plot_fov(self.im_ax, cmap='gray')

        # sources
        self.sources = []
        x = np.mean(self.im_ax.get_xlim())
        y = np.mean(self.im_ax.get_ylim())
        for i, (x, y) in enumerate(self.coordinates):
            if i == 0:
                color = '#9941c3'
            else:
                color = '#6bc341'
            source = Source(self, self.circ_radius,
                            annulus_radii=self.annulus_radii,
                            color_active=color)
            source.set_coordinates(x, y)
            source.set_active()
            source.plot_label(short_source_label(i),
                              color=color, fontsize='small')
            self.sources.append(source)

    def savefig(self, *args, **kwargs):
        try:
            self.fig
        except AttributeError:
            self.plot()
        self.fig.savefig(*args, **kwargs)

def running_average_kernel(name, width):
    boxcar = lambda d: float(np.abs(d) <= width / 2)
    triangle = lambda d: (width / 2 - np.abs(d)) * boxcar(d)
    kernels = {
        'boxcar': boxcar,
        'triangle': triangle,
        'none': None,
        }
    return kernels[name]

@np.vectorize
def total_seconds(timedelta):
    ''' Vectorised version of timedelta.total_seconds() '''
    return timedelta.total_seconds()

def weighted_running_average(arr, weight_func, x=None):
    ''' Compute a running average on arr, weighting the contribution of each
    term with weight-function.

    Parameters
    ==========
    arr : np.ndarray(ndim=1)
        Array of values on which to compute the running average.
    weight_func : function
        A function which takes a distance as input, and returns a weight.
        Weights don't have to be normalised.
    x : np.ndarray or None (default: None)
        If x is an array, use it to compute the distances before they are
        passed to weight_func. This allows to compute a running average on
        non-regularly sampled data.

    Returns
    =======
    ret : np.ndarray
        An array of the same shape as the input arr, equivalent to:
        - when x is None:
            $ret_i = \sum_{j=-n}^n arr_{i+j} × w(j) / (2n+1)$
        - when x is specified:
            $ret_i = \sum_{j=-n}^n arr_{i+j} × w(|x_i - x_{i+j}|) / (2n+1)$.
    '''
    n = len(arr)
    if x is None:
        x = np.arange(n)

    distances = x.repeat(n).reshape(-1, n).T
    distances = np.abs(np.array([d - d[i] for i, d in enumerate(distances)]))
    weights = np.vectorize(weight_func)(distances)
    norm = n / weights.sum(axis=1).repeat(n).reshape(-1, n)
    weights *= norm

    ret = arr.copy()
    ret = arr.repeat(n).reshape(-1, n).T
    ret *= weights
    ret = np.mean(ret, axis=1)

    return ret


def main():
    args = get_parsed_args()

    # data reduction and image alignment --------------------------------------
    cube_path = os.path.join(args.output_path, 'cube_prepared.fits')
    aligned_cube_path = os.path.join(args.output_path, 'cube_aligned.fits')
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    ## This section is optimized to use only a little more memory than the
    ## size of `images`. It is equivalent to:
    #   images, timestamps = open_or_compute(
    #       cube_path reduction,
    #       args.sci_path, args.sci_dark_path,
    #       args.flat_path, args.flat_dark_path)
    #   images, timestamps = open_or_compute(
    #       aligned_cube_path,
    #       files_handling.pass_timestamps(alignment.register_stars),
    #       images, timestamps)
    if os.path.exists(aligned_cube_path):
        images, timestamps = files_handling.load_fits_data(
            aligned_cube_path, norm_to_exptime=False, timestamps_hdu=1)
    else:
        images, timestamps = files_handling.open_or_compute(
            cube_path, preparation.prepare,
            args.sci_path, args.sci_dark_path,
            args.flat_path, args.flat_dark_path,
            )
        images, timestamps = files_handling.open_or_compute(
            aligned_cube_path,
            files_handling.pass_timestamps(alignment.register_stars),
            images, timestamps)
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -


    # choose stars ------------------------------------------------------------
    coordinates_path = os.path.join(args.output_path, 'sources_coord.txt')
    images_plot = ImageSeriesPlot(images, default_image=1)
    circ_radius = 10
    annulus_radii = None

    if os.path.exists(coordinates_path):
        print("Loading coordinates from '{}'.".format(coordinates_path))
        coordinates = np.loadtxt(coordinates_path)
    else:
        sources_selector = MultipleSourcesSelector(
            images_plot, circ_radius=circ_radius, annulus_radii=annulus_radii,
            mandatory_sources=1)
        coordinates = sources_selector.run()
        print("Saving coordinates to '{}'.".format(coordinates_path))
        np.savetxt(coordinates_path, coordinates)

    sources_map = SourcesMap(
        images_plot, coordinates,
        circ_radius=circ_radius, annulus_radii=annulus_radii)
    sources_map.savefig(os.path.join(args.output_path, 'sources_map.pdf'))

    # photometry --------------------------------------------------------------
    photometry_path = os.path.join(args.output_path, 'sources_photometry.npz')
    if os.path.exists(photometry_path):
        print("Loading photometry from '{}'.".format(photometry_path))
        npz = np.load(photometry_path)
        all_sources = npz['all_sources']
        distances = npz['distances']
    else:

        n_stars = len(coordinates)
        n_images = len(images)
        distances = np.ndarray((n_images, n_stars))
        sample_sources = photometry.sep_extract(images[0], threshold=args.sep_threshold)
        all_sources = np.ndarray((n_images, n_stars), dtype=sample_sources.dtype)
        for i, im in enumerate(tqdm(images, desc='Photometry')):
            sources = photometry.sep_extract(im, threshold=args.sep_threshold)
            sources, dist = photometry.find_closest_sources(sources, coordinates)
            distances[i, :] = dist
            all_sources[i, :] = sources

        print("Saving photometry to '{}'.".format(photometry_path))
        np.savez(photometry_path, all_sources=all_sources, distances=distances)

    # plot data ---------------------------------------------------------------

    # intensities
    intensities = all_sources['flux']
    intensities /= np.nanmean(intensities, axis=0)
    intensity_exo = intensities[:, 0]
    intensity_ref = np.nanmean(intensities[:, 1:], axis=1)

    # relative intensities
    rel_intensity_exo = intensity_exo / intensity_ref
    rel_intensity_exo /= np.mean(rel_intensity_exo)
    # relative times
    t_abs = timestamps
    t_ref = t_abs[0]
    t_rel = total_seconds(t_abs - t_ref)
    t_unit = 's'

    # apply a running average to the data
    kernel = running_average_kernel(args.running_average,
                                    args.running_average_width)
    if kernel is not None:
        rel_intensity_exo = weighted_running_average(
            rel_intensity_exo, kernel, x=t_rel)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # plot --------------------------------------------------------------------

    # Plot the non-normalized light curve of each source.
    plt.figure(1)
    plt.clf()
    plt.xlim(t_abs.min(), t_abs.max())
    label_x_pos = plt.xlim()[1] 
    label_x_pos *= (1 + np.std(plt.xlim()) / np.mean(plt.xlim()) * .01)
    for i, source_intensity in enumerate(intensities.T):
        plt.axhline(i + 1, color='k', linewidth=.5, alpha=.5)
        plt.plot(t_abs, source_intensity + i,
                 linestyle='', marker='.', markersize=1.1)
        plt.text(label_x_pos, i + 1, source_label(i),
                 horizontalalignment='left', verticalalignment='center')
    fmt = mpl.dates.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(fmt)
    plt.xlabel('Time on {}'.format(t_ref.date()))
    plt.ylabel('Normalized intensity')
    plt.ylim(.5, len(intensities.T) + .5)
    plt.title('Light curves normalized to their averages', fontsize='small')
    plt.savefig(os.path.join(args.output_path, 'light_curves_normalized.pdf'))

    # Plot the light curve of each source, normalized to the average of the
    # light curves of all other sources.
    plt.figure(1)
    plt.clf()
    plt.xlim(t_abs.min(), t_abs.max())
    label_x = plt.xlim()[1] * (1 + np.std(plt.xlim()) / np.mean(plt.xlim()) * .01)
    for i, source_intensity in enumerate(intensities.T):
        plt.axhline(i + 1, color='k', linewidth=.5, alpha=.5)
        intensities_mask = np.ones(len(intensities.T))
        intensities_mask[0] = np.nan
        intensities_mask[i] = np.nan
        normalized_intensity = source_intensity / \
            np.nanmean(intensities * intensities_mask, axis=1)
        plt.plot(t_abs, normalized_intensity + i,
                 linestyle='', marker='.', markersize=1.1,
                 label=source_label(i))
        plt.text(
            label_x, i + 1,
            source_label(i),
            horizontalalignment='left',
            verticalalignment='center',
            )
    fmt = mpl.dates.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(fmt)
    plt.xlabel('Time on {}'.format(t_ref.date()))
    plt.ylabel('Relative intensity')
    plt.ylim(.5, len(intensities.T) + .5)
    plt.title('Light curves relative to other reference stars', fontsize='small')
    plt.savefig(os.path.join(args.output_path, 'light_curves_relative.pdf'))

    # plot exoplanet light curves
    plt.figure(1)
    plt.clf()
    plt.plot_date(t_abs, rel_intensity_exo, ',')
    fmt = mpl.dates.DateFormatter('%H:%M')
    plt.gca().xaxis.set_major_formatter(fmt)
    plt.xlabel('Time on {}'.format(t_ref.date()))
    plt.ylabel('Relative intensity')
    plt.ylim(.9, 1.1)
    plt.title(args.name)
    plt.savefig(os.path.join(args.output_path, 'light_curve_exoplanet.pdf'))

if __name__ == '__main__':
    main()
