# Copyright 2014-2016 The ODL development group
#
# This file is part of ODL.
#
# ODL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ODL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ODL.  If not, see <http://www.gnu.org/licenses/>.

"""Partial objects for per-iterate actions in iterative methods."""

# Imports for common Python 2/3 codebase
from __future__ import print_function, division, absolute_import
from future import standard_library
standard_library.install_aliases()

from abc import ABCMeta, abstractmethod
import time

from odl.util.utility import with_metaclass


__all__ = ('Partial', 'StorePartial', 'ForEachPartial', 'PrintTimingPartial',
           'PrintIterationPartial', 'PrintNormPartial', 'ShowPartial')


class Partial(with_metaclass(ABCMeta, object)):

    """Abstract base class for sending partial results of iterations."""

    @abstractmethod
    def __call__(self, result):
        """Apply the partial object to result.

        Parameters
        ----------
        result : `LinearSpaceVector`
            Partial result after n iterations

        Returns
        -------
        None
        """

    def __and__(self, other):
        """Return ``self & other``

        Compose partials, calls both in sequence.

        Parameters
        ----------
        other : `Partial`
            The other partial to compose with

        Returns
        -------
        result : `Partial`
            A partial whose `__call__` method calls both constituents
            partials.

        Examples
        --------
        >>> store = StorePartial()
        >>> iter = PrintIterationPartial()
        >>> both = store & iter
        >>> both
        StorePartial() & PrintIterationPartial()
        """
        return AndPartial(self, other)

    def __repr__(self):
        """Return ``repr(self)``"""
        return '{}()'.format(self.__class__.__name__)


class AndPartial(Partial):

    """Partial used for combining several partials"""

    def __init__(self, *callbacks):
        """Initialize an instance.

        Parameters
        ----------
        *callbacks : callable or `Partial`'s
            Partials to be called in sequence as listed.
        """
        partials = [c if isinstance(c, Partial) else ForEachPartial(c)
                    for c in callbacks]

        self.partials = partials

    def __call__(self, result):
        """Apply all partials to result."""
        for p in self.partials:
            p(result)

    def __repr__(self):
        """Return ``repr(self)``"""
        return ' & '.join('{}'.format(p) for p in self.partials)


class StorePartial(Partial):

    """Simple object for storing all partial results of the solvers.

    Can optionally apply a function, for example the norm or calculating the
    residual.
    """

    def __init__(self, results=None, function=None):
        """Initialize an instance.

        Parameters
        ----------
        results : list, optional
            List in which to store the partial results.
            Default: new list (``[]``)
        results : callable, optional
            Function to be called on all incomming results before storage.
            Default: copy

        Examples
        --------
        Store results as is

        >>> partial = StorePartial()

        Provide list to store partial results in.

        >>> results = []
        >>> partial = StorePartial(results=results)

        Store the norm of the results

        >>> norm_function = lambda x: x.norm()
        >>> partial = StorePartial(function=norm_function)
        """
        self._results = [] if results is None else results
        self._function = function

    @property
    def results(self):
        """The partial results."""
        return self._results

    def __call__(self, result):
        """Append result to results list."""
        if self._function:
            self._results.append(self._function(result))
        else:
            self._results.append(result.copy())

    def __iter__(self):
        """Allow iteration over the results"""
        return iter(self.results)

    def __getitem__(self, index):
        """Get partial result."""
        return self.results[index]

    def __len__(self):
        """The number of results stored."""
        return len(self.results)

    def __str__(self):
        """Return ``str(self)``"""
        resultstr = '' if self.results == [] else str(self.results)
        return 'StorePartial({})'.format(resultstr)

    def __repr__(self):
        """Return ``repr(self)``"""
        resultrepr = '' if self.results == [] else repr(self.results)
        return 'StorePartial({})'.format(resultrepr)


class ForEachPartial(Partial):

    """Simple object for applying a function to each iterate."""

    def __init__(self, function):
        """Initialize an instance.

        Parameters
        ----------
        function : callable
            Function to call for each iteration
        """
        assert callable(function)
        self.function = function

    def __call__(self, result):
        """Apply function to result."""
        self.function(result)

    def __str__(self):
        """Return ``str(self)``"""
        return 'ForEachPartial({})'.format(self.function)

    def __repr__(self):
        """Return ``repr(self)``"""
        return 'ForEachPartial({!r})'.format(self.function)


class PrintIterationPartial(Partial):

    """Print the iteration count."""

    _default_text = 'iter ='

    def __init__(self, text=None):
        """Initialize an instance.

        Parameters
        ----------
        text : string
            Text to display before the iteration count. Default: 'iter ='
        """
        self.text = text if text is not None else self._default_text
        self.iter = 0

    def __call__(self, _):
        """Print the current iteration."""
        print("{} {}".format(self.text, self.iter))
        self.iter += 1

    def __repr__(self):
        textstr = '' if self.text == self._default_text else self.text
        return 'PrintIterationPartial({})'.format(textstr)


class PrintTimingPartial(Partial):

    """Print the time elapsed since the previous iteration."""

    def __init__(self):
        """Initialize an instance."""
        self.time = time.time()

    def __call__(self, _):
        """Print time elapsed from the previous iteration."""
        t = time.time()
        print("Time elapsed = {:<5.03f} s".format(t - self.time))
        self.time = t


class PrintNormPartial(Partial):

    """Print the current norm."""

    def __init__(self):
        """Initialize an instance."""
        self.iter = 0

    def __call__(self, result):
        """Print the current norm."""
        print("norm = {}".format(result.norm()))
        self.iter += 1


class ShowPartial(Partial):

    """Show the partial result."""

    def __init__(self, *args, **kwargs):
        """Initialize a new instance.

        Parameters are passed through to the vectors show method. Additional
        parameters included.

        Parameters
        ----------
        display_step : positive int, optional
            Number of iterations between plots. Default: 1

        Other Parameters
        ----------------
        kwargs :
            Optional arguments passed on to ``x.show``
        """
        self.args = args
        self.kwargs = kwargs
        self.fig = kwargs.pop('fig', None)
        self.display_step = kwargs.pop('display_step', 1)
        self.iter = 0

    def __call__(self, x):
        """Show the current iterate."""
        if (self.iter % self.display_step) == 0:
            self.fig = x.show(*self.args, fig=self.fig, **self.kwargs)

        self.iter += 1

    def __repr__(self):
        """Return ``repr(self)``"""

        return '{}(display_step={}, fig={}, *{!r}, **{!r})'.format(
            self.__class__.__name__,
            self.display_step,
            self.fig,
            self.args,
            self.kwargs)


if __name__ == '__main__':
    # pylint: disable=wrong-import-position
    from odl.util.testutils import run_doctests
    run_doctests()
