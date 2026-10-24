# Python API

Every step, and every IDL routine it is built from, can be called from Python. The modules
are named after the steps; each function names, in its docstring, the Spextool routine it
reproduces.

## Organizing a night and running the steps

```{eval-rst}
.. automodule:: triplespec_reducer.night
   :members: organize, Night, Sequence, Frame, lookup_standard, safe_name

.. automodule:: triplespec_reducer.pipeline
   :members: Options, Reducer
```

## The steps

```{eval-rst}
.. automodule:: triplespec_reducer.frames
   :members:

.. automodule:: triplespec_reducer.flat
   :members:

.. automodule:: triplespec_reducer.wavecal
   :members:

.. automodule:: triplespec_reducer.extract
   :members:

.. automodule:: triplespec_reducer.combineimgs
   :members:

.. automodule:: triplespec_reducer.combine
   :members:

.. automodule:: triplespec_reducer.telluric
   :members:

.. automodule:: triplespec_reducer.merge
   :members:

.. automodule:: triplespec_reducer.clean
   :members:
```

## IDL arithmetic

```{eval-rst}
.. automodule:: triplespec_reducer.idl
   :members:
```
