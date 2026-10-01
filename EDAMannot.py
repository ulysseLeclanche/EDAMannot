from importlib import import_module

from edamannot import _EXPORTS as _PACKAGE_EXPORTS


__all__ = sorted(_PACKAGE_EXPORTS)


def __getattr__(name: str):
    """Lazily resolve a public EDAMannot API attribute.
    """
    try:
        module_name, attribute_name = _PACKAGE_EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(
            f"module 'EDAMannot' has no attribute {name!r}"
        ) from exc

    module = import_module(f"edamannot.{module_name}")
    return getattr(module, attribute_name)


def __dir__():
    """Return standard module attributes plus the public EDAMannot API."""
    return sorted(set(globals()) | set(__all__))
