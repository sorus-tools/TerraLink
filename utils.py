from typing import Callable, Optional

from qgis.core import Qgis


def emit_progress(cb: Optional[Callable[[int, Optional[str]], None]], value: float, msg: Optional[str] = None) -> None:
    """Safe progress emitter that tolerates callback errors."""
    if cb is None:
        return
    try:
        cb(int(max(0, min(100, value))), msg)
    except Exception:
        pass


def log_error(iface, title: str, msg: str, level: int = Qgis.Critical) -> None:
    """Display an error via QGIS message bar if available, else print."""
    try:
        if iface and hasattr(iface, "messageBar"):
            iface.messageBar().pushMessage(title, msg, level=level)
        else:
            print(f"{title}: {msg}")
    except Exception:
        print(f"{title}: {msg}")
    raise RuntimeError(msg)
