from dataclasses import dataclass
from typing import Any

@dataclass(frozen=True)
class State:
    fluid: Any
    T: float
    P: float
    h: float
    s: float
    x: float | None = None
