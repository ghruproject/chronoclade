"""Exceptions shared across workflow modules."""


class WorkflowError(RuntimeError):
    """Raised when an external workflow stage fails."""
