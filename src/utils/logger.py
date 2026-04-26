"""
Logging configuration for the Primer Designer application.
Provides centralized logging with file and console output.
"""

import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Optional

# Global logger instance
_logger: Optional[logging.Logger] = None


def setup_logger(
    name: str = "primer_designer",
    log_dir: Optional[str] = None,
    level: int = logging.INFO,
    console_output: bool = True
) -> logging.Logger:
    """
    Set up and configure the application logger.

    Args:
        name: Logger name
        log_dir: Directory for log files. If None, uses default 'logs' directory.
        level: Logging level
        console_output: Whether to also output to console

    Returns:
        Configured logger instance
    """
    global _logger

    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Clear existing handlers
    logger.handlers.clear()

    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Set up log directory
    if log_dir is None:
        # Default to 'logs' directory relative to the package
        base_dir = Path(__file__).parent.parent.parent
        log_dir = base_dir / "logs"
    else:
        log_dir = Path(log_dir)

    log_dir.mkdir(parents=True, exist_ok=True)

    # Create file handler with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = log_dir / f"{name}_{timestamp}.log"
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Add console handler if requested
    if console_output:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    _logger = logger
    logger.info(f"Logger initialized. Log file: {log_file}")

    return logger


def get_logger() -> logging.Logger:
    """
    Get the application logger. Sets up a default logger if none exists.

    Returns:
        Logger instance
    """
    global _logger
    if _logger is None:
        _logger = setup_logger()
    return _logger


class LoggerMixin:
    """Mixin class that provides logging capabilities to any class."""

    @property
    def logger(self) -> logging.Logger:
        """Get the logger instance."""
        return get_logger()

    def log_info(self, message: str):
        """Log an info message."""
        self.logger.info(f"[{self.__class__.__name__}] {message}")

    def log_warning(self, message: str):
        """Log a warning message."""
        self.logger.warning(f"[{self.__class__.__name__}] {message}")

    def log_error(self, message: str):
        """Log an error message."""
        self.logger.error(f"[{self.__class__.__name__}] {message}")

    def log_debug(self, message: str):
        """Log a debug message."""
        self.logger.debug(f"[{self.__class__.__name__}] {message}")

    def log_exception(self, message: str):
        """Log an exception with traceback."""
        self.logger.exception(f"[{self.__class__.__name__}] {message}")
