"""
src/logger.py

Centralized logging setup for the Maillard Pipeline.
Replaces generic print() statements to provide timestamped,
severity-labeled, and reroutable output streams.
"""

import logging
import sys

# Define standard format
LOG_FORMAT = "%(asctime)s | %(name)-12s | %(levelname)-8s | %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """
    Returns a configured logger for the specified module name.
    """
    logger = logging.getLogger(name)
    
    # Only configure if no handlers are set to avoid duplicate logs in loops
    if not logger.handlers:
        logger.setLevel(level)
        
        # Configure console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)
        formatter = logging.Formatter(LOG_FORMAT, DATE_FORMAT)
        console_handler.setFormatter(formatter)
        
        logger.addHandler(console_handler)
        
        # Prevent propagation to root logger to avoid double-printing
        logger.propagate = False
        
    return logger
