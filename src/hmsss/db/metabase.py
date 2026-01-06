#!/usr/bin/python
import sqlite3
import os
import sys
import traceback
import time
from typing import List, Dict, Set, Any

from hmsss.core.logging import get_logger

logger = get_logger(__name__)

# TODO hier muss die database für die metadaten und reads der metagenome erstellt werden.
