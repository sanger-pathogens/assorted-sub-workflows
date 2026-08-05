#!/usr/bin/env python3

import sys
import joblib

model = joblib.load(sys.argv[1])
model_params = sys.argv[2:]

