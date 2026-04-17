#!/bin/bash

set -e

echo "=================================="
echo "CSP Visualizer Setup"
echo "=================================="

cd "$(dirname "$0")"

# create venv if needed
if [ ! -d "csp_env" ]; then
    echo "Creating virtual environment..."
    python3 -m venv csp_env
fi

# ALWAYS use venv python explicitly (IMPORTANT FIX)
VENV_PYTHON="./csp_env/bin/python"
VENV_PIP="./csp_env/bin/pip"

# upgrade pip inside venv
echo "Upgrading pip..."
$VENV_PYTHON -m pip install --upgrade pip

# install dependencies INSIDE venv
echo "Installing dependencies..."
$VENV_PIP install pandas biopandas gemmi openpyxl

# run script using venv python (CRITICAL FIX)
echo "Running script..."
$VENV_PYTHON b-factor.py