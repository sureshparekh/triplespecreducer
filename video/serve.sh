#!/bin/bash
# Serve the films' page for video/render.py (Vite, http://localhost:5241).
# Uses a Vite from PATH, or set VITE to its executable.
cd "$(dirname "$0")/film"
exec "${VITE:-vite}" --config vite.config.mjs
