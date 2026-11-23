@echo off
if "%~1"=="" (
  echo Usage: RUN_ALPHA.bat "C:\path\to\frozen_run"
  exit /b 2
)
python "%~dp0launch.py" --run-dir "%~1" --samples 128
exit /b %errorlevel%
